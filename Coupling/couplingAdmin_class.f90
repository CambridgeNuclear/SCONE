module couplingAdmin_class

  use numPrecision
  use universalVariables only : nameDensity, nameTemperature
  use dictionary_class,  only : dictionary
  use dictParser_func,   only : fileToDict
  use outputFile_class,  only : outputFile
  use errors_mod,        only : fatalError
  
  use geometryReg_mod,   only : gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr

  implicit none
  private

  character(nameLen), parameter, dimension(2) :: ALLOWABLE_FIELDS = ['temperature', 'density']
  character(nameLen), parameter :: END_SIGNAL = "SIGUSR2", &
                                   CONTINUE_SIGNAL = "SIGUSR1"

  !!
  !! Object to manage coupling between radiation and external physics solvers/driver scripts
  !!
  !! This object is responsible for managing the status of coupled calculations. It should do the following:
  !! - say whether coupling is happening (default: no)
  !! - know the means of communication by which coupling occurs (currently files are used as flags)
  !! - send a signal to the driver that it has finished a number of iterations
  !! - wait for the driver/other physics to finish their calculations
  !! - restart on receiving a signal from the driver
  !! - create any tallies necessary for coupling (such as power tallies)
  !! - output the results of those tallies to an appropriate location
  !! - know the location of fields provided by other physics and update them in SCONE
  !! 
  !! Private members:
  !! isCoupled       -> flag for whether coupling is being performed
  !! commFileIn      -> path to the file which triggers a SCONE calculation by the driver
  !! commFileOut     -> path to the file which SCONE outputs to the driver on completing a calculation
  !! outputFile      -> path to the file which SCONE outputs containing results to be used by the driver
  !! updateFreq      -> frequency with which SCONE pauses to exchange data with the driver
  !! maxIt           -> maximum number of iterations before finishing coupling
  !! tally           -> Tally admin containing necessary tallies for the driver
  !! fieldPaths      -> paths to files which SCONE reads to update the given fields
  !!
  !! Public interface:
  !! init          -> initialises the coupling setup from a dictionary
  !! kill          -> clean up
  !! couple        -> performs the main coupling logic
  !! endCoupling   -> deactivates coupling
  !!
  !! doCoupling    -> returns a logical if coupling is happening
  !! doUpdate      -> returns a logical if a physics update should be performed on a given iteration.
  !!
  !! makeTallies   -> returns a tally admin used for coupling
  !! outputTallies -> outputs tally data to be read by other physics
  !! updateFields  -> updates the relevant physics fields before resuming a calculation
  !!
  !! waitForSignal         -> pauses the calculation until a 'resume' signal is received
  !! signalIterationOver   -> signals that an iteration is complete
  !! signalCalculationOver -> signals that the coupled calculation is complete
  !!
  type, public :: couplingAdmin
    private
    logical(defBool)   :: isCoupled = .false.
    character(pathLen) :: commFileIn  = ''
    character(pathLen) :: commFileOut = ''
    character(pathLen) :: outputFile  = ''
    character(nameLen) :: outputFormat = ''
    integer(shortInt)  :: updateFreq = 0
    integer(shortInt)  :: maxIt = huge(shortInt)
    type(tallyAdmin)   :: tally
    character(pathLen), dimension(:), allocatable :: fieldPaths
  contains
    procedure :: init
    procedure :: kill
    procedure :: couple
    procedure :: endCoupling

    ! Status checks
    procedure :: doCoupling
    procedure, private :: doUpdate
    
    ! Results handling and manipulation
    procedure :: attachTally
    procedure :: outputTallies
    procedure, private :: updateFields

    ! Communication
    procedure, private :: waitForSignal
    procedure, private :: signalIterationOver
    procedure, private :: signalCalculationOver

  end type couplingAdmin

contains

  !!
  !! Initialises coupling from a dictionary
  !!
  subroutine init(self, dict)
    class(couplingAdmin), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    integer(shortInt)                   :: i, nFields
    logical(defBool)                    :: found, duplicate
    class(dictionary), pointer          :: tempDict
    type(outputFile)                    :: test_out
    character(100), parameter :: Here ='init (couplingAdmin_class.f90)'

    ! Read communication details
    call dict % get(self % commFileIn, 'receiveFile')
    call dict % get(self % commFileOut, 'sendFile')

    ! Read how often to exchange data
    call dict % get(self % updateFreq, 'updateFreq')
    if (self % updateFreq < 1) call fatalError(Here,&
            'Physics update frequency is silly: '//numToChar(self % updateFreq))

    ! Read maximum number of iterations over which to perform coupling
    call dict % getOrDefault(self % maxIt, 'maxIt', huge(shortInt))
    if (self % maxIt < 1) call fatalError(Here, &
            'Maximum iteration number is silly: '//numToChar(self % maxIt))

    ! Read tally info needed by other solvers
    tallyDict => dict % getDictPtr('tally')
    call self % tally % init(tallyDict)

    ! Read the output file name
    call dict % get(self % outputFile, 'outputFile')
    
    ! Initialise output file before calculation (so mistake in format will be caught early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Read the field update location(s)
    call dict % get(self % fieldPaths, 'fieldPaths')

    ! Read the field name(s)
    call dict % get(self % fieldNames, 'fieldNames')

    ! Ensure the paths correspond to the name
    if (size(self % fieldPaths) /= size(self % fieldNames)) call fatalError(Here,&
            'The number of field paths must match the number of field names.')

    ! Ensure the field names are unique and allowable
    nFields = size(self % fieldNames)
    do i = 1, nFields

      found = any(trim(self % fieldNames(i)) == trim(ALLOWABLE_FIELDS))
      if (.not. found) call fatalError(Here,'Field names must be one of: '//ALLOWABLE_FIELDS)

      if (i < nFields) then
        duplicate = any(trim(self % fieldNames(i)) == trim(self % fieldNames(i+1:nFields)))
        if (duplicate) call fatalError(Here,'Field names must be unique: '//trim(self % fieldNames(i)))
      end if
      
    end do
  
  end subroutine init

  !!
  !! Clean-up
  !!
  subroutine kill(self)
    class(couplingAdmin), intent(inout) :: self

    self % isCoupled = .false.
    self % updateFreq = 0
    self % commFileIn  = ''
    self % commFileOut = ''
    self % outputFile  = ''
    self % outputFormat  = ''
    self % maxIt = huge(shortInt)
    call self % tally % kill()
    if (allocated(fieldPaths)) deallocate(fieldPaths)

  end subroutine kill

  !!
  !! Subroutine which performs coupling logic during transport.
  !! Evaluates whether a physics update should be performed.
  !! If so, outputs tally information, cleans up tallies, signals and waits,
  !! and updates fields. Deactivates coupling if later than the maximum
  !! iteration.
  !!
  subroutine couple(self, it)
    class(couplingAdmin), intent(inout) :: self
    integer(shortInt), intent(in)       :: it

    if (.not. self % doCoupling()) return
    if (.not. self % doUpdate(it)) return

    call self % outputTallies()
    call self % tally % resetMemory()
    
    if (it < self % maxIt) then
      call self % signalIterationOver()
      call self % waitForSignal()
      call self % updateFields()
    else
      call self % signalCalculationOver()
      self % isCoupled = .false.
    end if

  end subroutine couple

  !!
  !! Returns whether coupling should be occurring or not
  !!
  function doCoupling(self) result(isCoupled)
    class(couplingAdmin), intent(in) :: self
    logical(defBool)                 :: isCoupled

    isCoupled = self % isCoupled
    
  end function doCoupling

  !!
  !! End coupling
  !!
  subroutine endCoupling(self)
    class(couplingAdmin), intent(inout) :: self

    if (self % isCoupled) then
      self % doCoupling = .false.
      call self % signalCalculationOver()
    end if

  end subroutine endCoupling

  !!
  !! Returns whether to perform a physics update based on iteration number
  !! and update frequency
  !!
  function doUpdate(self, it) result(update)
    class(couplingAdmin), intent(in) :: self
    integer(shortInt), intent(in)    :: it
    logical(defBool)                 :: update

    update = mod(it, self % updateFreq) == 0

  end function doUpdate

  !!
  !! Update the fields from the paths to the new field info
  !!
  subroutine updateFields(self)
    class(couplingAdmin), intent(in) :: self
    integer(shortInt)                :: i
    class(dictionary)                :: dict
    character(nameLen)               :: fieldName
    class(field), pointer            :: ptr
    character(100), parameter :: Here = 'updateFields (couplingAdmin_class.f90)'

    do i = 1, size(self % fieldPaths)
      call fileToDict(dict, self % fieldPaths(i))

      ! Obtain the pointer to the corresponding field.
      ! Reinitialise the field from the dictionary.
      select case(trim(self % fieldNames(i)))
        case('temperature')
          fieldName = nameTemperature

        case('density')
          fieldName = nameDensity

        case default
          call fatalError(Here, 'Unrecognised field type requested')
      end select

      ptr => gr_fieldPtr(gr_fieldIdx(fieldName))

      call ptr % kill()
      call ptr % init(dict)

    end do

  end subroutine updateFields

  !!
  !! Writes a file to signal that the current iteration is finished
  !!
  subroutine signalIterationOver(self)
    class(couplingAdmin), intent(in) :: self
    integer(shortInt)                :: unit

    open(newunit=unit, file=self % commFileOut, status="replace", action="write")
    write(unit, '(A)') CONTINUE_SIGNAL
    close(unit)

  end subroutine signalIterationOver

  !!
  !! Writes a file to signal that the entire calculation is finished
  !!
  subroutine signalCalculationOver(self)
    class(couplingAdmin), intent(in) :: self
    integer(shortInt)                :: unit

    open(newunit=unit, file=self % commFileOut, status="replace", action="write")
    write(unit, '(A)') END_SIGNAL
    close(unit)

  end subroutine signalCalculationOver

  !!
  !! Pauses the calculation while waiting for the existence
  !! of a signal file.
  !!
  subroutine waitForSignal(self)
    class(couplingAdmin), intent(in) :: self
    logical(defBool)                 :: exists
    integer(shortInt)                :: unit
    character(*)                     :: line

    do

      ! Check for the presence of a signal file
      inquire(file=self % commFileIn, exist=exists)

      if (exists) then
        ! Read and delete the communication file
        open(newunit=unit, file=self % commFileIn, status="old")

        read(unit, '(A)') line

        ! Check whether to continue iteration or end coupling
        if (trim(adjustl(line)) == END_SIGNAL) then
          self % isCoupled = .false.
        elseif (trim(adjustl(line)) /= CONTINUE_SIGNAL) then
          call fatalError(Here, 'Unrecognised signal: '//trim(adjustl(line)))
        end if

        close(unit, status="delete")
        exit
      end if
      
      call sleep(1)

    end do

  end subroutine waitForSignal

  !!
  !! Attach coupling tallies to another tally
  !!
  subroutine attachTally(self, tallyPtr)
    class(couplingAdmin), intent(in)          :: self
    class(tallyAdmin), pointer, intent(inout) :: tallyPtr
    
    if (self % doCoupling())
      call tallyPtr % push(self % tally)
    end if

  end subroutine attachTally

  !!
  !! Output tally results
  !!
  subroutine outputTallies(self)
    class(couplingAdmin), intent(in) :: self
    type(outputFile)                 :: out

    call out % init(self % outputFormat, filename = self % outputFile)
    call self % tally % print(out)

  end subroutine outputTallies

end module couplingAdmin_class
