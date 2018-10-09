module macroClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use outputFile_class,           only : outputFile

  use particle_class,             only : particle, phaseCoord
  use particleDungeon_class,      only : particleDungeon

  ! Tally Interfaces
  use tallyEstimator_class,       only : tallyEstimator
  use tallyClerk_inter,           only : tallyClerk
  use tallyMap_inter,             only : tallyMap
  use tallyMapSlot_class,         only : tallyMapSlot
  use tallyMapFactory_func,       only : new_tallyMap

  ! Nuclear Data
  use transportNuclearData_inter, only : transportNuclearData
  use xsMacroSet_class,           only : xsMacroSet_ptr

  implicit none
  private

  character(*),parameter :: CLASS_NAME = 'macroClerk'

  interface macroClerk
    module procedure macroClerk_fromDict
  end interface

  !!
  !! Simple clerk to tally macroscopic reaction rates
  !! For now supports only collision estimator
  !!
  type, public,extends(tallyClerk) :: macroClerk
    private
    character(nameLen)   :: name
    integer(shortInt)    :: cycleCount  = 0                ! Cycles counter
    real(defReal)        :: targetRelSD = 0.0

    ! Maping variables
    logical(defBool)     :: hasMap =.false.
    type(tallyMapSlot)   :: map

    ! Result variables
    integer(shortInt)                             :: MT    ! Response MT. MT = 0 is flux.
    integer(shortInt)                             :: Nbins
    type(tallyEstimator),dimension(:),allocatable :: bins  ! Result estimates

  contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display
    procedure :: isConverged
    procedure :: init
    procedure :: print

    ! Overwrite report procedures
    procedure :: reportInColl
    procedure :: reportCycleEnd
  end type macroClerk

contains

  !!
  !! Return codes for reports this clerk accepts
  !!
  function validReports(self) result(validCodes)
    class(macroClerk), intent(in)              :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, cycleEnd_CODE]

  end function validReports

  !!
  !! Display progress. Print total estimate and maximum relative STD s
  !!
  subroutine display(self)
    class(macroClerk), intent(in)         :: self
    real(defReal),dimension(self % Nbins) :: res
    real(defReal),dimension(self % Nbins) :: std
    real(defReal)                         :: resSum, maxSTD

    ! Obtain current estimate of total rate and maximum std
    call self % bins % getEstimate(res, std, self % cycleCount)
    resSum = sum(res)

    ! Find maximum relative std
    where(res > FP_REL_TOL)
      std = std / res
    elsewhere
      std = ZERO
    end where
    maxSTD = maxval(std)

    ! Print estimates to a console
    print '(A,ES12.5,A,ES12.5)', 'Total rate: ', resSum, ' Maximum Rel. STD: ', maxSTD

  end subroutine display

  !!
  !! Perform convergance check in the Clerk
  !!
  function isConverged(self) result(isIt)
    class(macroClerk), intent(in)         :: self
    logical(defBool)                      :: isIt
    real(defReal),dimension(self % Nbins) :: res
    real(defReal),dimension(self % Nbins) :: std
    real(defReal)                         :: maxSTD

    ! Obtain current maximum STD
    call self % bins % getEstimate(res, std, self % cycleCount)

    ! Find maximum relative std
    where(res > FP_REL_TOL)
      std = std / res
    elsewhere
      std = ZERO
    end where
    maxSTD = maxval(std)

    ! Compare against target
    isIt = (maxSTD < self % targetRelSD)

  end function isConverged

  !!
  !! Initialise macroClerk from dictionary
  !!
  subroutine init(self,dict,name)
    class(macroClerk),intent(inout)             :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    character(nameLen)                          :: type
    character(nameLen),dimension(:),allocatable :: maps
    type(dictionary)                            :: tempDict
    character(100),parameter :: Here ='init (keffActiveClerk_class.f90)'

    ! Assign name
    self % name = name

    ! Verify if it is convergance trigger
    call dict % getOrDefault(type,'trigger','no')

    ! Read convergance target
    if( type == 'yes') then
      call dict % get(self % targetRelSD,'RelSDtarget')

    end if

    ! Read binning map
    call dict % keysDict(maps)

    if(size(maps) == 0) then
      self % hasMap = .false.

    else if(size(maps) == 1) then
      call dict % get(tempDict,maps(1))
      self % map  = new_tallyMap(tempDict)
      self % hasMap = .true.

    else
      call fatalError(Here,'Multiple maps are not yet supported')

    end if

    ! Load number of bins
    if( self % hasMap) then
      self % Nbins = self % map % bins()

    else
      self % Nbins = 1
    end if

    ! Allocate estimators
    allocate(self % bins(self % Nbins))

    ! Load response
    call dict % getOrDefault(self % MT,'MT',0)

  end subroutine init

  !!
  !! Write contents of the macroClerk to output file
  !!
  subroutine print(self,outFile)
    class(macroClerk), intent(in)    :: self
    class(outputFile), intent(inout) :: outFile
    real(defReal)                    :: val, std
    integer(shortInt)                :: i
    character(nameLen)               :: name

    ! Begin block
    call outFile % startBlock(self % name)

    if (self % hasMap) then
      ! Write axis data
      call self % map % print(outFile)

      ! Write axis descriptor
      name = 'axisDescriptor'
      call outFile % startArray(name,[1])
      call outFile % addValue(self % map % getAxisName())
      call outFile % endArray()
    end if

    ! Write results
    name = 'Res'
    call outFile % startArray(name,[self % Nbins])

    do i=1,self % Nbins
      call self % bins(i) % getEstimate(val, std, self % cycleCount)
      call outFile % addResult(val,std)
    end do

    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

  !!
  !! Process in collision report
  !!
  subroutine reportInColl(self,p)
    class(macroClerk), intent(inout) :: self
    class(particle), intent(in)      :: p
    type(xsMacroSet_ptr)             :: XSs
    real(defReal)                    :: totalXS, score
    integer(shortInt)                :: idx
    character(100), parameter  :: Here = 'reportInColl (macroClerk_class.f90)'

    ! Find index of the bin
    if (self % hasMap) then
      idx = self % map % map(p)
    else
      idx = 1
    end if

    if( idx == 0) return

    ! Obtain XSs
    ! Check if it dynamic type is supported
    ! If it is obtain macroscopic XSs
    ! It it isn't throw error
    associate (xsData => p % xsData)
      select type(xsData)
        class is (transportNuclearData)
          call xsData % getMatMacroXS(XSs, p, p % matIdx())

        class default
          call fatalError(Here,'Dynamic type of XS data attached to particle is not transportNuclearData')

      end select
    end associate

    totalXS  = XSs % totalXS()

    ! Calculate flux and scores
    score = p % w / totalXS
    if ( self % MT /= 0) score = score * XSs % xsOf(self % MT)
    ! Increment score
    call self % bins(idx) % add(score)

  end subroutine reportInColl

  !!
  !! Process end of cycle report
  !!
  subroutine reportCycleEnd(self,end)
    class(macroClerk), intent(inout)   :: self
    class(particleDungeon), intent(in) :: end

    self % cycleCount = self % cycleCount + 1
    call self % bins % closeBatch(1)

  end subroutine reportCycleEnd

  !!
  !! Return an instance of mactoClerk from dictionary and name
  !!
  function macroClerk_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(macroClerk)               :: new

    call new % init(dict,name)

  end function macroClerk_fromDict


end module macroClerk_class
