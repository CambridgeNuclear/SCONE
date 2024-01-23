module tallyAdmin_class

  use numPrecision
  use tallyCodes
  use genericProcedures,      only : fatalError, charCmp
  use dictionary_class,       only : dictionary
  use dynArray_class,         only : dynIntArray
  use charMap_class,          only : charMap
  use particle_class,         only : particle, particleState
  use particleDungeon_class,  only : particleDungeon
  use tallyClerk_inter,       only : tallyClerk
  use tallyClerkSlot_class,   only : tallyClerkSlot
  use tallyResult_class,      only : tallyResult, tallyResultEmpty
  use scoreMemory_class,      only : scoreMemory
  use outputFile_class,       only : outputFile

  ! Nuclear Data Interface
  use nuclearDataReg_mod,     only : ndReg_get => get
  use nuclearDatabase_inter,  only : nuclearDatabase

  implicit none
  private


  !! Parameters
  integer(longInt), parameter :: NO_NORM = -17_longInt

  !!
  !! TallyAdmin is responsible for:
  !!   1) Routing event reports to all tallyClerk that accept them
  !!   2) Performing normalisation of results
  !!   3) Printing selected results to screen to monitor calculation progress
  !!   4) Returning tallyResult objects from selected clerks for interaction with
  !!      Physics packages
  !!   5) Printing all results to a file at the end of calculation
  !!
  !! TallyAdmin can contain pointer to another tally admin (called attachment or atch).
  !! All reports send to tallyAdmin are also processed by an attachment. This allows to create
  !! linked lists of tallyAdmins. In predicted usage there should be at maximum 1 attachment
  !! related to data required by a physics package for a calculation.
  !!
  !! Private Members:
  !!   atch -> Pointer to an attachment tallyClerk (implements linked-list)
  !!   normBinAddr   -> Address of a bin used for normalisation
  !!   normValue     -> Target Value for normalisation
  !!   normClerkName -> Name of a Clerk used for normalisation
  !!   tallyClerks   -> Array of all defined tally Clerks
  !!   clerksNameMap -> CharMap that maps Clerk Name to its index in tallyClerks
  !!   inCollClerks     -> List of indices of all Clerks that require inCollReport
  !!   outCollClerks    -> List of indices of all Clerks that require outCollReport
  !!   pathClerks       -> List of indices of all Clerks that require pathReport
  !!   transClerks      -> List of indices of all Clerks that require transReport
  !!   spawnClerks      -> List of indices of all Clerks that require spawnReport
  !!   histClerks       -> List of indices of all Clerks that require histReport
  !!   cycleStartClerks -> List of indices of all Clerks that require cycleStartReport
  !!   cycleEndClerks   -> List of indices of all Clerks that require cycleEndReport
  !!   displayList      -> List of indices of all Clerks that are registered for display
  !!   mem              -> Score Memory for all defined Clerks
  !!
  !! Interface:
  !!   init   -> Initialise from dictionary
  !!   kill   -> Return to uninitialised state
  !!   push   -> Push new attachment to the end of the list
  !!   pop    -> Remove attachment from the end of the list
  !!   getEnd -> Returns pointer to the ned of the list (does not remove it like pop)
  !!   reportInColl     -> Process pre-collision reports in all clerks
  !!   reportOutColl    -> Process post-collision reports in all clerks
  !!   reportPath       -> Process pathlength reports in all clerks
  !!   reportTrans      -> Process transition reports in all clerks
  !!   reportSpawn      -> Process fission reports in all clerks
  !!   reportHist       -> Process History reports in all clerks
  !!   reportCycleStart -> Process Start Of Cycle reports in all clerks
  !!   reportCycleEnd   -> Process End of Cycle reports in all clerks
  !!   getResult        -> Return tallyResult object from a named Clerk
  !!   display     -> Call "display" on all Clerks registered to display
  !!   isConverged -> Return .true. if all convergance targets have been reached
  !!   print       -> Prints results to an output file object
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! mytallyAdmin {
  !!   #display (clerk2 clerk1); #
  !!   #norm    clerk3;          #       ! Clerk should be size 1 (first bin of clerk is normalised)
  !!   #normVal 13.0;            #       ! Must be present if "norm" is present
  !!   #batchSize   4;           #       ! Default value 1
  !!   clerk1 { <clerk definition here> }
  !!   clerk2 { <clerk definition here> }
  !!   clerk3 { <clerk definition here> }
  !!   clerk4 { <clerk definition here> }
  !! }
  !!
  type, public :: tallyAdmin
    private
    ! Attachment
    type(tallyAdmin),pointer :: atch => null()  ! Pointer to tallyAdmin attachment

    ! Normalisation data
    integer(longInt)   :: normBinAddr  = NO_NORM
    real(defReal)      :: normValue
    character(nameLen) :: normClerkName

    ! Clerks and clerks name map
    type(tallyClerkSlot),dimension(:),allocatable :: tallyClerks
    type(charMap)                                 :: clerksNameMap

    ! Lists of Clerks to be executed for each report
    type(dynIntArray)  :: inCollClerks
    type(dynIntArray)  :: outCollClerks
    type(dynIntArray)  :: pathClerks
    type(dynIntArray)  :: transClerks
    type(dynIntArray)  :: spawnClerks
    type(dynIntArray)  :: histClerks
    type(dynIntArray)  :: cycleStartClerks
    type(dynIntArray)  :: cycleEndClerks

    ! List of clerks to display
    type(dynIntArray)  :: displayList

    ! Score memory
    type(scoreMemory)  :: mem
  contains

    ! Build procedures
    procedure :: init
    procedure :: kill

    ! Attachment procedures
    procedure :: push   ! Add attachment to the end of the list
    procedure :: pop    ! Remove attachment from the end of the list
    procedure :: getEnd ! Copy pointer to the end of the list

    ! Report Interface
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportSpawn
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    ! Interaction procedures
    procedure :: getResult

    ! Display procedures
    procedure :: display

    ! Convergance check
    procedure :: isConverged

    ! File writing procedures
    procedure :: print

    procedure,private :: addToReports

  end type tallyAdmin

contains

  !!
  !! Initialise tallyAdmin form dictionary
  !!
  !! Args:
  !!   dict [in] -> Dictionary with definition
  !!
  !! Errors:
  !!   fatalError if there are mistakes in definition
  !!
  subroutine init(self,dict)
    class(tallyAdmin), intent(inout)            :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen),dimension(:),allocatable :: names
    integer(shortInt)                           :: i, j, cyclesPerBatch
    integer(longInt)                            :: memSize, memLoc
    character(100), parameter :: Here ='init (tallyAdmin_class.f90)'

    ! Clean itself
    call self % kill()

    ! Obtain clerks dictionary names
    call dict % keys(names,'dict')

    ! Allocate space for clerks
    allocate(self % tallyClerks(size(names)))

    ! Load clerks into slots and clerk names into map
    do i=1,size(names)
      call self % tallyClerks(i) % init(dict % getDictPtr(names(i)), names(i))
      call self % clerksNameMap % add(names(i),i)

    end do

    ! Register all clerks to recive their reports
    do i=1,size(self % tallyClerks)
      associate( reports => self % tallyClerks(i) % validReports() )
        do j=1,size(reports)
          call self % addToReports(reports(j), i)

        end do
      end associate
    end do

    ! Obtain names of clerks to display
    if( dict % isPresent('display')) then
      call dict % get(names,'display')

      ! Register all clerks to display
      do i=1,size(names)
        call self % displayList % add( self % clerksNameMap % get(names(i)))
      end do
    end if

    ! Read batching size
    call dict % getOrDefault(cyclesPerBatch,'batchSize',1)

    ! Initialise score memory
    ! Calculate required size.
    memSize = sum( self % tallyClerks % getSize() )
    call self % mem % init(memSize, 1, batchSize = cyclesPerBatch)

    ! Assign memory locations to the clerks
    memLoc = 1
    do i=1,size(self % tallyClerks)
      call self % tallyClerks(i) % setMemAddress(memLoc)
      memLoc = memLoc + self % tallyClerks(i) % getSize()

    end do

    ! Verify that final memLoc and memSize are consistant
    if(memLoc - 1 /= memSize) then
      call fatalError(Here, 'Memory addressing failed.')
    end if

    ! Read name of normalisation clerks if present
    if(dict % isPresent('norm')) then
      call dict % get(self % normClerkName,'norm')
      call dict % get(self % normValue,'normVal')
      i = self % clerksNameMap % get(self % normClerkName)
      self % normBinAddr = self % tallyClerks(i) % getMemAddress()
    end if

  end subroutine init

  !!
  !! Deallocates all content, returns to default values
  !!
  recursive subroutine kill(self)
    class(tallyAdmin), intent(inout) :: self

    ! Kill attchment
    if(associated(self % atch)) call self % atch % kill()

    ! Return parameters to default
    self % normBinAddr = NO_NORM
    self % atch => null()

    ! Kill clerks slots
    if(allocated(self % tallyClerks)) then
      call self % tallyClerks % kill()
      deallocate(self % tallyClerks)
    end if

    ! Kill processing lists
    call self % displayList % kill()

    call self % inCollClerks % kill()
    call self % outCollClerks % kill()
    call self % pathClerks % kill()
    call self % transClerks % kill()
    call self % spawnClerks % kill()
    call self % histClerks % kill()
    call self % cycleStartClerks % kill()
    call self % cycleEndClerks % kill()

    ! Kill score memory
    call self % mem % kill()

  end subroutine kill

  !!
  !! Add an attachment to the end of the list
  !!
  !! Args:
  !!   atch [in] -> Pointer to TallyAdmin to add at the end of the list
  !!
  !! Errors:
  !!   None
  !!   (There is some danger of infinate recursion if null initialisation of self % atch fails
  !!   for some reason. If that happens blame the compiler. Code is standard conforming)
  !!
  recursive subroutine push(self, atch)
    class(tallyAdmin), intent(inout)      :: self
    type(tallyAdmin), pointer, intent(in) :: atch

    if(associated(self % atch)) then
      call self % atch % push(atch)

    else
      self % atch => atch
    end if

  end subroutine push

  !!
  !! Remove attachment from the end of the list
  !!
  !! Note that any content pointed to by atch on entry will be lost!
  !!
  !! Args:
  !!   atch [out] -> Will be pointed to last element in linked list
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine pop(self, atch)
    class(tallyAdmin), intent(inout)       :: self
    type(tallyAdmin), pointer, intent(out) :: atch

    if(.not. associated(self % atch)) then ! Single element list
      atch => null()

    elseif( associated(self % atch % atch)) then ! Go down the list
      call self % atch % pop(atch)

    else ! Remove last element
      atch => self % atch
      self % atch => null()

    end if

  end subroutine pop

  !!
  !! Get pointer to the end element in the list
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Pointer to the attachment at the end of the linked-list
  !!
  !! Errors:
  !!   None
  !!
  recursive function getEnd(self) result(atch)
    class(tallyAdmin), intent(in)    :: self
    type(tallyAdmin),pointer         :: atch

    if(.not. associated(self % atch)) then
      atch => null()

    elseif( associated(self % atch % atch)) then
      atch => self % atch % getEnd()

    else
      atch => self % atch

    end if
  end function getEnd

  !!
  !! Display convergance progress of selected tallies on the console
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine display(self)
    class(tallyAdmin), intent(in) :: self
    integer(shortInt)             :: i
    integer(shortInt)             :: idx

    ! Call attachment
    if(associated(self % atch)) then
      call display(self % atch)
    end if

    ! Go through all clerks marked as part of the display
    do i=1,self % displayList % getSize()
      idx = self % displayList % get(i)
      call self % tallyClerks(idx) % display(self % mem)

    end do

  end subroutine display

  !!
  !! Perform convergence check in selected clerks
  !!
  !! NOT IMPLEMENTED YET
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   True is all convergance targets have been reached. False otherwise.
  !!
  !! Errors:
  !!   None
  !!
  function isConverged(self) result(isIt)
    class(tallyAdmin), intent(in)    :: self
    logical(defBool)                 :: isIt
    !integer(shortInt)                :: i,N

    isIt = .false.

  end function isConverged

  !!
  !! Add all results to outputfile
  !!
  !! Args:
  !!   output [inout] -> Output file object
  !!
  !! Errors:
  !!   None
  !!
  subroutine print(self,output)
    class(tallyAdmin), intent(in)        :: self
    class(outputFile), intent(inout)     :: output
    integer(shortInt)                    :: i
    character(nameLen)                   :: name

    ! Print tallyAdmin settings
    name = 'batchSize'
    call output % printValue(self % mem % getBatchSize(), name)

    ! Print Clerk results
    do i=1,size(self % tallyClerks)
      call self % tallyClerks(i) % print(output, self % mem)
    end do

  end subroutine print

  !!
  !! Process pre-collision report
  !!
  !! Assumptions:
  !!   Particle is provided just after transition. Before any implicit treatment.
  !!
  !! Args:
  !!   p [in]       -> Particle
  !!   virtual [in] -> Flag indicating virtual collision
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine reportInColl(self, p, virtual)
    class(tallyAdmin), intent(inout) :: self
    class(particle), intent(in)      :: p
    logical(defBool), intent(in)     :: virtual
    integer(shortInt)                :: i, idx
    class(nuclearDatabase),pointer   :: xsData
    character(100), parameter :: Here = "reportInColl (tallyAdmin_class.f90)"

    ! Call attachment
    if(associated(self % atch)) then
      call reportInColl(self % atch, p, virtual)
    end if

    ! Get Data
    xsData => ndReg_get(p % getType(), where = Here)

    ! Go through all clerks that request the report
    do i=1,self % inCollClerks % getSize()
      idx = self % inCollClerks % get(i)
      call self % tallyClerks(idx) % reportInColl(p, xsData, self % mem, virtual)

    end do

  end subroutine reportInColl

  !!
  !! Process post-collision report
  !!
  !! Assumptions:
  !!   PreCollision state in particle is set to just before this collision
  !!
  !! Args:
  !!   p [in]   -> Particle
  !!   MT [in]  -> MT number of the reaction particle has underwent
  !!   muL [in] -> Cosine of the collision deflection angle in LAB frame
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine reportOutColl(self, p, MT, muL)
    class(tallyAdmin), intent(inout)      :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    integer(shortInt)                     :: i, idx
    class(nuclearDatabase),pointer        :: xsData
    character(100), parameter :: Here = "reportOutColl (tallyAdmin_class.f90)"

    ! Call attachment
    if(associated(self % atch)) then
      call reportOutColl(self % atch, p, MT, muL)
    end if

    ! Get Data
    xsData => ndReg_get(p % getType(), where = Here)

    ! Go through all clerks that request the report
    do i=1,self % outCollClerks % getSize()
      idx = self % outCollClerks % get(i)
      call self % tallyClerks(idx) % reportOutColl(p, MT, muL, xsData, self % mem)

    end do

  end subroutine reportOutColl

  !!
  !! Process pathlength report
  !!
  !! Assumptions:
  !!   The entire path took place in a single material or unique cell
  !!   Report is called before particle moves to a new material or unique cell
  !!
  !! Args:
  !!   p [in] -> Particle
  !!   L [in] -> Length of the Path [cm]
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine reportPath(self, p, L)
    class(tallyAdmin), intent(inout)     :: self
    class(particle), intent(in)          :: p
    real(defReal), intent(in)            :: L
    integer(shortInt)                    :: i, idx
    class(nuclearDatabase),pointer       :: xsData
    character(100), parameter :: Here = "reportPath (tallyAdmin_class.f90)"

    ! Call attachment
    if(associated(self % atch)) then
      call reportPath(self % atch, p, L)
    end if

    ! Get Data
    xsData => ndReg_get(p % getType(), where = Here)

    ! Go through all clerks that request the report
    do i=1,self % pathClerks % getSize()
      idx = self % pathClerks % get(i)
      call self % tallyClerks(idx) % reportPath(p, L, xsData, self % mem)

    end do

  end subroutine reportPath

  !!
  !! Process transition report
  !!
  !! Assumptions:
  !!   TODO: Decide on the details of this report
  !!
  !! Args:
  !!   p [in] -> Particle
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine reportTrans(self, p)
    class(tallyAdmin), intent(inout) :: self
    class(particle), intent(in)      :: p
    integer(shortInt)                :: i, idx
    class(nuclearDatabase),pointer   :: xsData
    character(100), parameter :: Here = "reportTrans (tallyAdmin_class.f90)"

    ! Call attachment
    if(associated(self % atch)) then
      call reportTrans(self % atch, p)
    end if

    ! Get Data
    xsData => ndReg_get(p % getType(), where = Here)

    ! Go through all clerks that request the report
    do i=1,self % transClerks % getSize()
      idx = self % transClerks % get(i)
      call self % tallyClerks(idx) % reportTrans(p, xsData, self % mem)

    end do

  end subroutine reportTrans

  !!
  !! Process fission report
  !!
  !! Assumptions:
  !!   TODO: Decide on the details of this report
  !!
  !! Args:
  !!   pOld [in] -> Particle that caused the fission event
  !!   pNew [in] -> Particle state of the fission neutron
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine reportSpawn(self, pOld, pNew)
    class(tallyAdmin), intent(inout) :: self
    class(particle), intent(in)      :: pOld
    class(particleState), intent(in) :: pNew
    integer(shortInt)                :: i, idx
    class(nuclearDatabase),pointer   :: xsData
    character(100), parameter :: Here = "reportSpwan (tallyAdmin_class.f90)"

    ! Call attachment
    if(associated(self % atch)) then
      call reportSpawn(self % atch, pOld, pNew)
    end if

    ! Get Data
    xsData => ndReg_get(pOld % getType(), where = Here)

    ! Go through all clerks that request the report
    do i=1,self % spawnClerks % getSize()
      idx = self % spawnClerks % get(i)
      call self % tallyClerks(idx) % reportSpawn(pOld, pNew, xsData, self % mem)
    end do

  end subroutine reportSpawn

  !!
  !! Process history report
  !!
  !! Assumptions:
  !!   Particle has been killed and "fate" in the particle has been set.
  !!
  !! Args:
  !!   p [in] -> Particle
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine reportHist(self, p)
    class(tallyAdmin), intent(inout)  :: self
    class(particle), intent(in)       :: p
    integer(shortInt)                 :: i, idx
    class(nuclearDatabase),pointer    :: xsData
    character(100), parameter :: Here = "reportHist (tallyAdmin_class.f90)"

    ! Call attachment
    if(associated(self % atch)) then
      call reportHist(self % atch, p)
    end if

    ! Get Data
    xsData => ndReg_get(p % getType(), where = Here)

    ! Go through all clerks that request the report
    do i=1,self % histClerks % getSize()
      idx = self % histClerks % get(i)
      call self % tallyClerks(idx) % reportHist(p, xsData, self % mem)

    end do


  end subroutine reportHist

  !!
  !! Process beginning of a cycle
  !!
  !! Assumptions:
  !!   No particle have been transported yet in the cycle.
  !!   Report is called for the first time or after reportCycleEnd
  !!   All normalisations and modifications have already been applied to the particleDungeoin
  !!
  !! Args:
  !!   start [in] -> Particle Dungeon at the start of a cycle (after all normalisations)
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine reportCycleStart(self, start)
    class(tallyAdmin), intent(inout)   :: self
    class(particleDungeon), intent(in) :: start
    integer(shortInt)                  :: i
    integer(shortInt), save            :: idx
    !$omp threadprivate(idx)

    ! Call attachment
    if(associated(self % atch)) then
      call reportCycleStart(self % atch, start)
    end if

    ! Go through all clerks that request the report
    !$omp parallel do
    do i=1,self % cycleStartClerks % getSize()
      idx = self % cycleStartClerks % get(i)
      call self % tallyClerks(idx) % reportCycleStart(start, self % mem)
    end do
    !$omp end parallel do

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  !! Assumptions:
  !!   All particles given in "reportCycleStart" have been already transported
  !!   It is called after "reportCycleStart"
  !!   No modification or normalisation was applied to "end" particle Dungeon
  !!   "k_eff" member of end is set to criticality used to adjust fission source (implicit
  !!     fission site generation)
  !!
  !! Args:
  !!   end [in] -> Particle Dungeon at the end of a cycle (before any normalisations)
  !!
  !! Errors:
  !!   None
  !!
  recursive subroutine reportCycleEnd(self,end)
    class(tallyAdmin), intent(inout)   :: self
    class(particleDungeon), intent(in) :: end
    integer(shortInt)                  :: i
    integer(shortInt), save            :: idx
    real(defReal)                      :: normFactor, normScore
    character(100), parameter :: Here ='reportCycleEnd (tallyAdmin)class.f90)'
    !$omp threadprivate(idx)

    ! Call attachment
    if(associated(self % atch)) then
      call reportCycleEnd(self % atch, end)
    end if

    ! Go through all clerks that request the report
    !$omp parallel do
    do i=1,self % cycleEndClerks % getSize()
      idx = self % cycleEndClerks % get(i)
      call self % tallyClerks(idx) % reportCycleEnd(end, self % mem)
    end do
    !$omp end parallel do

    ! Calculate normalisation factor
    if( self % normBInAddr /= NO_NORM ) then
      normScore  = self % mem % getScore(self % normBinAddr)
      if (normScore == ZERO) then
        call fatalError(Here, 'Normalisation score from clerk:' // self % normClerkName // 'is 0')

      end if
      normFactor = self % normValue / normScore

    else
      normFactor = ONE
    end if

    ! Close cycle multipling all scores by multiplication factor
    call self % mem % closeCycle(normFactor)

  end subroutine reportCycleEnd

  !!
  !! Get result from the clerk defined by name
  !!
  !! Args:
  !!   res [inout] -> Allocatable Tally Result
  !!   name [in]   -> Name of the Clerk to get Result from
  !!
  !! Errors:
  !!   If "name" is not present in the admin res is allocated to tallyResultEmpty
  !!
  pure subroutine getResult(self, res, name)
    class(tallyAdmin), intent(in)                 :: self
    class(tallyResult),allocatable, intent(inout) :: res
    character(*), intent(in)                      :: name
    character(nameLen)                            :: name_loc
    integer(shortInt)                             :: idx
    integer(shortInt),parameter                   :: NOT_PRESENT = -3

    ! Deallocate if allocated result
    if(allocated(res)) deallocate(res)

    ! Copy name to character with nameLen
    name_loc = name

    ! Find clerk index
    idx = self % clerksNameMap % getOrDefault(name_loc, NOT_PRESENT)

    if(idx == NOT_PRESENT) then ! Return empty result
      allocate(res, source = tallyResultEmpty() )

    else ! Return result from the clerk named == name
      call self % tallyClerks(idx) % getResult(res, self % mem)

    end if

  end subroutine getResult

  !!
  !! Append sorting array identified with the code with tallyClerk idx
  !!
  !! Private helper subroutine
  !!
  !! Args:
  !!   reportCode [in] -> One of the report codes from tallyCodes
  !!   idx [in]        -> Clerk Index to add to approperiate List
  !!
  !! Errors:
  !!   fatalError if tallyCode is invalid
  !!
  subroutine addToReports(self, reportCode, idx)
    class(tallyAdmin),intent(inout) :: self
    integer(shortInt), intent(in)   :: reportCode
    integer(shortInt), intent(in)   :: idx
    character(100),parameter  :: Here='addToReports (tallyAdmin_class.f90)'

    select case(reportCode)
      case(inColl_CODE)
        call self % inCollClerks % add(idx)

      case(outColl_CODE)
        call self % outCollClerks % add(idx)

      case(path_CODE)
        call self % pathClerks % add(idx)

      case(trans_CODE)
        call self % transClerks % add(idx)

      case(spawn_CODE)
        call self % spawnClerks % add(idx)

      case(hist_CODE)
        call self % histClerks % add(idx)

      case(cycleStart_CODE)
        call self % cycleStartClerks % add(idx)

      case(cycleEnd_CODE)
        call self % cycleEndClerks % add(idx)

      case default
        call fatalError(Here, 'Undefined reportCode')
    end select

  end subroutine addToReports

end module tallyAdmin_class
