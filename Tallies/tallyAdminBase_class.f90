module tallyAdminBase_class

  use numPrecision
  use tallyCodes
  use genericProcedures,      only : fatalError, charCmp
  use dictionary_class,       only : dictionary
  use particle_class,         only : particle, phaseCoord
  use particleDungeon_class,  only : particleDungeon
  use tallyClerk_inter,       only : tallyClerk
  use tallyClerkSlot_class,   only : tallyClerkSlot
  use tallyClerkFactory_func, only : new_tallyClerk
  use outputFile_class,       only : outputFile


  implicit none
  private
  !! **** MOST LIKLEY CHANGE INTERFACES FOR CLERKS TO INCLUDE FLUX FOR IN COLLISION AND PATH
  !! **** PRECALCULATE FLUX HERE SO THERE IS NO NEED TO WARY ABOUT DYNAMIC TYPE OF XSDATA IN CLERKS
  !!
  !! Base class for the tallies black box.
  !! Its responsibilities are as flolow:
  !! 1) Accept events reports and routes then to individual tallyClerks
  !! 2) Returns k-eff estimate for a current cycle
  !! 3) Controls end of calculation
  !! 4) Controls printing of calculation progress (to a console)
  !! 5) Controls printing of result estimators to a file (filePath and file  Format)
  !!
  !! This class will be extended by inheritance to provide additional functionality
  !! i.e. return mesh based weight windows based on fission matrix or similar
  !!
  !!
  type, public:: tallyAdminBase
    private
    logical(defBool), public :: checkConvergence = .false.

    type(tallyClerkSlot),dimension(:),allocatable :: tallyClerks

    ! Lists of Clerks to be executed for each procedure
    integer(shortInt),dimension(:),allocatable    :: inCollClerks
    integer(shortInt),dimension(:),allocatable    :: outCollClerks
    integer(shortInt),dimension(:),allocatable    :: pathClerks
    integer(shortInt),dimension(:),allocatable    :: transClerks
    integer(shortInt),dimension(:),allocatable    :: histClerks
    integer(shortInt),dimension(:),allocatable    :: cycleStartClerks
    integer(shortInt),dimension(:),allocatable    :: cycleEndClerks

    ! List of clerks to display
    integer(shortInt), dimension(:),allocatable   :: displayList
    integer(shortInt), dimension(:),allocatable   :: triggerList

  contains
    ! Report Interface
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    ! Display procedures
    procedure :: display

    ! Convergance check
    procedure :: isConverged

    ! File writing procedures
    procedure :: print

    ! Build procedures
    procedure :: init
    procedure :: addTallyClerk
    procedure :: kill

    procedure,private :: addToReports

  end type tallyAdminBase


  public :: reportInColl
  public :: reportOutColl
  public :: reportPath
  public :: reportTrans
  public :: reportHist
  public :: reportCycleStart
  public :: reportCycleEnd
  public :: display
  public :: isConverged
  public :: init
  public :: print
  public :: kill

contains
  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self,p)
    class(tallyAdminBase), intent(inout) :: self
    class(particle), intent(in)          :: p
    integer(shortInt)                    :: i, idx

!    ! Go through all clerks that request the report
!    do i=1,size(self % inCollClerks)
!      idx = self % inCollClerks(i)
!      call self % tallyClerks(idx) % reportInColl(p)
!
!    end do

  end subroutine reportInColl

  !!
  !! Display convergance progress of selected tallies on the console
  !!
  subroutine display(self)
    class(tallyAdminBase), intent(in) :: self
    integer(shortInt)                 :: i

!    ! Go through all clerks marked as part of the display
!    do i=1,size(self % displayList)
!      call self % tallyClerks(i) % display()
!
!    end do

  end subroutine display

  !!
  !! Perform convergence check in selected clerks
  !!
  function isConverged(self) result(isIt)
    class(tallyAdminBase), intent(in)    :: self
    logical(defBool)                     :: isIt
    integer(shortInt)                    :: i,N

!    N = size( self % triggerList)
!
!    if( N > 0 ) then
!      isIt = self % tallyClerks(1) % isConverged()
!      do i = 2,N
!      isIt = isIt .and. self % tallyClerks(i) % isConverged()
!
!      end do
!    else
!      isIt = .false.
!
!    end if

  end function isConverged

  !!
  !! Add all results to outputfile
  !!
  subroutine print(self,output)
    class(tallyAdminBase), intent(in)    :: self
    class(outputFile), intent(inout)     :: output
    integer(shortInt)                    :: i

!    do i=1,size(self % tallyClerks)
!      call self % tallyClerks(i) % print(output)
!    end do

  end subroutine print

  !!
  !! Process outgoing collision report
  !! Assume that pre is AFTER any implicit treatment (i.e. implicit capture)
  !!
  subroutine reportOutColl(self,p,MT,muL)
    class(tallyAdminBase), intent(inout)  :: self
    class(particle), intent(in)           :: p
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    integer(shortInt)                     :: i, idx

!    ! Go through all clerks that request the report
!    do i=1,size(self % outCollClerks)
!      idx = self % outCollClerks(i)
!      call self % tallyClerks(idx) % reportOutColl(p,MT,muL)
!
!    end do

  end subroutine reportOutColl

  !!
  !! Process pathlength report
  !! ASSUMPTIONS:
  !! Pathlength must be contained within a single cell and material
  !!
  subroutine reportPath(self,p,L)
    class(tallyAdminBase), intent(inout) :: self
    class(particle), intent(in)          :: p
    real(defReal), intent(in)            :: L
    integer(shortInt)                    :: i, idx

!    ! Go through all clerks that request the report
!    do i=1,size(self % pathClerks)
!      idx = self % pathClerks(i)
!      call self % tallyClerks(idx) % reportPath(p,L)
!
!    end do

  end subroutine reportPath

  !!
  !! Process transition report
  !! ASSUMPTIONS:
  !! Transition must be a straight line
  !! Pre and Post direction is assumed the same (aligned with r_pre -> r_post vector)
  !!
  subroutine reportTrans(self,p)
    class(tallyAdminBase), intent(inout) :: self
    class(particle), intent(in)          :: p
    integer(shortInt)                    :: i, idx

!    ! Go through all clerks that request the report
!    do i=1,size(self % transClerks)
!      idx = self % transClerks(i)
!      call self % tallyClerks(idx) % reportTrans(p)
!
!    end do

  end subroutine reportTrans

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !!
  subroutine reportHist(self,p)
    class(tallyAdminBase), intent(inout) :: self
    class(particle), intent(in)          :: p
    integer(shortInt)                    :: i, idx

!    ! Go through all clerks that request the report
!    do i=1,size(self % histClerks)
!      idx = self % histClerks(i)
!      call self % tallyClerks(idx) % reportHist(p)
!
!    end do


  end subroutine reportHist

  !!
  !! Process beginning of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(tallyAdminBase), intent(inout) :: self
    class(particleDungeon), intent(in)   :: start
    integer(shortInt)                    :: i, idx

!    ! Go through all clerks that request the report
!    do i=1,size(self % cycleStartClerks)
!      idx = self % cycleStartClerks(i)
!      call self % tallyClerks(idx) % reportCycleStart(start)
!
!    end do

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(tallyAdminBase), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end
    integer(shortInt)                    :: i, idx

!    ! Go through all clerks that request the report
!    do i=1,size(self % cycleEndClerks)
!      idx = self % cycleEndClerks(i)
!      call self % tallyClerks(idx) % reportCycleEnd(end)
!
!    end do

  end subroutine reportCycleEnd

  !!
  !! Initialise tallyAdminBase form dictionary
  !!
  subroutine init(self,dict)
    class(tallyAdminBase), intent(inout)        :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen),dimension(:),allocatable :: clerks
    type(dictionary)                            :: locDict
    character(nameLen)                          :: entry
    logical(defBool)                            :: partOfDisplay, partOfTriggers
    integer(shortInt)                           :: i


!    ! Deallocate
!    call self % kill()
!
!    ! Allocate contents
!    allocate(self % inCollClerks(0)     )
!    allocate(self % outCollClerks(0)    )
!    allocate(self % pathClerks(0)       )
!    allocate(self % transClerks(0)      )
!    allocate(self % histClerks(0)       )
!    allocate(self % cycleStartClerks(0) )
!    allocate(self % cycleEndClerks(0)   )
!
!    allocate(self% displayList(0)       )
!    allocate(self% triggerList(0)       )
!
!    allocate(self % tallyClerks(0))
!
!    ! Read all dictionaries
!    call dict % keysDict(clerks)
!
!    ! Load all dictionaries as clerks
!    do i=1,size(clerks)
!      ! Copy dictionary to local copy
!      call dict % get(locDict,clerks(i))
!
!      ! Check if it is part of the display
!      call locDict % getOrDefault(entry,'display','no')
!      partOfDisplay = charCmp(entry,'yes')
!
!      ! Check if it is part of the convergance triggers
!      call locDict % getOrDefault(entry,'trigger','no')
!      partOfTriggers = charCmp(entry,'yes')
!
!      ! Get new clerk from factory and store it ina aslot
!      call self % addTallyClerk( new_tallyClerk(locDict,clerks(i)),partOfDisplay,partOfTriggers)
!
!    end do


  end subroutine init

  !!
  !! Attach new tally
  !!
  subroutine addTallyClerk(self, clerk, partOfDisplay, partOfTriggers)
    class(tallyAdminBase), intent(inout)          :: self
    class(tallyClerk), intent(in)                 :: clerk
    logical(defBool),intent(in)                   :: partOfDisplay
    logical(defBool),intent(in)                   :: partOfTriggers
    type(tallyClerkSlot)                          :: localSlot
    integer(shortInt),dimension(:),allocatable    :: reportCodes
    integer(shortInt)                             :: N, i

    character(100),parameter  :: Here = 'addTallyClerk (tallyAdminBase_class.f90)'

!    ! Check if provided clerk is a slot. Give error if it is
!    select type(clerk)
!      type is (tallyClerkSlot)
!        call fatalError(Here,'tallyCleakSlot was passed. It is forbidden to avoid nested slots.')
!    end select
!
!    ! Check if the tallyAdminBase is initialised
!    if( .not. allocated(self % tallyClerks) ) then
!      call fatalError(Here,'tallyAdminBase is uninitialised')
!    end if
!
!    ! Append tally Clerks. Automatic reallocation on assignment. F2008 feature
!    localSlot = clerk
!    self % tallyClerks = [self % tallyClerks, localSlot]
!
!    ! Obtain list of reports requested by the loaded clerk
!    N = size(self % tallyClerks)
!    reportCodes = self % tallyClerks(N) % validReports()
!
!    ! Append report sorting arrays with index of new tallyClerk
!    do i=1,size(reportCodes)
!      call self % addToReports( reportCodes(i), N )
!    end do
!
!    ! If clerk is partOfDisplay append displayList
!    if (partOfDisplay) self % displayList = [self % displayList, N ]
!
!    ! If clerk is partOfTriggers append triggerList
!    if (partOfTriggers) self % triggerList = [self % triggerList, N ]
!    if (partOfTriggers) self % checkConvergence = .true.

  end subroutine addTallyClerk

  !!
  !! Deallocates all content
  !!
  subroutine kill(self)
    class(tallyAdminBase), intent(inout) :: self

!    if(allocated(self % tallyClerks)) deallocate( self % tallyClerks )
!
!    if(allocated(self % displayList)) deallocate( self % displayList)
!    if(allocated(self % triggerList)) deallocate( self % triggerList)
!
!    if(allocated(self % inCollClerks))     deallocate( self % inCollClerks)
!    if(allocated(self % outCollClerks))    deallocate( self % outCollClerks )
!    if(allocated(self % pathClerks))       deallocate( self % pathClerks )
!    if(allocated(self % transClerks))      deallocate( self % transClerks )
!    if(allocated(self % histClerks))       deallocate( self % histClerks )
!    if(allocated(self % cycleStartClerks)) deallocate( self % cycleStartClerks )
!    if(allocated(self % cycleEndClerks))   deallocate( self % cycleEndClerks )


  end subroutine kill

  !!
  !! Append sorrting array identified with the code with tallyClerk idx
  !!
  subroutine addToReports(self,reportCode,idx)
    class(tallyAdminBase),intent(inout) :: self
    integer(shortInt), intent(in)       :: reportCode
    integer(shortInt), intent(in)       :: idx
    character(100),parameter  :: Here='addToReports (tallyAdminBase_class.f90)'

!    select case(reportCode)
!      case(inColl_CODE)
!        self % inCollClerks = [self % inCollClerks, idx]
!
!      case(outColl_CODE)
!        self % outCollClerks = [ self % outCollClerks, idx]
!
!      case(path_CODE)
!        self % pathClerks = [ self % pathClerks, idx]
!
!      case(trans_CODE)
!        self % transClerks = [ self % transClerks, idx]
!
!      case(hist_CODE)
!        self % histClerks = [ self % histClerks, idx]
!
!      case(cycleStart_CODE)
!        self % cycleStartClerks = [ self % cycleStartClerks, idx]
!
!      case(cycleEnd_CODE)
!        self % cycleEndClerks = [ self % cycleEndClerks, idx]
!
!      case default
!        call fatalError(Here, 'Undefined reportCode')
!    end select

  end subroutine addToReports


end module tallyAdminBase_class
