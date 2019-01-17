module collProbClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use outputFile_class,           only : outputFile

  use particle_class,             only : particle, phaseCoord, particleState
  use particleDungeon_class,      only : particleDungeon

  ! Tally Interfaces
  use tallyEstimator_class,       only : tallyEstimator
  use tallyClerk_inter,           only : tallyClerk
  use materialMap_class,          only : materialMap

  ! Nuclear Data

  implicit none
  private

  character(*),parameter :: CLASS_NAME = 'collProbClerk'

  !!
  !! Constructor
  !!
  interface collProbClerk
    module procedure collProbClerk_fromDict
  end interface

  !!
  !! Simple clerk to tally macroscopic reaction rates
  !! For now supports only collision estimator
  !!
  type, public,extends(tallyClerk) :: collProbClerk
    private
    character(nameLen)   :: name
    integer(shortInt)    :: cycleCount  = 0                ! Cycles counter
    real(defReal)        :: targetRelSD = 0.0

    ! Maping variables
    type(materialMap)    :: map

    ! Result variables
    integer(shortInt)                             :: N_map
    type(tallyEstimator),dimension(:),allocatable :: bins  ! Result estimates
    real(defReal),dimension(:),allocatable        :: out_w

  contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display
    procedure :: isConverged
    procedure :: init
    procedure :: print

    ! Overwrite report procedures
    procedure :: reportTrans
    procedure :: reportCycleEnd

  end type collProbClerk

contains

  !!
  !! Return codes for reports this clerk accepts
  !!
  function validReports(self) result(validCodes)
    class(collProbClerk), intent(in)              :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [trans_CODE, cycleEnd_CODE]

  end function validReports

  !!
  !! Display progress. Print total estimate and maximum relative STD s
  !!
  subroutine display(self)
    class(collProbClerk), intent(in)         :: self
    real(defReal)                            :: resSum, maxSTD

    print '(A)', 'collProbClerk does not support display'

  end subroutine display

  !!
  !! Perform convergance check in the Clerk
  !!
  function isConverged(self) result(isIt)
    class(collProbClerk), intent(in)         :: self
    logical(defBool)                         :: isIt
    real(defReal)                            :: maxSTD

    ! This functionality is not yet implemented
    isIt = .false.

  end function isConverged

  !!
  !! Initialise collProbClerk from dictionary
  !!
  subroutine init(self,dict,name)
    class(collProbClerk),intent(inout)          :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    integer(shortInt)                           :: N
    character(100),parameter :: Here ='init (collProbClerk_class.f90)'

    ! Assign name
    self % name = name

    ! Initialise maps
    self % map = materialMap(dict)

    ! Find size of map
    N = self % map % bins()
    self % N_map = N

    ! Create space for results
    allocate(self % bins(N*N))
    allocate(self % out_w(N))
    self % out_w = ZERO

  end subroutine init

  !!
  !! Write contents of the collProbClerk to output file
  !!
  subroutine print(self,outFile)
    class(collProbClerk), intent(in)   :: self
    class(outputFile), intent(inout)   :: outFile
    real(defReal)                      :: val, std
    integer(shortInt)                  :: i
    character(nameLen)                 :: name

    ! Begin block
    call outFile % startBlock(self % name)

    ! Write map data
    call self % map % print(outFile)

    ! Write axis descriptor
    name = 'axisDescriptor'
    call outFile % startArray(name,[2])
    name = 'StartingBin'
    call outFile % addValue(name)
    name = 'EndBin'
    call outFile % addValue(name)
    call outFile % endArray()

    ! Write results
    name = 'Res'
    call outFile % startArray(name,[self % N_map, self % N_map ])

    do i=1,size(self % bins)
      call self % bins(i) % getEstimate(val, std, self % cycleCount)
      call outFile % addResult(val,std)
    end do

    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

  !!
  !! Process transition report
  !!
  subroutine reportTrans(self,p)
    class(collProbClerk), intent(inout)  :: self
    class(particle), intent(in)          :: p
    type(particle)                       :: p_temp
    integer(shortInt)                    :: B_end, B_start, B
    real(defReal)                        :: w_end, w_start

    ! Find current bin
    B_end = self % map % map(p)

    ! Find starting bin *** BIT DIRTY
    p_temp  = p % preTransition
    B_start = self % map % map(p_temp)

    ! Return if tranbsistion is outside matrix
    if (B_end == 0 .or. B_start == 0 ) return

    ! Obtain starting and ending weights
    w_start = p % preTransition % wgt
    w_end   = p % w

    ! Calculate bin index
    B = self % N_map * (B_start - 1) + B_end

    ! Store result
    call self % bins(B) % add(w_end)

    ! Store outgoing weight
    self % out_w(B_start) = self % out_w(B_start) + w_start

  end subroutine reportTrans


  !!
  !! Process end of cycle report
  !!
  subroutine reportCycleEnd(self,end)
    class(collProbClerk), intent(inout)   :: self
    class(particleDungeon), intent(in)    :: end
    integer(shortInt)                     :: i, j, B

    self % cycleCount = self % cycleCount + 1

    ! Loop over start and end bins
    do i=1,self % N_map
      do j=1, self % N_map
        B = self % N_map * (i - 1) + j
        call self % bins(B) % closeBatch(self % out_w(i))
      end do
    end do

    self % out_w = ZERO

  end subroutine reportCycleEnd

  !!
  !! Return an instance of mactoClerk from dictionary and name
  !!
  function collProbClerk_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(collProbClerk)               :: new

    call new % init(dict,name)

  end function collProbClerk_fromDict

end module collProbClerk_class
