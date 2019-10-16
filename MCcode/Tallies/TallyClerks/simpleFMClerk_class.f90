module simpleFMClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile

  ! Basic tally modules
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk
  use tallyResult_class,          only : tallyResult

  ! Nuclear Data
  use nuclearDataReg_mod,         only : ndReg_get => get
  use nuclearDatabase_inter,      only : nuclearDatabase

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Tally Response
  use macroResponse_class,        only : macroResponse

  implicit none
  private

  !!
  !! Simple 1-D fission matrix
  !! This is a prototype implementation
  !! Uses collision estimator only
  !! Contains only a single map for discretisation
  !!
  !! Notes:
  !!    -> If collision particle has invalid nuclear data type collision is ignored
  !!    -> Collisions in non-fissile materials are ignored
  !!    -> FM is stored in column-major order [prodBin, startBin]
  !! Sample dictionary input:
  !!
  !!  clerkName {
  !!      type simpleFMClerk;
  !!      map { <TallyMapDef> }
  !!  }
  !!
  type, public, extends(tallyClerk) :: simpleFMClerk
    private
    !! Map defining the discretisation
    class(tallyMap), allocatable           :: map
    type(macroResponse)                    :: resp
    real(defReal),dimension(:),allocatable :: startWgt
    integer(shortInt)                      :: N = 0 !! Number of bins

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportCycleStart
    procedure  :: reportInColl
    procedure  :: reportCycleEnd

    ! Overwrite default run-time result procedure
    procedure  :: getResult

    ! Output procedures
    procedure  :: display
    procedure  :: print

    ! Deconstructor
    procedure  :: kill
  end type simpleFMClerk

  !!
  !! Fission matrix result class
  !!   Stored in column first order
  !!    dim1 -> target bin
  !!    dim2 -> orgin bin
  !!    dim3 -> 1 is values; 2 is STDs
  !!
  type,public, extends( tallyResult) :: FMresult
    integer(shortInt)                           :: N  = 0 ! Size of FM
    real(defReal), dimension(:,:,:),allocatable :: FM  ! FM proper
  end type FMResult

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  subroutine init(self, dict, name)
    class(simpleFMClerk), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    character(nameLen), intent(in)      :: name

    ! Assign name
    call self % setName(name)

    ! Read map
    call new_tallyMap(self % map, dict % getDictPtr('map'))

    ! Read size of the map
    self % N = self % map % bins(0)

    ! Allocate space for starting weights
    allocate(self % startWgt(self % N))

    ! Initialise response
    call self % resp % build(macroNuFission)

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(simpleFMClerk),intent(in)            :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, cycleStart_Code ,cycleEnd_Code]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  elemental function getSize(self) result(S)
    class(simpleFMClerk), intent(in) :: self
    integer(shortInt)                :: S

    S = self % N * self % N

  end function getSize

  !!
  !! Process start of the cycle
  !! Calculate starting weights in each bin
  !!
  subroutine reportCycleStart(self, start, mem)
    class(simpleFMClerk), intent(inout) :: self
    class(particleDungeon), intent(in)  :: start
    type(scoreMemory), intent(inout)    :: mem
    integer(shortInt)                   :: idx, i

    self % startWgt = ZERO

    ! Loop through a population and calculate starting weight in each bin
    do i=1,start % popSize()
      associate( state => start % get(i) )
        idx = self % map % map(state)
        if(idx > 0) self % startWgt(idx) = self % startWgt(idx) + state % wgt
      end associate
    end do

  end subroutine reportCycleStart

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self, p, mem)
    class(simpleFMClerk), intent(inout)  :: self
    class(particle), intent(in)          :: p
    type(scoreMemory), intent(inout)     :: mem
    class(transportNuclearData), pointer :: xsDat
    type(particleState)                  :: state
    integer(shortInt)                    :: sIdx, cIdx
    integer(longInt)                     :: addr
    real(defReal)                        :: score
    class(nuclearDatabase),pointer       :: xsData
    character(100), parameter :: Here = 'reportInColl simpleFMClear_class.f90'

    ! Get nuclear data or return if it is not transportNuclearData
    select type(xs => p % xsData)
      class is(transportNuclearData)
        xsDat => xs
      class default
        return
    end select

    ! Return if material is not fissile
    if(.not.xsDat % isFissileMat(p % matIdx())) return

    ! Find starting index in the map
    sIdx = self % map % map( p % preHistory)

    ! Find collision index in the map
    state = p
    cIdx = self % map % map(state)

    ! Defend against invalid collision or starting bin
    if(cIdx == 0 .or. sIdx == 0 ) return

    ! Calculate fission neutron production
    score = self % resp % get(p) * p % w / xsDat % getTotalMatXS(p, p % matIdx())

    ! Score element of the matrix
    addr = self % getMemAddress() + (sIdx - 1) * self % N + cIdx - 1
    call mem % score(score, addr)

  end subroutine reportInColl

  !!
  !! Process cycle end
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(simpleFMClerk), intent(inout) :: self
    class(particleDungeon), intent(in)  :: end
    type(scoreMemory), intent(inout)    :: mem
    integer(shortInt)                   :: i, j
    integer(longInt)                    :: addrFM
    real(defReal)                       :: normFactor

    if(mem % lastCycle()) then
      ! Set address to the start of Fission Matrix
      ! Decrease by 1 to get correct addres on the fisrt iteration of the loop
      addrFM  = self % getMemAddress() - 1

      ! Normalise and accumulate estimates
      do i=1,self % N
        ! Calculate normalisation factor
        normFactor = self % startWgt(i)
        if(normFactor /= ZERO) normFactor = ONE / normFactor

        do j=1,self % N
          ! Normalise FM column
          addrFM = addrFM + 1
          call mem % closeBin(normFactor, addrFM)
        end do
      end do
    end if

  end subroutine reportCycleEnd

  !!
  !! Return result from the clerk for interaction with Physics Package
  !!  Returns FMresult defined in this module
  !!   If res is already allocated to a FM of fitting size it reuses already allocated space
  !!    This should improve performance when updating estimate of FM each cycle
  !!
  pure subroutine getResult(self, res, mem)
    class(simpleFMClerk), intent(in)               :: self
    class(tallyResult),allocatable, intent(inout)  :: res
    type(scoreMemory), intent(in)                  :: mem
    integer(shortInt)                              :: i, j
    integer(longInt)                               :: addr
    real(defReal)                                  :: val, STD

    ! Allocate result to FMresult
    ! Do not deallocate if already allocated to FMresult
    ! Its not to nice -> clean up
    if(allocated(res)) then
      select type(res)
        class is (FMresult)
          ! Do nothing

        class default ! Reallocate
          deallocate(res)
          allocate( FMresult :: res)
     end select

    else
      allocate( FMresult :: res)

    end if

    ! Load data inti the FM
    select type(res)
      class is(FMresult)
        ! Check size and reallocate space if needed
        if( any(shape(res % FM) /= [self % N, self % N, 2])) then
          if(allocated(res % FM)) deallocate(res % FM)
          allocate(res % FM(self % N, self % N, 2))
        end if

        ! Set size of the FM
        res % N = self % N

        ! Load entries
        addr = self % getMemAddress() - 1
        do i = 1,self % N
          do j=1, self % N
            addr = addr + 1
            call mem % getResult(val, STD, addr)
            res % FM(j, i, 1) = val
            res % FM(j, i, 2) = STD
          end do
        end do

    end select
  end subroutine getResult

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self, mem)
    class(simpleFMClerk), intent(in) :: self
    type(scoreMemory), intent(in)    :: mem

    print *, 'simpleFMClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  subroutine print(self, outFile, mem)
    class(simpleFMClerk), intent(in) :: self
    class(outputFile), intent(inout) :: outFile
    type(scoreMemory), intent(in)    :: mem
    integer(shortInt)                :: i
    integer(longInt)                 :: addr
    real(defReal)                    :: val, std
    character(nameLen)               :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Print map information
    call self % map % print(outFile)

    ! Print fission matrix
    name = 'FM'
    addr = self % getMemAddress() - 1

    call outFile % startArray(name, [self % N, self % N])

    do i=1,self % N * self % N
      addr = addr + 1
      call mem % getResult(val, std, addr)
      call outFile % addResult(val, std)
    end do
    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

  !!
  !! Returns to uninitialised state
  !!
  elemental subroutine kill(self)
    class(simpleFMClerk), intent(inout) :: self

    if(allocated(self % map)) deallocate(self % map)
    if(allocated(self % startWgt)) deallocate(self % startWgt)
    self % N = 0
    call self % resp % kill()

  end subroutine kill

end module simpleFMClerk_class
