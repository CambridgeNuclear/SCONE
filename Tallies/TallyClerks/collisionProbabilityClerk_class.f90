module collisionProbabilityClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile

  ! Basic tally modules
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill
  use tallyResult_class,          only : tallyResult

  ! Nuclear Data
  use nuclearDatabase_inter,      only : nuclearDatabase
  use neutronMaterial_inter,      only : neutronMaterial, neutronMaterial_CptrCast

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  implicit none
  private

  !!
  !! Collision probability matrix tally
  !!
  !! Uses collision estimator only
  !! Generates collision probability matrices for a generic multi-map.
  !! In principle this could be generalised to, e.g., scatter probability matrices.
  !! Doing this generalisation would require uncommenting a few lines within this
  !! file and allowing different responses to be fed to the clerk on initialisation.
  !!
  !! Notes:
  !!    -> The resulting collision probability matrix will account for regions outside
  !!       of the provided map, corresponding to the first index of the matrix.
  !!    -> This should be used with caution when there is a vacuum boundary: the particle
  !!       will never 'collide' there and so CP becomes a bit ambiguous.
  !!    -> If collision particle has invalid nuclear data type collision is ignored
  !!    -> CPM is stored in column-major order [prodBin, startBin].
  !!    -> CPs are only non-zero within an energy group. It may be more efficient to
  !!       define a slightly different CPClerk which has a separate map for energy.
  !!       With a fine energy and space discretisation, a large sparse matrix will be
  !!       produced by the current clerk which could be avoided with a slightly different
  !!       implementation.
  !!
  !! Private Members:
  !!   map      -> Map to divide phase-space into bins
  !!   resp     -> Response for transfer function
  !!               (Not presently used, would be macroTotal by default)
  !!   N        -> Number of Bins
  !!
  !! Interface:
  !!   tallyClerk Interface
  !!
  !! Sample dictionary input:
  !!
  !!  clerkName {
  !!      type collisionProbabilityClerk;
  !!      map { <TallyMapDef> }
  !!  }
  !!
  type, public, extends(tallyClerk) :: collisionProbabilityClerk
    private
    !! Map defining the discretisation
    class(tallyMap), allocatable           :: map
    integer(shortInt)                      :: N = 0 !! Number of bins
    !type(macroResponse)                    :: resp

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportInColl
    procedure  :: reportCycleEnd

    ! Overwrite default run-time result procedure
    procedure  :: getResult

    ! Output procedures
    procedure  :: display
    procedure  :: print

    ! Deconstructor
    procedure  :: kill

  end type collisionProbabilityClerk

  !!
  !! Collision probability matrix result class
  !!   Stored in column first order
  !!    dim1 -> target bin
  !!    dim2 -> orgin bin
  !!    dim3 -> 1 is values; 2 is STDs
  !!
  type,public, extends( tallyResult) :: CPMResult
    integer(shortInt)                           :: N  = 0 ! Size of CPM
    real(defReal), dimension(:,:,:),allocatable :: CPM    ! CPM proper
  end type CPMResult

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(collisionProbabilityClerk), intent(inout) :: self
    class(dictionary), intent(in)                   :: dict
    character(nameLen), intent(in)                  :: name

    ! Assign name
    call self % setName(name)

    ! Read maps
    call new_tallyMap(self % map, dict % getDictPtr('map'))

    ! Read size of the map
    ! Add 1 to allow for neutrons colliding in or originating from
    ! outside the map. This permits the CP tally to make sense when
    ! there is leakage or when the entire phase space is not covered
    self % N = self % map % bins(0) + 1

    ! Initialise response
    ! Not used at present - can be used if we desire to generalise
    ! collision probability estimators in future
    !call self % resp % build(macroTotal)

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(collisionProbabilityClerk),intent(in) :: self
    integer(shortInt),dimension(:),allocatable  :: validCodes

    validCodes = [inColl_CODE, cycleEnd_Code]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(collisionProbabilityClerk), intent(in) :: self
    integer(shortInt)                            :: S

    S = self % N * self % N

  end function getSize

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, xsData, mem, virtual)
    class(collisionProbabilityClerk), intent(inout) :: self
    class(particle), intent(in)                     :: p
    class(nuclearDatabase),intent(inout)            :: xsData
    type(scoreMemory), intent(inout)                :: mem
    logical(defBool), intent(in)                    :: virtual
    type(particleState)                             :: state
    integer(shortInt)                               :: sIdx, cIdx
    integer(longInt)                                :: addr
    real(defReal)                                   :: score
    class(neutronMaterial), pointer                 :: mat
    character(100), parameter :: Here = 'reportInColl (collisionProbabilityClerk_class.f90)'

    ! This clerk does not handle virtual scoring yet
    if (virtual) return

    ! Get material or return if it is not a neutron
    mat => neutronMaterial_CptrCast(xsData % getMaterial(p % matIdx()))

    if (.not.associated(mat)) return

    ! Find starting index in the map
    ! It is important that preCollision is not changed by a collisionProcessor
    ! before the particle is fed to the tally, otherwise results will be meaningless
    sIdx = self % map % map(p % preCollision)

    ! Find collision index in the map
    state = p
    cIdx = self % map % map(state)

    ! Invalid indices are allowed given that CPs must sum to one - this will include
    ! neutrons which collide outside the mapped region of phase space
    ! These correspond to index = 0

    ! Calculate collision probability
    ! Used the simple estimator - the commented line can allow CP to generalise to
    ! other responses
    ! For collision probability, top and bottom will cancel -- for other probabilities,
    ! this need not be the case
    !score = self % resp % get(p, xsData) * p % w / xsData % getTotalMatXS(p, p % matIdx())
    score = p % w

    ! I think this is right but I need to double check!
    ! Score element of the matrix
    addr = self % getMemAddress() + sIdx * self % N + cIdx
    call mem % score(score, addr)

  end subroutine reportInColl

  !!
  !! Process cycle end
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(collisionProbabilityClerk), intent(inout) :: self
    class(particleDungeon), intent(in)              :: end
    type(scoreMemory), intent(inout)                :: mem
    integer(shortInt)                               :: i, j
    integer(longInt)                                :: addrCPM, addrCPM0
    real(defReal)                                   :: val, normFactor

    if (mem % lastCycle()) then
      ! Set address to the start of collision probability matrix
      ! Decrease by 1 to get correct address on the first iteration of the loop
      addrCPM  = self % getMemAddress() - 1

      ! Normalise and accumulate estimates
      do i = 1, self % N

        ! Calculate normalisation factor
        ! by summing collision scores
        normFactor = ZERO
        addrCPM0 = addrCPM
        do j = 1, self % N
          addrCPM0 = addrCPM0 + 1
          val = mem % getScore(addrCPM0)
          normFactor = normFactor + val
        end do

        if (normFactor /= ZERO) normFactor = ONE / normFactor

        ! Normalise CPM column
        do j = 1, self % N
          addrCPM = addrCPM + 1
          call mem % closeBin(normFactor, addrCPM)
        end do

      end do
    end if

  end subroutine reportCycleEnd

  !!
  !! Return result from the clerk for interaction with Physics Package
  !! Returns CPMresult defined in this module
  !! If res is already allocated to a CPM of fitting size it reuses already allocated space
  !! This should improve performance when updating estimate of CPM each cycle
  !!
  !! See tallyClerk_inter for details
  !!
  pure subroutine getResult(self, res, mem)
    class(collisionProbabilityClerk), intent(in)  :: self
    class(tallyResult),allocatable, intent(inout) :: res
    type(scoreMemory), intent(in)                 :: mem
    integer(shortInt)                             :: i, j
    integer(longInt)                              :: addr
    real(defReal)                                 :: val, STD

    ! Allocate result to CPMresult
    ! Do not deallocate if already allocated to CPMresult
    ! Its not too nice -> clean up
    if (allocated(res)) then
      select type(res)
        class is (CPMresult)
          ! Do nothing

        class default ! Reallocate
          deallocate(res)
          allocate( CPMresult :: res)
     end select

    else
      allocate( CPMresult :: res)

    end if

    ! Load data into the CPM
    select type(res)
      class is(CPMresult)
        ! Check size and reallocate space if needed
        ! This is horrible. Hove no time to polish. Blame me (PMC this time)
        if (allocated(res % CPM)) then
          if( any(shape(res % CPM) /= [self % N, self % N, 2])) then
            deallocate(res % CPM)
            allocate(res % CPM(self % N, self % N, 2))
          end if
        else
          allocate(res % CPM(self % N, self % N, 2))
        end if

        ! Set size of the CPM
        res % N = self % N

        ! Load entries
        addr = self % getMemAddress() - 1
        do i = 1, self % N
          do j = 1, self % N
            addr = addr + 1
            call mem % getResult(val, STD, addr)
            res % CPM(j, i, 1) = val
            res % CPM(j, i, 2) = STD
          end do
        end do

    end select
  end subroutine getResult

  !!
  !! Display convergence progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(collisionProbabilityClerk), intent(in) :: self
    type(scoreMemory), intent(in)                :: mem

    print *, 'collisionProbabilityClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(collisionProbabilityClerk), intent(in) :: self
    class(outputFile), intent(inout)             :: outFile
    type(scoreMemory), intent(in)                :: mem
    integer(shortInt)                            :: i
    integer(longInt)                             :: addr
    real(defReal)                                :: val, std
    character(nameLen)                           :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Print map information
    call self % map % print(outFile)

    ! Print collision probability matrix
    name = 'CPM'
    addr = self % getMemAddress() - 1

    call outFile % startArray(name, [self % N, self % N])

    do i = 1, self % N * self % N
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
  !! See tallyClerk_inter for details
  !!
  elemental subroutine kill(self)
    class(collisionProbabilityClerk), intent(inout) :: self

    ! Call superclass
    call kill_super(self)

    if(allocated(self % map)) deallocate(self % map)
    self % N = 0
    !call self % resp % kill()

  end subroutine kill

end module collisionProbabilityClerk_class
