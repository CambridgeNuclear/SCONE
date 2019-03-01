module simpleFM_class

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

  ! Nuclear Data
  use transportNuclearData_inter, only : transportNuclearData

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
  !!      type simpleFM;
  !!      map { <TallyMapDef> }
  !!  }
  !!
  type, public, extends(tallyClerk) :: simpleFM
    private
    !! Map defining the discretisation
    class(tallyMap), allocatable :: map
    type(macroResponse)          :: resp
    integer(shortInt)            :: N = 0 !! Number of bins

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportInColl
    procedure  :: reportCycleEnd

    ! Output procedures
    procedure  :: display
    procedure  :: print

    ! Deconstructor
    procedure  :: kill
  end type simpleFM

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  subroutine init(self, dict, name)
    class(simpleFM), intent(inout) :: self
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name

    ! Read map
    call new_tallyMap(self % map, dict % getDictPtr('map'))

    ! Read size of the map
    self % N = self % map % bins(0)

    ! Initialise response
    call self % resp % build(macroNuFission)

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(simpleFM),intent(in)                 :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, cycleEnd_Code]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  elemental function getSize(self) result(S)
    class(simpleFM), intent(in) :: self
    integer(shortInt)           :: S

    S = self % N + self % N * self % N

  end function getSize

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self, p, mem)
    class(simpleFM), intent(inout)       :: self
    class(particle), intent(in)          :: p
    type(scoreMemory), intent(inout)     :: mem
    class(transportNuclearData), pointer :: xsDat
    type(particleState)                  :: state
    integer(shortInt)                    :: sIdx, cIdx
    integer(longInt)                     :: addr
    real(defReal)                        :: score

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

    ! Defend against invalid collision bin
    if(cIdx == 0) return

    ! Calculate fission neutron production
    score = self % resp % get(p)

    ! Score neutron production in the collision bin
    addr =  self % getMemAddress() + cIdx - 1
    call mem % score(score, addr)

    ! Defend againt invalid starting bin
    if(sIdx == 0) return

    ! Score element of the matrix
    addr = self % getMemAddress() + self % N + (sIdx - 1) * self % N + cIdx - 1
    call mem % score(score, addr)

  end subroutine reportInColl

  !!
  !! Process cycle end
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(simpleFM), intent(inout)     :: self
    class(particleDungeon), intent(in) :: end
    type(scoreMemory), intent(inout)   :: mem
    integer(shortInt)                  :: i, j
    integer(longInt)                   :: addrFM, addrPow
    real(defReal)                      :: normFactor

    if(mem % lastCycle()) then
      ! Set address to the start of Fission Matrix and Power vector
      ! Decrease by 1 to get correct addres on the fisrt iteration of the loop
      addrPow = self % getMemAddress() - 1
      addrFM  = addrPow + self % N -1

      ! Normalise and accumulate estimates
      do i=1,self % N
        ! Calculate power normalisation
        addrPow = addrPow + 1
        normFactor = mem % getScore(addrPow)
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
  !! Display convergance progress on the console
  !!
  subroutine display(self, mem)
    class(simpleFM), intent(in)   :: self
    type(scoreMemory), intent(in) :: mem

    print *, 'simpleFM does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  subroutine print(self, outFile, mem)
    class(simpleFM), intent(in)      :: self
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

    ! Print neutron production vector
    name = 'prodVector'
    addr = self % getMemAddress() - 1

    call outFile % startArray(name, [self % N])

    do i=1, self % N
      addr = addr + 1
      call mem % getResult(val, std, addr)
      call outFile % addResult(val, std)
    end do

    call outFile % endArray()

    ! Print fission matrix
    name = 'FM'

    call outFile % startArray(name, [self % N, self % N])

    do i=1,self % N * self % N
      addr = addr + 1
      call mem % getResult(val, std, addr)
      call outFile % addResult(val, std)
    end do

    call outFile % endBlock()

  end subroutine print

  !!
  !! Returns to uninitialised state
  !!
  elemental subroutine kill(self)
    class(simpleFM), intent(inout) :: self

    if(allocated(self % map)) deallocate(self % map)
    self % N = 0
    call self % resp % kill()

  end subroutine kill

end module simpleFM_class
