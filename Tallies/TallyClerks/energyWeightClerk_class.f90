module energyWeightClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState, P_PHOTON
  use outputFile_class,           only : outputFile
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill

  ! Nuclear Data interface
  use nuclearDatabase_inter,      only : nuclearDatabase

  ! Tally Filters
  use tallyFilter_inter,          only : tallyFilter
  use tallyFilterFactory_func,    only : new_tallyFilter

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Tally Responses
  use tallyResponseSlot_class,    only : tallyResponseSlot


  use tallyResult_class,     only : tallyResult, tallyResultEmpty

  implicit none
  private

  !!
  !! Colision estimator of reaction rates
  !! Calculates flux weighted integral from collisions
  !!
  !! Private Members:
  !!   filter   -> Space to store tally Filter
  !!   map      -> Space to store tally Map
  !!   response -> Array of responses
  !!   width    -> Number of responses (# of result bins for each map position)
  !!
  !! Interface
  !!   tallyClerk Interface
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myAbsorptionClerk {
  !!   type energyWeightClerk;
  !!   # filter { <tallyFilter definition> } #
  !!   # map    { <tallyMap definition>    } #
  !!   response (resName1 #resName2 ... #)
  !!   resName1 { <tallyResponse definition> }
  !!   #resNamew { <tallyResponse definition #
  !! }
  !!
  type, public, extends(tallyClerk) :: energyWeightClerk
    private
    ! Filter, Map & Vector of Responses
    class(tallyFilter), allocatable                  :: filter
    class(tallyMap), allocatable                     :: map
    type(tallyResponseSlot),dimension(:),allocatable :: response

    ! Useful data
    integer(shortInt)  :: width = 0

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportHist
    procedure  :: reportTrans

    ! Output procedures
    procedure  :: display
    procedure  :: print
    procedure  :: getResult

  end type energyWeightClerk


  type, public, extends(tallyResult)         :: energyWeightClerkResult
    real(defReal), dimension(:), allocatable :: materialEnergy
    real(defReal), dimension(:), allocatable :: radiationEnergy
  end type energyWeightClerkResult


  integer(shortInt), parameter :: MATERIAL_ENERGY  = 1
  integer(shortInt), parameter :: RADIATION_ENERGY = 2


contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(energyWeightClerk), intent(inout)     :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    character(nameLen),dimension(:),allocatable :: responseNames
    integer(shortInt)                           :: i
    type(dictionary)                            :: locDict
    character(100), parameter :: Here = 'init (energyWeightClerk_class.f90)'

    ! Assign name
    call self % setName(name)

    ! Load filter
    if( dict % isPresent('filter')) then
      call new_tallyFilter(self % filter, dict % getDictPtr('filter'))
    end if

    ! Load map
    if( dict % isPresent('map')) then
      call new_tallyMap(self % map, dict % getDictPtr('map'))
    end if

    ! Call error if responses are given
    if( dict % isPresent('response')) then
      call fatalError(Here, 'Warning: response not needed for energyWeightClerk')
    end if

    ! Initialise weight response automatically
    call locDict % init(1)
    call locDict % store('type','weightResponse')
    allocate(self % response(1))
    call self % response(1) % init(locDict)

    self % width = 2

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(energyWeightClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill and deallocate filter
    if(allocated(self % filter)) then
      deallocate(self % filter)
    end if

    ! Kill and deallocate map
    if(allocated(self % map)) then
      call self % map % kill()
      deallocate(self % map)
    end if

    ! Kill and deallocate responses
    if(allocated(self % response)) then
      deallocate(self % response)
    end if

    self % width = 0

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(energyWeightClerk),intent(in)        :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [trans_CODE,hist_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(energyWeightClerk), intent(in) :: self
    integer(shortInt)                 :: S

    S = self % width
    if(allocated(self % map)) S = S * self % map % bins(0)

  end function getSize




  subroutine reportHist(self, p, xsData, mem)
    class(energyWeightClerk), intent(inout)  :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    type(particleState)                   :: state
    integer(shortInt)                     :: binIdx, i
    integer(longInt)                      :: adrr
    real(defReal)                         :: scoreVal, flx
    character(100), parameter :: Here =' reportHist (energyWeightClerk_class.f90)'

    ! Get current particle state
    state = p

    if (p % fate == LEAK_FATE) return

    ! Check if within filter
    if(allocated( self % filter)) then
      if(self % filter % isFail(state)) return
    end if

    ! Find bin index
    if(allocated(self % map)) then
      binIdx = self % map % map(state)
    else
      binIdx = 1
    end if

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Calculate bin address
    adrr = self % getMemAddress() + self % width * (binIdx -1)  - 1

    ! Calculate flux sample 1/totXs
    flx = ONE / xsData % getTotalMatXS(p, p % matIdx())

    ! Append all bins
    scoreVal = self % response(1) % get(p, xsData) * p % w * flx
    ! Deal with infinite cross sections - may not be the right solution for generality
    if (scoreVal /= scoreVal) then
      scoreVal = p % w
    end if
    call mem % score(scoreVal, adrr + MATERIAL_ENERGY)

  end subroutine reportHist





  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportTrans(self, p, xsData, mem)
    class(energyWeightClerk), intent(inout)  :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    type(particleState)                   :: state
    integer(shortInt)                     :: binIdx, i
    integer(longInt)                      :: adrr
    real(defReal)                         :: scoreVal, flx
    character(100), parameter :: Here = 'reportTrans (energyWeightClerk_class.f90)'

    ! Get current particle state
    state = p

    ! Consider only radiation particles (those surviving timestep)
    if (p % fate /= AGED_FATE) return
    if (p % type /= P_PHOTON) return

    ! Check if within filter
    if(allocated( self % filter)) then
      if(self % filter % isFail(state)) return
    end if

    ! Find bin index
    if(allocated(self % map)) then
      binIdx = self % map % map(state)
    else
      binIdx = 1
    end if

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Calculate bin address
    adrr = self % getMemAddress() + self % width * (binIdx -1)  - 1

    ! Calculate flux sample 1/totXs
    flx = ONE / xsData % getTotalMatXS(p, p % matIdx())

    ! Append all bins
    scoreVal = self % response(1) % get(p, xsData) * p % w * flx
    ! Deal with infinite cross sections - may not be the right solution for generality
    if (scoreVal /= scoreVal) then
      scoreVal = p % w
    end if
    call mem % score(scoreVal, adrr + RADIATION_ENERGY)

  end subroutine reportTrans

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(energyWeightClerk), intent(in)  :: self
    type(scoreMemory), intent(in)      :: mem

    print *, 'energyWeightClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(energyWeightClerk), intent(in)       :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    real(defReal)                              :: val, std
    integer(shortInt)                          :: i
    integer(shortInt),dimension(:),allocatable :: resArrayShape
    character(nameLen)                         :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If collision clerk has map print map information
    if( allocated(self % map)) then
      call self % map % print(outFile)
    end if

    ! Write results.
    ! Get shape of result array
    if(allocated(self % map)) then
      resArrayShape = [self % width, self % map % binArrayShape()]
    else
      resArrayShape = [self % width]
    end if

    ! Start array
    name ='Res'
    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i=1,product(resArrayShape)
      call mem % getResult(val, std, self % getMemAddress() - 1 + i)
      call outFile % addResult(val,std)

    end do

    call outFile % endArray()
    call outFile % endBlock()

  end subroutine print

  !!
  !! Return result for interaction with Physics Package
  !! from the clerk in the slot
  !!
  !! See tallyClerk_inter for details
  !!
  pure subroutine getResult(self, res, mem)
    class(energyWeightClerk), intent(in)            :: self
    class(tallyResult), allocatable, intent(inout)  :: res
    type(scoreMemory), intent(in)                   :: mem
    real(defReal), dimension(:), allocatable        :: mat, rad
    integer(shortInt)                               :: i, N

    N = self % getSize() / 2
    allocate(mat(N))
    allocate(rad(N))

    ! Get result value for each material
    do i = 1, N 
      call mem % getResult(mat(i), self % getMemAddress()+2*i-2)
      call mem % getResult(rad(i), self % getMemAddress()+2*i-1)
    end do

    allocate(res, source = energyWeightClerkResult(mat,rad))

  end subroutine getResult

end module energyWeightClerk_class
