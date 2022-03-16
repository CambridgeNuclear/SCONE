module imcWeightClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use outputFile_class,           only : outputFile
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill

  ! Nuclear Data interface
  use nuclearDatabase_inter,      only : nuclearDatabase
  use materialMenu_mod,           only : mm_nMat => nMat

  ! Tally Filters
  use tallyFilter_inter,          only : tallyFilter
  use tallyFilterFactory_func,    only : new_tallyFilter

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap
  use materialMap_class,          only : materialMap

  ! Tally Responses
  use tallyResponseSlot_class,    only : tallyResponseSlot

  use tallyResult_class,     only : tallyResult, tallyResultEmpty
  use scoreMemory_class,     only : scoreMemory


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
  !! myImcWeightClerk {
  !!   type imcWeightClerk;
  !!   # filter { <tallyFilter definition> } #
  !!   # map    { <tallyMap definition>    } #
  !!   response (resName1 #resName2 ... #)
  !!   resName1 { <tallyResponse definition> }
  !!   #resNamew { <tallyResponse definition #
  !! }
  !!
  type, public, extends(tallyClerk) :: imcWeightClerk
    private
    ! Filter, Map & Vector of Responses
    class(tallyFilter), allocatable                  :: filter
    class(tallyMap), allocatable                     :: map
    type(tallyResponseSlot),dimension(:),allocatable :: response

    ! Useful data
    integer(shortInt)  :: width = 0
    integer(shortInt)  :: nMat

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! File reports and check status -> run-time procedures
    procedure  :: reportHist

    ! Output procedures
    procedure  :: display
    procedure  :: print
    procedure  :: getResult

  end type imcWeightClerk

  type,public, extends(tallyResult)           :: imcWeightResult
    real(defReal), dimension(:), allocatable  :: imcWeight
  end type imcWeightResult

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(imcWeightClerk), intent(inout)        :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    character(nameLen),dimension(:),allocatable :: responseNames
    integer(shortInt)                           :: i

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

    ! Store number of materials
    self % nMat = mm_nMat()

    ! Get names of response dictionaries
    !call dict % get(responseNames,'response')

    ! Load responses
    !allocate(self % response(size(responseNames)))
    !do i=1, size(responseNames)
    !  call self % response(i) % init(dict % getDictPtr( responseNames(i) ))
    !end do

    ! Set width
    !self % width = size(responseNames)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(imcWeightClerk), intent(inout) :: self

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
    class(imcWeightClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [hist_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(imcWeightClerk), intent(in) :: self
    integer(shortInt)                 :: S

    S = 1 !size(self % response)
    if(allocated(self % map)) S = S * self % map % bins(0)

  end function getSize

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportHist(self, p, xsData, mem)
    class(imcWeightClerk), intent(inout)  :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    type(particleState)                   :: state
    integer(shortInt)                     :: binIdx, i
    integer(longInt)                      :: adrr
    real(defReal)                         :: scoreVal, flx
    character(100), parameter :: Here =' reportHist (imcWeightClerk_class.f90)'

    ! Get current particle state
    state = p

    ! Check if within filter
    if(allocated( self % filter)) then
      if(self % filter % isFail(state)) return
    end if

    ! Find bin index
    if(allocated(self % map)) then
      binIdx = self % map % map(state)
    else
      binIdx = 1 !p % matIdx()
    end if

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Calculate bin address
    adrr = self % getMemAddress() + binIdx - 1

    ! Append to required bin
    if( p % isDead .and. p % fate /= LEAK_FATE ) then
      scoreVal = p % w
      call mem % score(scoreVal, adrr)
    end if

  end subroutine reportHist

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(imcWeightClerk), intent(in)  :: self
    type(scoreMemory), intent(in)      :: mem

    print *, 'imcWeightClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(imcWeightClerk), intent(in)          :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    real(defReal)                              :: val, std
    integer(shortInt)                          :: i
    integer(shortInt),dimension(:),allocatable :: resArrayShape
    character(nameLen)                         :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If imcWeight clerk has map print map information
    if( allocated(self % map)) then
      call self % map % print(outFile)
    end if

    ! Write results.
    ! Get shape of result array
    if(allocated(self % map)) then
      resArrayShape = [size(self % response), self % map % binArrayShape()]
    else
      resArrayShape = [size(self % response)]
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
    class(imcWeightClerk), intent(in)               :: self
    class(tallyResult), allocatable, intent(inout)  :: res
    type(scoreMemory), intent(in)                   :: mem
    real(defReal), dimension(:), allocatable        :: w
    integer(shortInt)                               :: i, N

    N = self % nMat
    allocate( w(N) )

    ! Get result value for each material
    do i = 1, N 
      call mem % getResult(w(i), self % getMemAddress()+i-1)
    end do

    allocate(res, source = imcWeightResult(w))

  end subroutine getResult

end module imcWeightClerk_class
