module collisionClerk_class

  use numPrecision
  use tallyCodes
  use dictionary_class,        only : dictionary
  use particle_class,          only : particle, particleState
  use outputFile_class,        only : outputFile
  use scoreMemory_class,       only : scoreMemory
  use tallyClerk_inter,        only : tallyClerk

  ! Tally Filters
  use tallyFilter_inter,       only : tallyFilter
  use tallyFilterFactory_func, only : new_tallyFilter

  ! Tally Maps
  use tallyMap_inter,          only : tallyMap
  use tallyMapFactory_func,    only : new_tallyMap

  ! Tally Responses
  use tallyResponseSlot_class, only : tallyResponseSlot

  implicit none
  private

  !!
  !! Colision estimator of reaction rates
  !! Calculates flux weighted integral from collisions
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myCollisionClerk {
  !!   type collisionClerk;
  !!   # filter { <tallyFilter definition> } #
  !!   # map    { <tallyMap definition>    } #
  !!   response (resName1 #resName2 ... #}
  !!   resName1 { <tallyResponse definition> }
  !!   #resNamew { <tallyResponse definition #
  !! }
  !!
  type, public, extends(tallyClerk) :: collisionClerk
    private
    ! Filter, Map & Vector of Responses
    class(tallyFilter), allocatable                  :: filter
    class(tallyMap), allocatable                     :: map
    type(tallyResponseSlot),dimension(:),allocatable :: response

    ! Usefull data
    integer(shortInt)  :: width

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: validReports

    ! File reports and check status -> run-time procedures
    procedure  :: reportInColl

    ! Output procedures
    procedure  :: display
    procedure  :: print

  end type collisionClerk

contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  subroutine init(self, dict, name)
    class(collisionClerk), intent(inout)        :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen), intent(in)              :: name
    character(nameLen),dimension(:),allocatable :: responseNames
    integer(shortInt)                           :: i

    ! Assign name
    call self % setName(name)

    ! Load filetr
    if( dict % isPresent('filter')) then
      call new_tallyFilter(self % filter, dict % getDictPtr('filter'))
    end if

    ! Load map
    if( dict % isPresent('map')) then
      call new_tallyFilter(self % filter, dict % getDictPtr('filter'))
    end if

    ! Get names of response dictionaries
    call dict % get(responseNames,'response')

    ! Load responses
    allocate(self % response(size(responseNames)))
    do i=1, size(responseNames)
      call self % response(i) % init(dict % getDictPtr( responseNames(i) ))
    end do

  end subroutine init

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  function validReports(self) result(validCodes)
    class(collisionClerk),intent(in)           :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE]

  end function validReports

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self, p, mem)
    class(collisionClerk), intent(inout)  :: self
    class(particle), intent(in)           :: p
    type(scoreMemory), intent(inout)      :: mem
    type(particleState)                   :: state
    integer(shortInt)                     :: binIdx, i
    integer(longInt)                      :: adrr
    real(defReal)                         :: scoreVal

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
      binIdx = 0
    end if

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Calculate bin address
    adrr = self % getMemAddress() + self % width * (binIdx - 1) - 1

    do i=1,self % width
      ! shortInt will be replace with long soon
      scoreVal = self % response(i) % get(p) * p % w ! TODO: ADD 1/totXS
      call mem % score(scoreVal, int(adrr + i,shortInt))
    end do

  end subroutine reportInColl

  !!
  !! Display convergance progress on the console
  !!
  subroutine display(self, mem)
    class(collisionClerk), intent(in)  :: self
    type(scoreMemory), intent(inout)   :: mem

    print *, 'collisionClerk does not support display yet'

  end subroutine display

  !!
  !! Write contents of the clerk to output file
  !!
  subroutine print(self, outFile, mem)
    class(collisionClerk), intent(in) :: self
    class(outputFile), intent(inout)  :: outFile
    type(scoreMemory), intent(inout)  :: mem
    real(defReal)                     :: val, std
    integer(shortInt)                 :: i
    character(nameLen)                :: name

    ! Begin block
    call outFile % startBlock(self % getName())

    ! If collision clerk has map print map information
    if( allocated(self % map)) then
      call self % map % print(outFile)
    end if

    ! Write results. Associate resArrayShape with Result Array Shape Vector
    associate (resArrayShape => [size(self % response), self % map % binArrayShape()])
      name ='Res'
      call outFile % startArray(name, resArrayShape)

      ! Print results to the file
      do i=1,product(resArrayShape)
        ! TODO: Change interface of score memory to allow 64-bit adresses (longInt)
        call mem % getResult(val, std, int(self % getMemAddress() - 1 + i,shortInt))
        call outFile % addResult(val,std)

      end do

      call outFile % endArray()
      call outFile % endBlock()
    end associate

  end subroutine print

end module collisionClerk_class
