module collNumMap_test
  use numPrecision
  use funit
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use dictParser_func,         only : charToDict
  use outputFile_class,        only : outputFile

  use collNumMap_class,        only : collNumMap

  implicit none


@testCase
  type, extends(TestCase) :: test_collNumMap
    private
    type(collNumMap) :: map
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_collNumMap

  !!
  !! Test parameters
  !!
  integer(shortInt), dimension(*), parameter :: COLL_NUMS = [0, 1, 2, 5, 10, 50, 81]


contains

  !!
  !! Sets up test_collNumMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_collNumMap), intent(inout) :: this
    type(dictionary)                      :: dict

    ! Initialise dictionary and build map
    call dict % init(1)

    ! Build material map definition
    call dict % store('collNumbers', COLL_NUMS)
    call this % map % init(dict)

  end subroutine setUp

  !!
  !! Kills test_collNumMap object
  !!
  subroutine tearDown(this)
    class(test_collNumMap), intent(inout) :: this

    call this % map % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Mapping test
  !!
@Test
  subroutine testMapping(this)
    class(test_collNumMap), intent(inout)    :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(5)           :: bins
    integer(shortInt),dimension(5),parameter :: EXPECTED_BINS = [2, 3, 0, 0, 4]

    do i = 1,5
      state % collisionN = i
      bins(i) = this % map % map(state)
    end do

    @assertEqual(EXPECTED_BINS,bins)

  end subroutine testMapping

  !!
  !! Test number of bins inquiry
  !!
@Test
  subroutine testNumberOfBinsInquiry(this)
    class(test_collNumMap), intent(inout) :: this

    @assertEqual(7, this % map % bins(0), 'Total number of bins')
    @assertEqual(7, this % map % bins(1), 'Number of bins in dimension 1')
    @assertEqual(0, this % map % bins(2), 'Number of bins in higher dimension')

  end subroutine testNumberOfBinsInquiry

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequence is without errors
  !!
@Test
  subroutine testPrint(this)
    class(test_collNumMap), intent(inout) :: this
    type(outputFile)                      :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map % print(out)
    @assertTrue(out % isValid(),'For number of collisions map ')
    call out % reset()

  end subroutine testPrint


end module collNumMap_test
