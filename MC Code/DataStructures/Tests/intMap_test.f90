module intMap_test
  use numPrecision
  use intMap_class, only : intMap
  use pFUnit_mod

  implicit none


@testCase
  type, extends(TestCase) :: test_intMap
    private
    type(intMap) :: map
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_intMap

  ! Parameters
  integer(shortInt), parameter :: VAL1 = 0
  integer(shortInt), parameter :: VAL2 = 1
  integer(shortInt), parameter :: VAL3 = -1
  integer(shortInt), parameter :: VAL4 =  huge(VAL1)
  integer(shortInt), parameter :: VAL5 = -huge(VAL1)
  integer(shortInt), parameter :: VAL6 =  133316666 ! Number of deamons of hell
  integer(shortInt), parameter :: VAL7 = -133316666

  integer(shortInt), parameter :: KEY1 = 0
  integer(shortInt), parameter :: KEY2 = 1
  integer(shortInt), parameter :: KEY3 = -1
  integer(shortInt), parameter :: KEY4 =  huge(VAL1)
  integer(shortInt), parameter :: KEY5 = -huge(VAL1)
  integer(shortInt), parameter :: KEY6 =  133316666 ! Number of deamons of hell
  integer(shortInt), parameter :: KEY7 = -133316666




contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp
    integer(shortInt)                 :: N, i

    ! Initialise
    call this % map % init(2)

    ! Load entries
     call this % map % add(KEY2, VAL2)
     call this % map % add(KEY3, VAL3)
     call this % map % add(KEY4, VAL4)
     call this % map % add(KEY5, VAL5)
     call this % map % add(KEY6, VAL6)
     call this % map % add(KEY7, VAL7)

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_intMap), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test retrieving values from the map.
  !! All values are present
  !!
@Test
  subroutine testRetrieval(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp

    temp = this % map % get(KEY2)
    @assertEqual(KEY2, temp)

    temp = this % map % get(KEY3)
    @assertEqual(KEY3, temp)

    temp = this % map % get(KEY4)
    @assertEqual(KEY4, temp)

    temp = this % map % get(KEY5)
    @assertEqual(KEY5, temp)

    temp = this % map % get(KEY6)
    @assertEqual(KEY6, temp)

    temp = this % map % get(KEY7)
    @assertEqual(KEY7, temp)

  end subroutine testRetrieval

@Test
  subroutine testOverwriting(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp

    temp = 7

    ! Store over existing keyword
    call this % map % add(KEY4, temp)

    ! Verify correctness
    @assertEqual(temp, this % map % get(KEY4))

  end subroutine testOverwriting


end module intMap_test
