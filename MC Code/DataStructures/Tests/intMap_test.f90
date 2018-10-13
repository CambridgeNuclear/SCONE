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

    ! Load entries
     call this % map % add(KEY1, VAL1)
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

    call this % map % kill()

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

    temp = this % map % get(KEY1)
    @assertEqual(VAL1, temp)

    temp = this % map % get(KEY2)
    @assertEqual(VAL2, temp)

    temp = this % map % get(KEY3)
    @assertEqual(VAL3, temp)

    temp = this % map % get(KEY4)
    @assertEqual(VAL4, temp)

    temp = this % map % get(KEY5)
    @assertEqual(VAL5, temp)

    temp = this % map % get(KEY6)
    @assertEqual(VAL6, temp)

    temp = this % map % get(KEY7)
    @assertEqual(VAL7, temp)

  end subroutine testRetrieval

  !!
  !! Test getting entery with default for absent keys
  !!
@Test
  subroutine testGetOrDefault(this)
    class(test_intMap), intent(inout) :: this
    integer(shortInt)                 :: temp
    integer(shortInt)                 :: default = 8

    ! Google says its Chinese Lucky number. We do need luck
    default = 8

    ! Retrieve present entries. Do all cases to spot possible
    ! small differences between get and getOrDefault implementations
    temp = this % map % getOrDefault(KEY1, default)
    @assertEqual(VAL1, temp)

    temp = this % map % getOrDefault(KEY2, default)
    @assertEqual(VAL2, temp)

    temp = this % map % getOrDefault(KEY3, default)
    @assertEqual(VAL3, temp)

    temp = this % map % getOrDefault(KEY4, default)
    @assertEqual(VAL4, temp)

    temp = this % map % getOrDefault(KEY5, default)
    @assertEqual(VAL5, temp)

    temp = this % map % getOrDefault(KEY6, default)
    @assertEqual(VAL6, temp)

    temp = this % map % getOrDefault(KEY7, default)
    @assertEqual(VAL7, temp)

    ! Get default. (Yes I love Tolkien - MAK)
    @assertEqual(default, this % map % getOrDefault(1973,default))

  end subroutine testGetOrDefault

  !!
  !! Test putting new entry under exoisting key
  !!
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
