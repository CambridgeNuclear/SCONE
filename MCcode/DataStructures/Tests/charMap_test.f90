module charMap_test
  use numPrecision
  use charMap_class, only : charMap
  use pFUnit_mod

  implicit none


@testCase
  type, extends(TestCase) :: test_charMap
    private
    type(charMap) :: map
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_charMap

  ! Parameters
  integer(shortInt), parameter :: VAL1 = 0
  integer(shortInt), parameter :: VAL2 = 1
  integer(shortInt), parameter :: VAL3 = -1
  integer(shortInt), parameter :: VAL4 =  huge(VAL1)
  integer(shortInt), parameter :: VAL5 = -huge(VAL1)
  integer(shortInt), parameter :: VAL6 =  133316666 ! Number of deamons of hell
  integer(shortInt), parameter :: VAL7 = -133316666

  character(nameLen), parameter :: KEY1 = 'Since you are reluctant '
  character(nameLen), parameter :: KEY2 = 'to provide us with the'
  character(nameLen), parameter :: KEY3 = 'location of the rebel base,'
  character(nameLen), parameter :: KEY4 = 'I have chosen to test this'
  character(nameLen), parameter :: KEY5 = 'station destructive power'
  character(nameLen), parameter :: KEY6 = 'on your home planet of'
  character(nameLen), parameter :: KEY7 = 'Alderaan'

contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_charMap), intent(inout) :: this
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
    class(test_charMap), intent(inout) :: this

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
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp

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
  !! Test Retrivial with deletions
  !!
@Test
  subroutine testRetrievalWithDel(this)
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp

    ! Delate Elements
    call this % map % del(KEY1)
    call this % map % del(KEY2)
    call this % map % del(KEY3)
    call this % map % del(KEY4)
    call this % map % del(KEY5)

    ! Obtain some elements
    temp = this % map % get(KEY6)
    @assertEqual(VAL6, temp)

    temp = this % map % get(KEY7)
    @assertEqual(VAL7, temp)

  end subroutine testRetrievalWithDel

  !!
  !! Test getting entery with default for absent keys
  !!
@Test
  subroutine testGetOrDefault(this)
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp
    integer(shortInt)                  :: default = 8
    character(nameLen)                 :: tempChar

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

    ! Get default.
    tempChar = 'this key is not present'
    @assertEqual(default, this % map % getOrDefault(tempChar, default))

  end subroutine testGetOrDefault

  !!
  !! Test deletions
  !!
@Test
  subroutine testDel(this)
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp

    ! Before deletion
    temp = this % map % get(KEY6)
    @assertEqual(VAL6, temp)

    call this % map % del(KEY6)

    ! After deletion
    temp = this % map % getOrDefault(KEY6, VAL1)
    @assertEqual(VAL1, temp)

  end subroutine testDel

  !!
  !! Test Looping over map
  !!
@Test
  subroutine testLooping(this)
    class(test_charMap), intent(inout)        :: this
    integer(shortInt),dimension(6),parameter  :: VALS = [VAL1, VAL2, VAL3, VAL4, VAL5, VAL7]
    character(nameLen),dimension(6),parameter :: KEYS = [KEY1, KEY2, KEY3, KEY4, KEY5, KEY7]
    character(nameLen),dimension(6)           :: KEYS_PAST
    integer(shortInt)                         :: counter, it, tVal
    character(nameLen)                        :: tKey

    ! Initialise parameters
    KEYS_PAST = "This is not a Key. It's picture of a key!"
    counter = 0

    ! Delete Entry 6 from map
    call this % map % del(KEY6)

    ! Loop over remaining elements
    it = this % map % begin()
    do while( it /= this % map % end() )
      tVal = this % map % atVal(it)
      tKey = this % map % atKey(it)

      @assertTrue(any(VALS == tVal),"Wrong Value")
      @assertTrue(any(KEYS == tKey),"Wrong Key")

      ! Make shure that KEY is not getting repeated
      @assertFalse(any(KEYS_PAST(1:counter) == tKey),"Repeated KEY")
      counter = counter + 1
      KEYS_PAST(counter) = tKey

      it = this % map % next(it)
    end do

    @assertEqual(6, counter)

  end subroutine testLooping

  !!
  !! Test Looping Edge Cases
  !!
@Test
  subroutine testLoopingEdgeCases(this)
    class(test_charMap), intent(inout) :: this
    type(charMap)                      :: locMap
    integer(shortInt)                  :: it

    ! Loop over uninitialised map
    it = locMap % begin()
    do while(it /= locMap % end())
      @assertTrue(.false.,"Should not enter the loop")
      it = locMap % next(it)
    end do

    ! Loop over empty map
    call locMap % init(8)
    it = locMap % begin()
    do while(it /= locMap % end())
      @assertTrue(.false.,"Should not enter the loop")
      it = locMap % next(it)
    end do

    ! Loop over map with deleted elements
    call locMap % add(KEY1,3)
    call locMap % add(KEY2,3)
    call locMap % add(KEY3,3)

    call locMap % del(KEY1)
    call locMap % del(KEY2)
    call locMap % del(KEY3)

    it = locMap % begin()
    !print *, it, locMap % atKey(it), locMap % map(:) % status
    do while(it /= locMap % end())
      @assertTrue(.false.,"Should not enter the loop")
      it = locMap % next(it)
    end do

  end subroutine testLoopingEdgeCases

  !!
  !! Test putting new entry under exoisting key
  !!
@Test
  subroutine testOverwriting(this)
    class(test_charMap), intent(inout) :: this
    integer(shortInt)                  :: temp

    temp = 7

    ! Store over existing keyword
    call this % map % add(KEY4, temp)

    ! Verify correctness
    @assertEqual(temp, this % map % get(KEY4))

  end subroutine testOverwriting

  !!
  !! Test getting length
  !!
@Test
  subroutine testGetLength(this)
    class(test_charMap), intent(inout) :: this

    @assertEqual(7, this % map % length())

  end subroutine testGetLength

  !!
  !! Test getOrDefaoult from uninitialised map
  !!
@Test
  subroutine testGetOrDefaultUninitialised(this)
    class(test_charMap), intent(inout) :: this
    character(nameLen)                 :: key
    type(charMap)                      :: locMap


    key = 'Key of admiral Dodanna'
    @assertEqual(7, locMap % getOrDefault(key, 7))

  end subroutine testGetOrDefaultUninitialised

end module charMap_test
