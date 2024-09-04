module dictionary_test
  use numPrecision
  use dictionary_class, only : dictionary
  use pFUnit_mod
  implicit none

@TestCase
  type, extends(TestCase) :: test_dictionary
    type(dictionary) :: dict
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_dictionary

  !! Parameters
  real(defReal),parameter                  :: realVal     = 3.3_defReal
  integer(shortInt),parameter              :: boolVal     = 1
  integer(shortInt), parameter             :: intVal      = 1_shortInt
  character(nameLen),parameter             :: charNameLen = 'GoFortran_DownWithCpp'
  character(pathLen), parameter            :: charPathLen ='/home/KyloRen/VaderFanFic'
  real(defReal), dimension(2), parameter   :: realArray = [-1.0E-17_defReal, 14.7_defReal]
  integer(shortInt),dimension(3),parameter :: intArray =[-6475_shortInt, 13_shortInt, 14_shortInt]
  character(nameLen), dimension(1), parameter :: charNameLenArray = ['TK-421']
  character(pathLen), dimension(2), parameter :: charPathLenArray = ['C:\User\Tarkin\DeathStarPlans              ', &
                                                                     '/home/Dodonna/Articles/whyRebelsUseUNIX.odt']


contains

  !!
  !! Sets up test_dictionary object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_dictionary), intent(inout) :: this
    type(dictionary) :: tempDict

    call this % dict % init(1)
    call this % dict % store('myReal', realVal)
    call this % dict % store('myInt', intVal )
    call this % dict % store('myCharNameLen', charNameLen)
    call this % dict % store('myCharPathLen', charPathLen)
    call this % dict % store('realArray', realArray)
    call this % dict % store('intArray', intArray)
    call this % dict % store('charNameLenArray', charNameLenArray)
    call this % dict % store('charPathLenArray', charPathLenArray)
    call this % dict % store('myBool',boolVal)

    tempDict = this % dict
    call this % dict % store('nestedDict', tempDict)

  end subroutine setUp

  !!
  !! SKills test_dictionary object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_dictionary), intent(inout) :: this

    call this % dict % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!!  TESTS PROPER BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!!
!! Test extracting Real value
!!
@test
  subroutine testGettingReal(this)
    class(test_dictionary), intent(inout)    :: this
    real(defReal)                            :: tempReal

    call this % dict % get(tempReal,'myReal')
    @assertEqual(realVal, tempReal, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(tempReal,'myReal',7.0_defReal)
    @assertEqual(realVal, tempReal, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(tempReal,'invalid',7.0_defReal)
    @assertEqual(7.0_defReal, tempReal, 'Get or Default Retrival Failed for Absent Keyword')

    ! Int as real
    call this % dict % get(tempReal,'myInt')
    @assertEqual(real(intVal,defReal), tempReal, 'Retrival of int into real')

    call this % dict % getOrDefault(tempReal,'myInt',7.0_defReal)
    @assertEqual(real(intVal,defReal), tempReal, 'Get or Default Retrival of Int for Present Keyword')

  end subroutine testGettingReal

!!
!! Test extracting Real Array
!!
@test
  subroutine testGettingRealArray(this)
    class(test_dictionary), intent(inout)    :: this
    real(defReal),dimension(:),allocatable   :: tempReal
    real(defReal),dimension(:),pointer       :: tempReal_ptr => null()

    call this % dict % get(tempReal,'realArray')
    @assertEqual(realArray, tempReal, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(tempReal,'realArray', [7.0_defReal])
    @assertEqual(realArray, tempReal, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(tempReal,'myRealNon',[7.0_defReal])
    @assertEqual([7.0_defReal], tempReal, 'Get or Default Retrival Failed for Absent Keyword[ptr]')

    ! With pointer attribute
    call this % dict % get(tempReal_ptr,'realArray')
    @assertEqual(realArray, tempReal_ptr, 'Ordinary Retrival Failed[ptr]')

    call this % dict % getOrDefault(tempReal_ptr,'realArray', [7.0_defReal])
    @assertEqual(realArray, tempReal_ptr, 'Get or Default Retrival Failed for Present Keyword[ptr]')

    call this % dict % getOrDefault(tempReal_ptr,'myRealNon',[7.0_defReal])
    @assertEqual([7.0_defReal], tempReal_ptr, 'Get or Default Retrival Failed for Absent Keyword[ptr]')

    ! Retrival of int array into real array
    call this % dict % get(tempReal,'intArray')
    @assertEqual(real(intArray,defReal), tempReal, 'Ordinary Retrival of int Failed')

    call this % dict % getOrDefault(tempReal,'intArray', [7.0_defReal])
    @assertEqual(real(intArray,defReal), tempReal, 'Get or Default int Retrival for Present Keyword')

    call this % dict % get(tempReal_ptr,'intArray')
    @assertEqual(real(intArray,defReal), tempReal_ptr, 'Ordinary Retrival of int Failed')

    call this % dict % getOrDefault(tempReal_ptr,'intArray', [7.0_defReal])
    @assertEqual(real(intArray,defReal), tempReal_ptr, 'Get or Default int Retrival for Present Keyword')

  end subroutine testGettingRealArray

!!
!! Test extracting Integer Value
!!
@test
  subroutine testGettingInt(this)
    class(test_dictionary), intent(inout)    :: this
    integer(shortInt)                        :: temp

    call this % dict % get(temp,'myInt')
    @assertEqual(intVal, temp, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'myInt',7_shortInt)
    @assertEqual(intVal, temp, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', 7_shortInt)
    @assertEqual(7_shortInt, temp, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingInt


!!
!! Test extracting logical Value
!!
@test
  subroutine testGettingBool(this)
    class(test_dictionary), intent(inout)    :: this
    logical(defBool)                         :: temp

    call this % dict % get(temp,'myBool')
    @assertTrue(temp, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'myBool',.false.)
    @assertTrue( temp, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', .false.)
    @assertFalse(temp, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingBool

!!
!! Test extracting Integer Array
!!
@test
  subroutine testGettingIntArray(this)
    class(test_dictionary), intent(inout)      :: this
    integer(shortInt),dimension(:),allocatable :: temp
    integer(shortInt),dimension(:),pointer     :: temp_ptr => null()

    call this % dict % get(temp,'intArray')
    @assertEqual(intArray, temp, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'intArray',[7_shortInt])
    @assertEqual(intArray, temp, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', [7_shortInt])
    @assertEqual([7_shortInt], temp, 'Get or Default Retrival Failed for Absent Keyword')

    ! With pointer attribute
    call this % dict % get(temp_ptr,'intArray')
    @assertEqual(intArray, temp_ptr, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp_ptr,'intArray',[7_shortInt])
    @assertEqual(intArray, temp_ptr, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp_ptr,'invalid', [7_shortInt])
    @assertEqual([7_shortInt], temp_ptr, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingIntArray

!!
!! Test extracting nameLen long character
!!
@test
  subroutine testGettingNameLenChar(this)
    class(test_dictionary), intent(inout) :: this
    character(nameLen)                    :: temp
    character(nameLen),parameter :: default = 'Mes Que Nada'

    call this % dict % get(temp,'myCharNameLen')
    @assertEqual(charNameLen, temp, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'myCharNameLen',default)
    @assertEqual(charNameLen, temp, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', default)
    @assertEqual(default , temp, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingNameLenChar

!!
!! Test extracting pathLen long character
!!
@test
  subroutine testGettingPathLenChar(this)
    class(test_dictionary), intent(inout) :: this
    character(pathLen)                    :: temp
    character(pathLen),parameter  :: default = 'Mes Que Nada'

    call this % dict % get(temp,'myCharPathLen')
    @assertEqual(charPathLen, temp, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'myCharPathLen', default)
    @assertEqual(charPathLen, temp, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', default)
    @assertEqual(default , temp, 'Get or Default Retrival Failed for Absent Keyword')

  end subroutine testGettingPathLenChar

!!
!! Test extracting nameLen long array of Chars
!!
@test
  subroutine testGettingNameLenCharArray(this)
    class(test_dictionary), intent(inout)        :: this
    character(nameLen),dimension(:),allocatable  :: temp
    character(nameLen),dimension(:),pointer      :: temp_ptr => null()
    character(nameLen),dimension(1),parameter    :: default = ['Brasil, meu Brasil Brasileiro']
    logical(defBool)                             :: isSame

    call this % dict % get(temp,'charNameLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charNameLenArray == temp)
    @assertTrue(isSame, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'charNameLenArray', default)
    isSame = all(charNameLenArray == temp)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', default)
    isSame = all(default == temp)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

    !* With  pointer attribute
    call this % dict % get(temp_ptr,'charNameLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charNameLenArray == temp_ptr)
    @assertTrue(isSame, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp_ptr,'charNameLenArray', default)
    isSame = all(charNameLenArray == temp_ptr)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp_ptr,'invalid', default)
    isSame = all(default == temp_ptr)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

  end subroutine testGettingNameLenCharArray

!!
!! Test extracting pathLen long array of Chars
!!
@test
  subroutine testGettingPathLenCharArray(this)
    class(test_dictionary), intent(inout)        :: this
    character(pathLen),dimension(:),allocatable  :: temp
    character(pathLen),dimension(:),pointer      :: temp_ptr => null()
    character(pathLen),dimension(1),parameter    :: default = ['Brasil, meu Brasil Brasileiro']
    logical(defBool)                             :: isSame

    call this % dict % get(temp,'charPathLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charPathLenArray == temp)
    @assertTrue(isSame, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp,'charPathLenArray', default)
    isSame = all(charPathLenArray == temp)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp,'invalid', default)
    isSame = all(default == temp)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

    !* With pointer attribute
    call this % dict % get(temp_ptr,'charPathLenArray')
    ! Fun Fact. pFUnit does not support character arrays comparisons.
    ! Let Fortran handle comparisons
    isSame = all(charPathLenArray == temp_ptr)
    @assertTrue(isSame, 'Ordinary Retrival Failed')

    call this % dict % getOrDefault(temp_ptr,'charPathLenArray', default)
    isSame = all(charPathLenArray == temp_ptr)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

    call this % dict % getOrDefault(temp_ptr,'invalid', default)
    isSame = all(default == temp_ptr)
    @assertTrue(isSame, 'Get or Default Retrival Failed for Present Keyword')

  end subroutine testGettingPathLenCharArray

@test
  subroutine testGettingNestedDictionary(this)
    class(test_dictionary), intent(inout) :: this
    type(dictionary)                      :: temp
    real(defReal)                         :: tempReal
    integer(shortInt)                     :: tempInt
    character(nameLen)                    :: tempCharNameLen
    character(pathLen)                    :: tempCharPathLen
    real(defReal),dimension(:),allocatable        :: tempRealArray
    integer(shortInt), dimension(:), allocatable  :: tempIntArray
    character(nameLen), dimension(:), allocatable :: tempCharArrayNameLen
    character(pathLen), dimension(:), allocatable :: tempCharArrayPathLen
    logical(defBool)                      :: isSame

    call this % dict % get(temp,'nestedDict')

    ! Get all contents of the dictionary
    call temp % get(tempReal, 'myReal')
    call temp % get(tempint, 'myInt')
    call temp % get(tempCharNameLen, 'myCharNameLen')
    call temp % get(tempCharPathLen, 'myCharPathLen')
    call temp % get(tempRealArray, 'realArray')
    call temp % get(tempIntArray, 'intArray')
    call temp % get(tempCharArrayNameLen, 'charNameLenArray')
    call temp % get(tempCharArrayPathLen, 'charPathLenArray')

    ! Verify that content was not deformed
    isSame = tempReal == realVal
    isSame = isSame .and. tempInt == intVal
    isSame = isSame .and. tempCharNameLen == charNameLen
    isSame = isSame .and. tempCharPathLen == charPathLen
    isSame = isSame .and. all(tempIntArray == intArray)
    isSame = isSame .and. all(tempRealArray == realArray)
    isSame = isSame .and. all(tempCharArrayNameLen == charNameLenArray)
    isSame = isSame .and. all(tempCharArrayPathLen == charPathLenArray)

    @assertTrue(isSame, 'Contents of nested dictionary were changed')

  end subroutine testGettingNestedDictionary

  !!
  !!
  !!
@Test
  subroutine testPointerPassing(this)
    class(test_dictionary), intent(inout)        :: this
    type(dictionary)                      :: temp
    real(defReal)                         :: tempReal
    integer(shortInt)                     :: tempInt
    character(nameLen)                    :: tempCharNameLen
    character(pathLen)                    :: tempCharPathLen
    real(defReal),dimension(:),allocatable        :: tempRealArray
    integer(shortInt), dimension(:), allocatable  :: tempIntArray
    character(nameLen), dimension(:), allocatable :: tempCharArrayNameLen
    character(pathLen), dimension(:), allocatable :: tempCharArrayPathLen
    logical(defBool)                      :: isSame

    ! Retrieve pointer in a subroutine
    ! It is to test whether pointer to a nested dictionary
    ! going out of scope upon subroutine termination
    ! causes deallocation of memory.
    !
    call getAndFinalPointer(this % dict)

    call this % dict % get(temp,'nestedDict')

    ! Get all contents of the dictionary
    call temp % get(tempReal, 'myReal')
    call temp % get(tempint, 'myInt')
    call temp % get(tempCharNameLen, 'myCharNameLen')
    call temp % get(tempCharPathLen, 'myCharPathLen')
    call temp % get(tempRealArray, 'realArray')
    call temp % get(tempIntArray, 'intArray')
    call temp % get(tempCharArrayNameLen, 'charNameLenArray')
    call temp % get(tempCharArrayPathLen, 'charPathLenArray')

    ! Verify that content was not deformed
    isSame = tempReal == realVal
    isSame = isSame .and. tempInt == intVal
    isSame = isSame .and. tempCharNameLen == charNameLen
    isSame = isSame .and. tempCharPathLen == charPathLen
    isSame = isSame .and. all(tempIntArray == intArray)
    isSame = isSame .and. all(tempRealArray == realArray)
    isSame = isSame .and. all(tempCharArrayNameLen == charNameLenArray)
    isSame = isSame .and. all(tempCharArrayPathLen == charPathLenArray)

    @assertTrue(isSame, 'Contents of nested dictionary were changed')

  end subroutine testPointerPassing

  !!
  !!
  !!
  subroutine getAndFinalPointer(dict)
    class(dictionary), intent(in) :: dict
    class(dictionary),pointer     :: ptr
    logical(defBool)              :: myBool

    ptr => dict % getDictPtr('nestedDict')
    myBool = ptr % isPresent('bla')

  end subroutine getAndFinalPointer

!!
!! Test keys retrival
!!
@test
  subroutine testKeys(this)
    class(test_dictionary), intent(inout)        :: this
    character(nameLen),dimension(:),allocatable  :: tempAll
    character(nameLen),dimension(:),allocatable  :: tempReal
    character(nameLen),dimension(:),allocatable  :: tempInt
    character(nameLen),dimension(:),allocatable  :: tempChar
    character(nameLen),dimension(:),allocatable  :: tempRealArray
    character(nameLen),dimension(:),allocatable  :: tempIntArray
    character(nameLen),dimension(:),allocatable  :: tempCharArray
    character(nameLen),dimension(:),allocatable  :: tempDict
    character(nameLen)                           :: keyword
    logical(defBool) :: isValid

    ! Obtain all keys arrays
    call this % dict % keys(tempAll)
    call this % dict % keys(tempReal,'real')
    call this % dict % keys(tempInt,'int')
    call this % dict % keys(tempChar,'char')
    call this % dict % keys(tempDict,'dict')

    call this % dict % keys(tempRealArray,'realArray')
    call this % dict % keys(tempIntArray,'intArray')
    call this % dict % keys(tempCharArray,'charArray')

    ! Verify keys for real
    keyword = 'myReal'
    isValid = contains(tempAll, keyword) .and. contains(tempReal, keyword)
    @assertTrue(isValid, 'Keywords failed for real')

    ! Verify keys for int
    keyword = 'myInt'
    isValid = contains(tempAll, keyword) .and. contains(tempInt, keyword)
    @assertTrue(isValid,' Keywords failed for int')

    ! Verify keys for char
    keyword = 'myCharNameLen'
    isValid = contains(tempAll, keyword) .and. contains(tempChar, keyword)
    keyword = 'myCharPathLen'
    isValid = isValid .and. contains(tempAll, keyword) .and. contains(tempChar, keyword)
    @assertTrue(isValid,' Keywords failed for char')

    ! Verify keys for dict
    keyword = 'nestedDict'
    isValid = contains(tempAll, keyword) .and. contains(tempDict, keyword)
    @assertTrue(isValid,' Keywords failed for nested dictionary')

    ! Verify keys for realArray
    keyword = 'realArray'
    isValid = contains(tempAll, keyword) .and. contains(tempRealArray, keyword)
    @assertTrue(isValid, 'Keywords failed for real array ')

    ! Verify keys for int
    keyword = 'intArray'
    isValid = contains(tempAll, keyword) .and. contains(tempIntArray, keyword)
    @assertTrue(isValid,' Keywords failed for int array')

    ! Verify keys for char
    keyword = 'charNameLenArray'
    isValid = contains(tempAll, keyword) .and. contains(tempCharArray, keyword)
    keyword = 'charPathLenArray'
    isValid = isValid .and. contains(tempAll, keyword) .and. contains(tempCharArray, keyword)
    @assertTrue(isValid,' Keywords failed for char array')


  contains
    function contains(array,keyword) result(doesIt)
      character(nameLen), dimension(:) :: array
      character(*)                     :: keyword
      logical(defBool)                 :: doesIt

      doesIt = count(array == keyword) == 1

    end function contains
  end subroutine testKeys

  !!
  !! Test isPresent function of a dictionary
  !!
@test
  subroutine testIsPresent(this)
    class(test_dictionary), intent(inout) :: this
    logical(defBool)                      :: isPresent

    isPresent = this % dict % isPresent('nestedDict')
    @assertTrue(isPresent)

    isPresent = this % dict % isPresent('invalid')
    @assertFalse(isPresent)

  end subroutine testIsPresent

  !!
  !! Test getSize function of a dictionary
  !!
@test
  subroutine testGetSize(this)
    class(test_dictionary), intent(inout) :: this

    ! Get size of scalar
    @assertEqual(1, this % dict % getSize('myInt'))

    ! Get size of int Array
    @assertEqual(3, this % dict % getSize('intArray'))

    ! Get size of realArray
    @assertEqual(2, this % dict % getSize('realArray'))

    ! Get size of word Array
    @assertEqual(1, this % dict % getSize('charNameLenArray'))

    ! Get length of the dictionary (number of entries)
    @assertEqual(10, this % dict % length())

  end subroutine testGetSize


end module dictionary_test
