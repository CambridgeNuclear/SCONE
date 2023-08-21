module dictParser_test
  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use funit
  implicit none

contains

  !!
  !! Sets up test_dictionary object we can use in a number of tests
  !!
@Test
  subroutine testFromChar()
    type(dictionary)   :: dict
    integer(shortInt)  :: tempInt
    real(defReal)      :: tempReal
    character(nameLen) :: tempChar
    class(dictionary), pointer :: dictPtr
    integer(shortInt), dimension(:), allocatable  :: tempIntArray
    real(defReal), dimension(:), allocatable      :: tempRealArray
    character(nameLen), dimension(:), allocatable :: tempCharArray
    character(*),parameter :: tape = " myInt 7;                                 &
                                       myChar my;                               &
                                       myReal 1.3;                              &
                                       weirdFloat 1E-11;                        &
                                       intArray (1 2 4 5);                      &
                                       realArray (1.1 2.2 3.4 1E-11);           &
                                       charArray (One element );                &
                                       subDict { myInt 3; myReal 3.2; }"


    ! Create dictionary
    call charToDict(dict, tape)

    ! Verify integer values
    call dict % get(tempInt, 'myInt')
    call dict % get(tempIntArray, 'intArray')

    @assertEqual(7, tempInt)
    @assertEqual([1, 2, 4, 5], tempIntArray)

    ! Verify real values
    call dict % get(tempReal, 'myReal')
    call dict % get(tempRealArray, 'realArray')

    @assertEqual(1.3_defReal, tempReal)
    @assertEqual([1.1_defReal, 2.2_defReal, 3.4_defReal, 1.0E-11_defReal], tempRealArray)

    ! Verify a problematic value
    ! It may not be parsed correclty with wrong fromat settings
    call dict % get(tempReal, 'weirdFloat')
    @assertEqual(1.0E-11_defReal, tempReal)

    ! Verify nested dictionary
    dictPtr => dict % getDictPtr('subDict')
    call dictPtr % get(tempInt, 'myInt')
    call dictPtr % get(tempReal, 'myReal')

    @assertEqual(3, tempInt)
    @assertEqual(3.2_defReal, tempReal)

  end subroutine testFromChar


end module dictParser_test
