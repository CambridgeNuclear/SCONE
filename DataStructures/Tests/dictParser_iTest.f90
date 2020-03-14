module dictParser_iTest
  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : fileToDict
  use pFUnit_mod
  implicit none

contains

  !!
  !! Test Reading a Dictionary from a File
  !!
@Test
  subroutine testFromFile()
    type(dictionary) :: dict
    integer(shortInt)  :: tempInt
    real(defReal)      :: tempReal
    character(nameLen) :: tempChar
    class(dictionary), pointer :: dictPtr
    integer(shortInt), dimension(:), allocatable  :: tempIntArray
    real(defReal), dimension(:), allocatable      :: tempRealArray
    character(nameLen), dimension(:), allocatable :: tempCharArray

    call fileToDict(dict,'./IntegrationTestFiles/testDictionary')

    ! Verify integer values
    call dict % get(tempInt, 'myInt')
    call dict % get(tempIntArray, 'intArray')

    @assertEqual(7, tempInt)
    @assertEqual([1, 2, 4, 5], tempIntArray)

    ! Verify real values
    call dict % get(tempReal, 'myReal')
    call dict % get(tempRealArray, 'realArray')

    @assertEqual(1.3_defReal, tempReal)
    @assertEqual([1.0_defReal, 2.2_defReal, 3.5_defReal], tempRealArray)

    ! Verify nested dictionary
    dictPtr => dict % getDictPtr('subDict')
    call dictPtr % get(tempInt, 'myInt')
    call dictPtr % get(tempReal, 'myReal')

    @assertEqual(3, tempInt)
    @assertEqual(3.2_defReal, tempReal)

  end subroutine testFromFile


end module dictParser_iTest
