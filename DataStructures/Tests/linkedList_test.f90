module linkedList_test
  use numPrecision
  use linkedList_class, only : linkedIntList
  use funit

  implicit none

contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour on unallocated linkedList
  !!
@Test
  subroutine testUnallocatedInt()
    type(linkedIntList) :: list

    @assertEqual(0, list % getSize(),'Size of unallocated list')

  end subroutine testUnallocatedInt

  !!
  !! Test correct behaviour on empty linkedIntList
  !!
@Test
  subroutine testEmptyInt()
    type(linkedIntlist) :: list

    call list % add(5)
    call list % add(3)
    call list % add(4)

    call list % kill()

    @assertEqual(0, list % getSize(),'Size of empty list')

  end subroutine testEmptyInt

  !!
  !! Test resizing
  !!
@Test
  subroutine testAddInt()
    type(linkedIntList) :: list

    call list % add(7)
    @assertEqual(1, list % getSize())

    call list % add(4)
    call list % add(6)
    @assertEqual(3, list % getSize())

    call list % kill()
    @assertEqual(0, list % getSize())

  end subroutine testAddInt

  !!
  !! Test usage
  !!
@Test
  subroutine testUsage()
    type(linkedIntList) :: list
    integer(shortInt)   :: i

    call list % add(8)
    call list % add(42)
    call list % add(666)

    ! Test getting elements by index
    @assertEqual(8, list % get(1))
    @assertEqual(42, list % get(2))
    @assertEqual(666, list % get(3))

    call list % kill()

    ! Build by elements
    do i=1,7
      call list % add(2 * i)
    end do

    ! Test getSize
    @assertFalse(list % getSize() == 0)
    @assertEqual(7, list % getSize())


  end subroutine testUsage

  !!
  !! Test kill of unallocated array
  !!
@Test
  subroutine testKillUnalloc()
    type(linkedIntList) :: list

    call list % kill()

  end subroutine testKillUnalloc

end module linkedList_test
