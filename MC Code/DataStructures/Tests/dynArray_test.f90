module dynArray_test
  use numPrecision
  use dynArray_class, only : dynIntArray
  use pFUnit_mod

  implicit none

contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour on unallocated dynArray
  !!
@Test
  subroutine testUnallocatedInt()
    type(dynIntArray) :: array

    @assertEqual(0, array % getSize(),'Size of unallocated array')
    @assertEqual(0, array % capacity(), 'Capacity of unallocated array')
    @assertTrue(array % isEmpty(),'isEmpty on unallocated array')

    ! Make shure no SEG ERR happens
    call array % shrink()

  end subroutine testUnallocatedInt

  !!
  !! Test correct behaviour on empty dynArray
  !!
@Test
  subroutine testEmptyInt()
    type(dynIntArray) :: array

    call array % add([1,2,3,4,5])
    call array % empty()

    @assertEqual(0, array % getSize(),'Size of empty array')
    @assertTrue(array % isEmpty(),'isEmpty on empty array')

    ! Make shure no SEG ERR happens
    call array % shrink()
    @assertEqual(0, array % capacity(),'Capacity of shrunk empty array')

  end subroutine testEmptyInt

  !!
  !! Test resizing
  !!
@Test
  subroutine testResizeInt()
    type(dynIntArray) :: array

    call array % resize(2)
    @assertLessThanOrEqual(2, array % capacity(),'Resize to 2 from 0')

    call array % resize(5)
    @assertLessThanOrEqual(5, array % capacity(),'Resize to 5 from 2')

    call array % resize(20)
    @assertLessThanOrEqual(20, array % capacity(),'Resize to 20 from 5 ')

    call array % kill()
    @assertEqual(0, array % capacity(),'Capacity of a killed array')

  end subroutine testResizeInt

  !!
  !! Test usage
  !!
@Test
  subroutine testUsage()
    type(dynIntArray) :: array
    integer(shortInt) :: i

    ! Build from vector by assignment
    array = [1,2,0,4,2]

    ! Test getting elements by index
    @assertEqual(1, array % get(1))
    @assertEqual(2, array % get(2))
    @assertEqual(0, array % get(3))
    @assertEqual(4, array % get(4))
    @assertEqual(2, array % get(5))

    call array % empty()

    ! Build by elements
    do i=1,7
      call array % add(2 * i)
    end do

    ! Test isEmpty and getSize
    @assertFalse(array % isEmpty())
    @assertEqual(7, array % getSize())

    ! Test popping
    do i=7,1,-1
      @assertEqual(2*i, array % pop())
    end do

   ! Test building by emelent and vector
   call array % add(1)
   call array % add([8,3])

   ! Test expose
   associate (ar => array % expose())
     @assertEqual(1, ar(1))
     @assertEqual(8, ar(2))
     @assertEqual(3, ar(3))
   end associate

   ! Shrink non-empty array
   call array % shrink()
   @assertEqual(3, array % capacity())

  end subroutine testUsage

end module dynArray_test
