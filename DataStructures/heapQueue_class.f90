module heapQueue_class

  use numPrecision
  use genericProcedures, only : swap
  use errors_mod,        only : fatalError

  implicit none
  private

  !!
  !! A fixed size heap queue designed to store N smallest values of a sample
  !!
  !! The queue is implemented as a binary heap, that is a binary tree with the
  !! heap property. The heap property is that the value of a parent node is
  !! larger than the value of both its children. This means that the largest
  !! value is always at the root of the tree.
  !!
  !! The queue is implemented as an array with the root at index 1. The children
  !! of a node at index i are at 2*i and 2*i + 1. Note that for the even number
  !! of elements the last node will only have one child.
  !!
  !! This data structure is intended to be used for sampling without replacement
  !! as it allows to find a value of a threshold that selects K values from a
  !! stream of the N random numbers.
  !!
  !! Interface:
  !!   init        -> Initialise the queue to a given size
  !!   pushReplace -> Add a value to the queue either growing the size or replacing the largest element
  !!                  if maximum size was reached
  !!   maxValue    -> Returns the largest value in the queue
  !!   getSize     -> Returns the current size of the queue
  !!
  type, public :: heapQueue
    private
    real(defReal), dimension(:) , allocatable :: heap
    integer(shortInt)                         :: size
  contains
    procedure :: init
    procedure :: pushReplace
    procedure :: maxValue
    procedure :: getSize

    procedure, private :: push
    procedure, private :: replace
  end type heapQueue

contains

  !!
  !! Initialise the queue to a given size
  !!
  !! Args:
  !!  maxSize [in] -> Maximum size of the queue
  !!
  subroutine init(self, maxSize)
    class(heapQueue), intent(out) :: self
    integer(shortInt), intent(in) :: maxSize

    self % size = 0
    allocate(self % heap(maxSize))

  end subroutine init

  !!
  !! Add a value to to queue either growing the size or replacing the largest
  !!
  !! Args:
  !!  value [in] -> Value to add to the queue
  !!
  subroutine pushReplace(self, val)
    class(heapQueue), intent(inout) :: self
    real(defReal), intent(in)       :: val

    if (self % size < size(self % heap)) then
      call self % push(val)
    else
      call self % replace(val)
    end if

  end subroutine pushReplace

  !!
  !! Add a value to the queue
  !!
  !! Assumes enough space is available
  !!
  !! Args:
  !!  value [in] -> Value to add to the queue
  !!
  subroutine push(self, val)
    class(heapQueue), intent(inout) :: self
    real(defReal), intent(in)       :: val
    integer(shortInt)               :: parent, child

    ! Increase the size of the queue and add the new value
    self % size = self % size + 1
    self % heap(self % size) = val

    ! If the heap is of size 1 there is no need to order it
    ! Also, avoid test fail in debug mode, since parent would be 0 and the
    ! code is looking for self % heap(parent) inside the while condition
    if (self % size == 1) return

    ! Shift the new value up the heap to restore the heap property
    child = self % size
    parent = child / 2

    do while (self % heap(parent) < self % heap(child))
      call swap(self % heap(parent), self % heap(child))
      child = parent
      parent = child / 2
      ! As above: avoid error in debug mode, caused by trying to access self % heap(0)
      if (parent == 0) return
    end do

  end subroutine push

  !!
  !! Replaces the largest value in the queue with a new value
  !!
  !! Args:
  !!  value [in] -> Value to add to the queue
  !!
  subroutine replace(self, val)
    class(heapQueue), intent(inout) :: self
    real(defReal), intent(in)       :: val
    integer(shortInt)               :: parent, child

    self % heap(1) = val

    parent = 1
    child = 2

    ! Shift down the new value until heap property is restored be comparing
    ! with the largest child and swapping if necessary
    do while (child <= self % size)

      ! We need to consider the case where there is only one child in a node
      ! (can happen when the size is even). If child is the last element it is the
      ! node with no sibling
      if (child /= self % size .and. self % heap(child) < self % heap(child + 1)) then
        child = child + 1
      end if

      ! If the parent is larger than the larger child we are done
      if (self % heap(parent) >= self % heap(child)) then
        return
      end if

      ! Otherwise swap the parent with the larger child and continue
      ! the recursion
      call swap(self % heap(parent), self % heap(child))
      parent = child

      ! Child points to the start of next level
      child = parent * 2
    end do

  end subroutine replace

  !!
  !! Returns the largest value in the queue
  !!
  !! Errors:
  !!   fatal error if the queue is empty
  !!
  function maxValue(self) result(val)
    class(heapQueue), intent(in) :: self
    real(defReal)                :: val
    character(100), parameter :: Here = "maxValue (heapQueue_class.f90)"

    if (self % size == 0) then
      call fatalError(Here, "The queue is empty!")
    end if

    val = self % heap(1)

  end function maxValue

  !!
  !! Get the current size of the queue
  !!
  pure function getSize(self) result(size)
    class(heapQueue), intent(in) :: self
    integer(shortInt)            :: size

    size = self % size

  end function getSize


end module heapQueue_class
