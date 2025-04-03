module linkedList_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  !! Presently only contains a linked list of shortInts. Can be easily duplicated
  !! for other variable types.

  !!
  !! Linked list node containing an int and pointer
  !! to the next node
  !!
  type :: intNode
    integer(shortInt) :: data
    type(intNode), pointer :: next => null()
  end type intNode


  !!
  !! Linked list for shortIntegers
  !!
  !! Contains the start of the list and end (for faster additions).
  !! Could be optimised by returning the node as an iterator. 
  !! The current implementation of traversal will start to become slow
  !! for searching long lists.
  !!
  !! Public interface:
  !! add     -> add a new entry to the list (omp critical)
  !! get     -> get a value in the list at a particular index
  !! getSize -> return the length of the list
  !! display -> output the contents of the list to display
  !! kill    -> clean up the list
  !!
  !! Private procedures:
  !! traverse -> move to a specified node in the list
  !!
  type, public :: linkedIntList
    private
    type(intNode), pointer :: head => null()
    type(intNode), pointer :: tail => null()
    integer(shortInt)      :: length = 0
  contains
    ! Public interface
    procedure  :: add => add_shortInt
    procedure  :: get => get_shortInt
    procedure  :: getSize => getSize_shortInt
    procedure  :: display => display_shortInt
    procedure  :: kill => kill_shortInt

    ! Private procedures
    procedure, private  :: traverse => traverse_shortInt

  end type linkedIntList

contains

  !!
  !! Add entry at the end of linked list
  !!
  !! Enclosed in a critical statement to prevent corruption.
  !!
  subroutine add_shortInt(self, entry)
    class(linkedIntList), intent(inout) :: self
    integer(shortInt), intent(in)       :: entry
    class(intNode), pointer             :: tail

    !$omp critical
    ! Go to the end of the list
    if (self % length == 0) then
      allocate(self % head)
      self % head % data = entry
      self % tail => self % head

    else
      tail => self % tail
      allocate(tail % next)
      self % tail => tail % next
      self % tail % data = entry
    end if

    self % length = self % length + 1
    !$omp end critical

  end subroutine add_shortInt

  !!
  !! Move down the linked list until a given index
  !!
  !! Errors:
  !!   Requested index exceeds the length of the list
  !!
  function traverse_shortInt(self, idx) result(resNode)
    class(linkedIntList), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    class(intNode), pointer          :: resNode
    integer(shortInt)                :: i
    character(100), parameter :: Here = 'traverse_shortInt (linkedList_Class.f90)'

    if (idx > self % length) call fatalError(Here,&
            'Requested index exceeds list length')

    resNode => self % head
    i = 1
    do while (i < idx)
      resNode => resNode % next
      i = i + 1
    end do

  end function traverse_shortInt

  !!
  !! Get size of the linked list
  !!
  pure function getSize_shortInt(self) result(S)
    class(linkedIntList), intent(in) :: self
    integer(shortInt)                :: S

    S = self % length

  end function getSize_shortInt

  !!
  !! Return a value from the list at a given index
  !!
  !! Errors:
  !!   The list has no elements
  !!
  function get_shortInt(self, idx) result(res)
    class(linkedIntList), intent(in) :: self
    integer(shortInt), intent(in)    :: idx
    integer(shortInt)                :: res
    class(intNode), pointer          :: resNode
    character(100), parameter :: Here = 'get_shortInt (linkedList_Class.f90)'

    if (self % getSize() == 0) call fatalError(Here,'Linked list is not allocated')
    
    resNode => self % traverse(idx)
    res = resNode % data

  end function get_shortInt

  !!
  !! Deallocate linked list
  !!
  subroutine kill_shortInt(self)
    class(linkedIntList), intent(inout) :: self
    integer(shortInt)                   :: i
    class(intNode), pointer             :: resNode

    ! Traverse the list and nullify pointers
    do while (self % length > 1)
      self % length = self % length - 1
      resNode => self % traverse(self % length)
      resNode % next => null()
    end do
    self % head => null()
    self % tail => null()
    self % length = 0

  end subroutine kill_shortInt

  !!
  !! Print contents of linked list
  !!
  subroutine display_shortInt(self)
    class(linkedIntList), intent(inout) :: self
    integer(shortInt)                   :: i
    class(intNode), pointer             :: resNode

    if (self % length == 0) then
      print *,'Empty list'
      return
    end if

    resNode => self % head
    i = 1
    do while (i < self % length)
      print *, resNode % data
      resNode => resNode % next
      i = i + 1
    end do

  end subroutine display_shortInt

end module linkedList_class
