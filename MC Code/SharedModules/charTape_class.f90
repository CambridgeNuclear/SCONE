module charTape_class

  use numPrecision

  implicit none
  private

  ! Growth ratio of allocated memory (MUST BE > 1 !!!)
  real(defReal), parameter :: GROWTH_RATIO =2.0_defReal


  !!
  !! Tape of characters (most likley) similar to C++ string.
  !! Avoids unnecessary memory reallocations
  !!
  type, public :: charTape
    private
    integer(shortInt)        :: end    = 0
    integer(shortInt)        :: leng   = 0
    character(:),allocatable :: tape
  contains
    generic    :: append => append_tape, append_char
    procedure  :: expose
    procedure  :: cut
    procedure  :: length

    procedure, private :: append_tape
    procedure, private :: append_char
    procedure, private :: resize
  end type charTape

contains

  !!
  !! Return contents of tape as a character
  !!
  function expose(self) result(tape)
    class(charTape), intent(in) :: self
    character(self % end)       :: tape

    tape = self % tape(1:self % end)

  end function expose

  !!
  !! Returns length of the string
  !!
  function length(self) result(L)
    class(charTape), intent(in) :: self
    integer(shortInt)           :: L

    L = self % end

  end function length

  !!
  !! Cut N last characters from the tape
  !!
  subroutine cut(self,N)
    class(charTape), intent(inout) :: self
    integer(shortInt)              :: N

    self % end = max(0, self % end - N)

  end subroutine cut

  !!
  !! Appends self with tape on the RHS
  !!
  subroutine append_tape(self,tape)
    class(charTape), intent(inout) :: self
    type(charTape), intent(in)     :: tape
    integer(shortInt)              :: end

    end = self % end

    ! Protect against degenerate tape
    if ( tape % end > 0 .and. allocated(tape % tape)) then
      ! Resize if necessary
      call self % resize(end + tape % end)

      ! Extend tape
      self % tape(end +1 : end + tape % end) = tape % tape(1:tape % end)

      ! Update current end
      self % end = end + tape % end
    end if

  end subroutine append_tape

  !!
  !! Appends tape with chara on the RHS
  !!
  subroutine append_char(self,chara)
    class(charTape), intent(inout) :: self
    character(*), intent(in)       :: chara
    integer(shortInt)              :: end

    end = self % end

    ! Resize if necessary
    call self % resize(end + len(chara) )

    ! Extend tape
    self % tape(end + 1:end + len(chara)) = chara

    ! Update current end
    self % end = end + len(chara)

  end subroutine append_char


  !!
  !! Given desired length N
  !!   Do nothing if N < len(self % tape)
  !!   Allocate to new length ceiling(N * GROWTH_RATIO) otherwise
  !!
  subroutine resize(self,N)
    class(charTape), intent(inout) :: self
    integer(shortInt), intent(in)  :: N
    character(:),allocatable       :: temp

    ! Return if length is OK
    if( N < self % leng) return

    self % leng = ceiling(N * GROWTH_RATIO)

    ! Extend or allocate tape
    if(allocated(self % tape)) then
      allocate( character(self % leng) :: temp)
      temp(1:self % end) = self % tape
      call move_alloc(temp, self % tape)

    else
      allocate( character(self % leng) :: self % tape)

    end if
  end subroutine resize

    
end module charTape_class
