module charTape_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar

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
    procedure  :: segment
    procedure  :: scanFrom
    procedure  :: get
    procedure  :: cut
    procedure  :: length
    procedure  :: clean

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
  !! Wrapper around scan intrinsic procedure for charTape
  !!
  !! Args:
  !!   start [in] -> Starting location for the search
  !!   set [in] -> Set of chacracters to find
  !!
  !! Result:
  !!   Location of the fisrs symbol in 'set' relative to start location
  !!   Returns 0 if no symbol is found
  !!
  !! Errors:
  !!   fatalError if start is not +ve or > then length
  !!
  function scanFrom(self, start, set) result(pos)
    class(charTape), intent(in)   :: self
    integer(shortInt), intent(in) :: start
    integer(shortInt)             :: pos
    character(*), intent(in)      :: set
    character(100), parameter :: Here = 'scan (charTape_class.f90)'

    ! Deal with errors
    if(start <= 0) then
      call fatalError(Here, 'Start cannot be -ve. Wa given: '//numToChar(start))
    elseif(start > self % length()) then
      call fatalError(Here,'Start: '//numToChar(start)//' is beyond length of tape: '//&
                            numToChar(self % length()))
    end if

    ! Return result
    pos = scan(self % tape(start:self % length()), set)

  end function scanFrom

  !!
  !! Copy return a segment of a charTape
  !!
  !! Args:
  !!   start [in] -> Starting location
  !!   end [in]   -> Ending location
  !!
  !! Result:
  !!   Character of length = end-start. Contains the section of the tape
  !!
  !! Errors:
  !!   fatalError if start is 0 or 0-ve
  !!   fatalErrror if end < start (Unless segments on creation of c with -ve len)
  !!
  function segment(self, start, end) result(c)
    class(charTape), intent(in)   :: self
    integer(shortInt), intent(in) :: start
    integer(shortInt), intent(in) :: end
    character(end-start + 1)      :: c
    character(100), parameter :: Here = 'segment (charTape_class.f90)'

    ! Process potential erros
    if (start <= 0) then
      call fatalError(Here, 'Start is not +ve: '//numToChar(start))
    elseif (end < start) then
      call fatalError(Here,'End: '//numToChar(end)//' is < to start: '//numToChar(start))
    elseif (end > self % length()) then
      call fatalError(Here,'End: '//numToChar(end)//' is beyond length of tape: '//&
                            numToChar(self % length()))
    end if

    ! Return segment
    c = self % tape(start:end)

  end function segment


  !!
  !! Return character at poistion pos
  !!
  !! Args:
  !!   pos [in] -> Integer. Position to get character from
  !!
  !! Result:
  !!   Character at position pos
  !!
  !! Errors:
  !!   fatalError if pos > length of the tape
  !!   fatalError if pox <= 0
  !!
  function get(self, pos) result(c)
    class(charTape), intent(in)   :: self
    integer(shortInt), intent(in) :: pos
    character(1)                  :: c
    character(100), parameter :: Here ='get (charTape_class.f90)'

     ! Error messages
     if (pos <= 0) then
       call fatalError(Here,'Position cannot be -ve. Was given '//numToChar(pos))
     elseif (pos > self % length()) then
       call fatalError(Here, 'Position: '//numToChar(pos)//&
                             ' is larger then tape length: '//numToChar(self % length()))
     end if

     c = self % tape(pos:pos)

  end function get


  !!
  !! Returns length of the string
  !!
  function length(self) result(L)
    class(charTape), intent(in) :: self
    integer(shortInt)           :: L

    L = self % end

  end function length

  !!
  !! Remove content without deallocation
  !!
  subroutine clean(self)
    class(charTape), intent(inout) :: self

    self % end = 0

  end subroutine clean

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
