module asciiOutput_inter

  use numPrecision

  implicit none
  private

  !!
  !! Abstract interface for the printing output in stream fashion to an ASCII file
  !!   It allows to swap diffrent printers to print to diffrent formats (i.e. MATLAB, CSV, JASON)
  !!
  !! Output is divided into blocks, that can be nested.
  !! Each block has an associated name execpt for root block that has no name.
  !! Root block is opened at initialisation and closed at finalisation.
  !! Each entry is associated with a "name" (keyword)
  !!  Following functionality is provided:
  !!    -> startBlock(name) - starts new block with the given name
  !!    -> endBlock         - ends block (return error if trying to  end root block)
  !!    -> startEntry(name) - start writing entry with name. Expect either array of single value
  !!    -> endEntry         - end writing entry
  !!    -> startArray(name,shape) - start writing array with given shape
  !!    -> endArray               - end writing array
  !!    -> print(char,mold)       - print char asuming it has the same type as mold(defReal, etc. )
  !!
  !! Assume that checks for name uniqueness are performed by the user class.
  !! Printer should return errors if calls for ends are out of correct sequence
  !!
  type, public,abstract :: asciiOutput
    private
  contains
    procedure(startBlock),deferred :: startBlock
    procedure(endBlock),deferred   :: endBlock
    procedure(startEntry),deferred :: startEntry
    procedure(endEntry),deferred   :: endEntry
    procedure(startArray),deferred :: startArray
    procedure(endArray),deferred   :: endArray
    generic                        :: print => print_defReal,  &
                                               print_shortInt, &
                                               print_longInt , &
                                               print_char

    procedure(print_defReal),deferred,private  :: print_defReal
    procedure(print_shortInt),deferred,private :: print_shortInt
    procedure(print_longInt),deferred,private  :: print_longInt
    procedure(print_char),deferred,private     :: print_char
  end type asciiOutput


  abstract interface
    !!
    !! Change state to writing new block with "name"
    !!
    subroutine startBlock(self,name)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
      character(*), intent(in)          :: name
    end subroutine startBlock

    !!
    !! End top level block and return to previous block
    !! Return error if this is called in root block
    !!
    subroutine endBlock(self)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
    end subroutine endBlock

    !!
    !! Change state to writing a new entry
    !! Can recive single value or array next
    !!
    subroutine startEntry(self,name)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
      character(*), intent(in)          :: name
    end subroutine startEntry

    !!
    !! End writing a new entry
    !!
    subroutine endEntry(self)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
    end subroutine endEntry

    !!
    !! Start writing array with shape & column-major order(leftmost index varies fastest)
    !!
    subroutine startArray(self,name,shape)
      import :: asciiOutput, &
                shortInt
      class(asciiOutput), intent(inout) :: self
      character(*), intent(in)           :: name
      integer(shortInt),dimension(:)     :: shape
    end subroutine startArray

    !!
    !! End writing array
    !!
    subroutine endArray(self)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
    end subroutine endArray

    !!
    !! Print val assuming it contains valid printed real with no leading or trailing blanks
    !!
    subroutine print_defReal(self,val,mold)
      import :: asciiOutput, &
                defReal
      class(asciiOutput), intent(inout) :: self
      character(*),intent(in)           :: val
      real(defReal), intent(in)         :: mold
    end subroutine print_defReal

    !!
    !! Print val assuming it contains valid printed shortInt with no leading or trailing blanks
    !!
    subroutine print_shortInt(self,val,mold)
      import :: asciiOutput, &
                shortInt
      class(asciiOutput), intent(inout) :: self
      character(*),intent(in)           :: val
      integer(shortInt), intent(in)     :: mold
    end subroutine print_shortInt

    !!
    !! Print val assuming it contains valid printed longInt with no leading or trailing blanks
    !!
    subroutine print_longInt(self,val,mold)
      import :: asciiOutput, &
                longInt
      class(asciiOutput), intent(inout) :: self
      character(*),intent(in)           :: val
      integer(longInt), intent(in)      :: mold
    end subroutine print_longInt

    !!
    !! Print val assuming it contains valid printed string with no leading or trailing blanks
    !!
    subroutine print_char(self,val,mold)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
      character(*),intent(in)           :: val
      character(*),intent(in)           :: mold
    end subroutine print_char
  end interface

end module asciiOutput_inter
