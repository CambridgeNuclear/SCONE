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
  !! Assume that checks for mixed arrays are performed by the user class.
  !! Printer should return errors if calls for ends are out of correct sequence
  !!
  type, public,abstract :: asciiOutput
    private
  contains
    procedure(writeToFile),deferred :: writeToFile
    procedure(startBlock),deferred  :: startBlock
    procedure(endBlock),deferred    :: endBlock
    procedure(startEntry),deferred  :: startEntry
    procedure(endEntry),deferred    :: endEntry
    procedure(startArray),deferred  :: startArray
    procedure(endArray),deferred    :: endArray
    procedure(printNum),deferred    :: printNum
    procedure(printChar),deferred   :: printChar
  end type asciiOutput


  abstract interface
    !!
    !! For now it prints to screen for debug
    !!
    subroutine writeToFile(self,unit)
      import :: asciiOutput, &
                shortInt
      class(asciiOutput), intent(inout) :: self
      integer(shortInt), intent(in)     :: unit
    end subroutine writeToFile

    !!
    !! Change state to writing new block with "name"
    !!
    subroutine startBlock(self,name)
      import :: asciiOutput, &
                nameLen
      class(asciiOutput), intent(inout) :: self
      character(nameLen), intent(in)    :: name
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
      import :: asciiOutput, &
                nameLen
      class(asciiOutput), intent(inout) :: self
      character(nameLen), intent(in)    :: name
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
    !! Name should alrady be provided by "startEntry"
    !!
    subroutine startArray(self,shape)
      import :: asciiOutput, &
                shortInt
      class(asciiOutput), intent(inout) :: self
      integer(shortInt),dimension(:)    :: shape
    end subroutine startArray

    !!
    !! End writing array
    !!
    subroutine endArray(self)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
    end subroutine endArray

    !!
    !! Print val assuming it contains valid printed number with no leading or trailing blanks
    !!
    subroutine printNum(self,val)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
      character(*),intent(in)           :: val
    end subroutine printNum

    !!
    !! Print val assuming it contains valid printed string with no leading or trailing blanks
    !!
    subroutine printChar(self,val)
      import :: asciiOutput
      class(asciiOutput), intent(inout) :: self
      character(*),intent(in)           :: val
    end subroutine printChar
  end interface

end module asciiOutput_inter
