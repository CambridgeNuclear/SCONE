module outputFile_inter

  use numPrecision
  use genericProcedures, only : fatalError
  use asciiOutput_inter, only : asciiOutput

  implicit none
  private

  character(*),parameter :: DEF_REAL_FORMAT = '(ES12.5)'
  character(*),parameter :: DEF_INT_FORMAT  = '(I12)'


  !!
  !! Interface for the file output.
  !!   It is used to allow easy swiching of output formats and
  !!   to allow writing to multiple files simultaneously
  !!
  !! Output is stream-like and not thread safe.
  !! Two kinds of outputs are allowed:
  !!   1) value  -> character, shortInt, longInt or defReal
  !!   2) result -> pair of defReals. 1st for result 2nd for STD
  !! Output can be scalar or multidimensional array (for now no limit on rank)
  !! Arrays are in column-major order (leftmost index varies fastest) (like in Fortran & Matlab)
  !!
  !! Output is divided into blocks, that can be nested.
  !! Each block has an associated name execpt for root block that has no name.
  !! Root block is opened at initialisation and closed at finalisation.
  !! Each value or result (scalar or array) is associated with a "name" (keyword)
  !! KEYWORDS AND SUB-BLOCK NAMES NEED TO BE UNIQUE WITHIN BLOCK!
  !! (***This rule may not be inforced in initial implementations-> requires "char set" )
  !!
  !! Output format is default(*) unless changed by the user (***this may change)
  !!
  type, public :: outputFile
    private
    class(asciiOutput),allocatable :: output

    ! State variables
    integer(shortInt)  :: blockLevel   = 0  ! Current depth in nested blocks
    integer(shortInt)  :: arrayElement = 0  ! Number of elements written to array
    integer(shortInt)  :: arrayLimit   = 0  ! Maximum number of elements that can be written to array
    character(:),allocatable :: current_block_name
    character(:),allocatable :: current_array_name


    ! Formats
    character(nameLen) :: real_format = DEF_REAL_FORMAT
    character(nameLen) :: int_format  = DEF_INT_FORMAT
  contains
    procedure :: init
    procedure :: startBlock


!   ! procedure :: init
!    ! Public interface
!    procedure(startBlock),deferred          :: startBlock
!    procedure(endBlock),deferred            :: endBlock
!    generic                                 :: printValue => printValue_defReal,  &
!                                                             printValue_shortInt, &
!                                                             printValue_longInt,  &
!                                                             printValue_char
!    procedure(printResult),deferred         :: printResult
!    procedure(startArray),deferred          :: startArray
!    procedure(endArray),deferred            :: endArray
!    procedure(addResult),deferred           :: addResult
!
!    ! Private procedures (not part of public interface)
!    procedure(printValue_defReal),deferred ,private :: printValue_defReal
!    procedure(printValue_shortInt),deferred,private :: printValue_shortInt
!    procedure(printValue_longInt),deferred,private  :: printValue_longInt
!    procedure(printValue_char),deferred,private     :: printValue_char
!
!    procedure(addValue_defReal),deferred ,private :: addValue_defReal
!    procedure(addValue_shortInt),deferred,private :: addValue_shortInt
!    procedure(addValue_longInt),deferred,private  :: addValue_longInt
!    procedure(addValue_char),deferred,private     :: addValue_char

  end type outputFile


!
!  abstract interface
!    !!
!    !!
!    !!
!    subroutine startBlock(self,name)
!      import :: outputFile
!      class(outputFile), intent(inout) :: self
!      character(*), intent(in)         :: name
!    end subroutine startBlock
!
!    !!
!    !!
!    !!
!    subroutine endBlock(self)
!      import :: outputFile
!      class(outputFile), intent(inout) :: self
!    end subroutine endBlock
!
!    !!
!    !!
!    !!
!    subroutine printValue_defReal(self,val,name)
!      import :: outputFile,&
!                defReal
!      class(outputFile), intent(inout) :: self
!      real(defReal),intent(in)         :: val
!      character(*), intent(in)         :: name
!    end subroutine printValue_defReal
!
!    !!
!    !!
!    !!
!    subroutine printValue_shortInt(self,val,name)
!      import :: outputFile,&
!                shortInt
!      class(outputFile), intent(inout) :: self
!      integer(shortInt),intent(in)     :: val
!      character(*), intent(in)         :: name
!    end subroutine printValue_shortInt
!
!    !!
!    !!
!    !!
!    subroutine printValue_longInt(self,val,name)
!      import :: outputFile,&
!                longInt
!      class(outputFile), intent(inout) :: self
!      integer(longInt),intent(in)      :: val
!      character(*), intent(in)         :: name
!    end subroutine printValue_longInt
!
!    !!
!    !!
!    !!
!    subroutine printValue_char(self,val,name)
!      import :: outputFile
!      class(outputFile), intent(inout) :: self
!      character(*), intent(in)         :: val
!      character(*), intent(in)         :: name
!    end subroutine printValue_char
!
!    !!
!    !!
!    !!
!    subroutine printResult(self,val,std,name)
!      import :: outputFile, &
!                defReal
!      class(outputFile), intent(inout) :: self
!      real(defReal), intent(in)        :: val
!      real(defReal), intent(in)        :: std
!      character(*), intent(in)         :: name
!    end subroutine printResult
!
!    !!
!    !!
!    !!
!    subroutine startArray(self,name,shape)
!      import :: outputFile, &
!                shortInt
!      class(outputFile), intent(inout) :: self
!      character(*), intent(in)         :: name
!      integer(shortInt),dimension(:)   :: shape
!    end subroutine startArray
!
!    !!
!    !!
!    !!
!    subroutine endArray(self)
!      import :: outputFile
!      class(outputFile), intent(inout) :: self
!    end subroutine endArray
!
!    !!
!    !!
!    !!
!    subroutine addValue_defReal(self,val)
!      import :: outputFile, &
!                defReal
!      class(outputFile), intent(inout) :: self
!      real(defReal), intent(in)        :: val
!    end subroutine addValue_defReal
!
!    !!
!    !!
!    !!
!    subroutine addValue_shortInt(self,val)
!      import :: outputFile, &
!                shortInt
!      class(outputFile), intent(inout) :: self
!      integer(shortInt), intent(in)    :: val
!    end subroutine addValue_shortInt
!
!    !!
!    !!
!    !!
!    subroutine addValue_longInt(self,val)
!      import :: outputFile, &
!                longInt
!      class(outputFile), intent(inout) :: self
!      integer(longInt), intent(in)    :: val
!    end subroutine addValue_longInt
!
!    !!
!    !!
!    !!
!    subroutine addValue_char(self,val)
!      import :: outputFile
!      class(outputFile), intent(inout) :: self
!      character(*), intent(in)         :: val
!    end subroutine addValue_char
!
!    !!
!    !!
!    !!
!    subroutine addResult(self,val,std)
!      import :: outputFile, &
!                defReal
!      class(outputFile), intent(inout) :: self
!      real(defReal), intent(in)        :: val
!      real(defReal),intent(in)         :: std
!    end subroutine addResult
!
!
!  end interface
!



contains
  !!
  !!
  !!
  subroutine init(self,type)
    class(outputFile), intent(inout) :: self
    character(*), intent(in)         :: type

  end subroutine init


  !!
  !! Start new nested block with name: "name"
  !!
  subroutine startBlock(self,name)
    class(outputFile), intent(inout) :: self
    character(*), intent(in)         :: name
    character(100), parameter :: Here ='startBlock (outputFile_inter.f90)'

    !*** Check that name is unique in current block
    !*** Add name to keywords occupied in current block

    ! Check that currently is not writing array
    if (self % arrayElement > 0 .and. self % arrayLimit > 0) then
      call fatalError(Here,'In block: '// self % current_block_name // &
                          ' Cannot begin new block: '//name//          &
                          ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Update state
    self % blockLevel = self % blockLevel + 1


    ! Start new block in output
    call self % output % startBlock(name)

  end subroutine startBlock

    
end module outputFile_inter
