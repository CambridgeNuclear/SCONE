module outputFile_class

  use numPrecision
  use genericProcedures,       only : fatalError, numToChar
  use stack_class,             only : stackChar
  use asciiOutput_inter,       only : asciiOutput
  use asciiOutputFactory_func, only : new_asciiOutput
  use iso_fortran_env,         only : OUTPUT_UNIT

  implicit none
  private

  character(*),parameter :: DEF_REAL_FORMAT = '(ES12.5)'
  character(*),parameter :: DEF_INT_FORMAT  = '(I12)'

  integer(shortInt), parameter :: NOT_ARRAY      = 0, &
                                  UNDEF_ARRAY     = 1, &
                                  RES_ARRAY      = 2, &
                                  VAL_ARRAY_REAL = 3, &
                                  VAL_ARRAY_INT  = 4, &
                                  VAL_ARRAY_CHAR = 5

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
  !! Names of the arrays and block need to fit into nameLen
  !!
  type, public :: outputFile
    private
    class(asciiOutput),allocatable :: output

    ! State variables
    integer(shortInt)  :: blockLevel    = 0  ! Current depth in nested blocks
    integer(shortInt)  :: arrayTop      = 0  ! Number of elements written to array
    integer(shortInt)  :: arrayLimit    = 0  ! Maximum number of elements that can be written to array
    integer(shortInt)  :: arrayType     = NOT_ARRAY
    character(nameLen) :: current_block_name =''
    character(nameLen) :: current_array_name =''
    type(stackChar)    :: block_name_stack

    ! Buffors
    integer(shortInt),dimension(:), allocatable :: shapeBuffer

    ! Formats
    character(nameLen) :: real_format = DEF_REAL_FORMAT
    character(nameLen) :: int_format  = DEF_INT_FORMAT

  contains
    procedure :: init
    procedure :: writeToFile
    procedure :: writeToConsole

    ! Block writing interface
    procedure :: startBlock
    procedure :: endBlock

    ! Single entry writing interface
    procedure :: printResult
    generic   :: printValue => printValue_defReal, printValue_shortInt, &
                               printValue_longInt, printValue_char
    ! Array writing interface
    procedure :: startArray
    procedure :: endArray
    generic   :: addResult => addResult_scalar, addResult_rank1
    generic   :: addValue  => addValue_defReal_scalar, addValue_defReal_rank1,  &
                              addValue_shortInt_scalar, addValue_shortInt_rank1,&
                              addValue_longInt_scalar, addValue_longInt_rank1,  &
                              addValue_char_scalar, addValue_char_rank1

    ! Private procedures to write single entry
    procedure,private :: printValue_defReal
    procedure,private :: printValue_shortInt
    procedure,private :: printValue_longInt
    procedure,private :: printValue_char

    ! Private array writing procedures
    procedure,private :: addResult_scalar
    procedure,private :: addResult_rank1
    procedure,private :: addValue_defReal_scalar
    procedure,private :: addValue_defReal_rank1
    procedure,private :: addValue_shortInt_scalar
    procedure,private :: addValue_shortInt_rank1
    procedure,private :: addValue_longInt_scalar
    procedure,private :: addValue_longInt_rank1
    procedure,private :: addValue_char_scalar
    procedure,private :: addValue_char_rank1

    ! Function to print numbers using format set in the given  output_file
    generic :: num2Char => num2Char_defReal, num2Char_shortInt, num2Char_longInt
    procedure, private :: num2Char_defReal
    procedure, private :: num2Char_shortInt
    procedure, private :: num2Char_longInt
  end type outputFile

contains
  !!
  !! Initialise output file by providing format type
  !!
  subroutine init(self,type)
    class(outputFile), intent(inout) :: self
    character(nameLen),intent(in)    :: type

    allocate( self % output, source = new_asciiOutput(type))

  end subroutine init

  !!
  !! Dump collected results to file
  !!
  subroutine writeToFile(self,file)
    class(outPutFile), intent(inout) :: self
    character(pathLen), intent(in)   :: file
    integer(shortInt)                :: unitNum
    integer(shortInt)                :: error
    character(99)                    :: errorMsg
    character(100), parameter :: Here ='writeToFile ( outputFile_class.f90)'

    ! Open file to write
    open ( newunit   = unitNum, &
           file      = file   , &
           action    = 'write', &
           iostat    = error  , &
           iomsg     = errorMsg )

    ! Catch error
    if (error /= 0 ) call fatalError(Here, errorMsg)

    ! Write to file
    call self % output % writeToFile(unitNum)

    ! Close file
    close(unitNum)

  end subroutine writeToFile

  !!
  !! Dump Collected Results to Console (default output)
  !!
  subroutine writeToConsole(self)
    class(outputFile), intent(inout) :: self

    call self % output % writeToFile(OUTPUT_UNIT)

  end subroutine writeToConsole


  !!
  !! Start new nested block with name: "name"
  !!
  subroutine startBlock(self,name)
    class(outputFile), intent(inout) :: self
    character(nameLen), intent(in)   :: name
    character(100), parameter :: Here ='startBlock (outputFile_class.f90)'

    !*** Check that name is unique in current block
    !*** Add name to keywords occupied in current block

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call fatalError(Here,'In block: '// self % current_block_name // &
                          ' Cannot begin new block: '//name//          &
                          ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Update state
    self % blockLevel = self % blockLevel + 1
    call self % block_name_stack % push( self % current_block_name)
    self % current_block_name = name

    ! Start new block in output
    call self % output % startBlock(name)

  end subroutine startBlock

  !!
  !! Exit from the current block.
  !!
  subroutine endBlock(self)
    class(outputFile), intent(inout) :: self
    character(100), parameter :: Here ='endBlock (outputFile_class.f90)'

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call fatalError(Here,'In block: '// self % current_block_name // &
                          ' Cannot exit from block: '//                &
                          ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Check that is not in root block
    if ( self % blockLevel == 0 ) call fatalError(Here,'Trying to exit from root block')

    ! Update state
    self % blockLevel = self % blockLevel - 1
    call self % block_name_stack % pop( self % current_block_name)

    ! End block in output
    call self % output % endBlock()

  end subroutine endBlock

  !!
  !! Prepare for writing array. Store shape information
  !! We don't know the type of content yet so we do not inform printer about shape
  !!
  subroutine startArray(self,name,shape)
    class(outputFile), intent(inout) :: self
    character(nameLen), intent(in)   :: name
    integer(shortInt),dimension(:)   :: shape
    character(100), parameter :: Here ='startArray (outputFile_class.f90)'

    !*** Check that name is unique in current block
    !*** Add name to keywords occupied in current block

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call fatalError(Here,'In block: '// self % current_block_name // &
                          ' Cannot start new array: '//name//          &
                          ' Becouse is writing array: '// self % current_array_name)
    end if

    ! * Check whether shape is degenerate

    ! Load shape information into buffer
    self % shapeBuffer = shape

    ! Load array name
    self % current_array_name = name

    ! Update State -> Array is present but has undefined type
    self % arrayType = UNDEF_ARRAY

    ! Print begining of new entry
    ! We don't know the shape of array becouse if has results its rank will increase
    call self % output % startEntry(name)

  end subroutine startArray

  !!
  !! Finish writing current array
  !!
  subroutine endArray(self)
    class(outputFile), intent(inout) :: self
    character(100), parameter :: Here = 'endArray ( outputFile_class.f90)'

    ! Check that all entries for current array were given
    if(self % arrayTop /= self % arrayLimit) then
      call fatalError(Here, 'Cannot close array: ' // trim(self % current_array_name) // &
                            ' in block: ' // trim(self % current_block_name)         //  &
                            ' Becouse only: '// numToChar(self % arrayTop)          //   &
                            ' of '// numToChar(self % arrayLimit)//' were provided')
    end if

    ! Update state
    self % arrayType  = NOT_ARRAY
    self % arrayLimit = 0
    self % arrayTop   = 0
    self % current_array_name = ''


    ! Close current array and entry
    call self % output % endArray()
    call self % output % endEntry()

  end subroutine endArray

  !!
  !! Add result to the array
  !!
  subroutine addResult_scalar(self,val,std)
    class(outputFile), intent(inout) :: self
    real(defReal), intent(in)        :: val
    real(defReal), intent(in)        :: std
    character(100),parameter :: Here ='addResult_scalar (outputFile_class.f90)'

    ! Check is array is undefined
    if( self % arrayType == UNDEF_ARRAY) then
      ! Update state
      self % arrayType = RES_ARRAY
      self % arrayTop  = 0
      self % arrayLimit = product([2,self % shapeBuffer])

      ! Inform printer about the array
      call self % output % startArray([2,self % shapeBuffer])
    end if

    ! Check for writing into non-existant array, array of mixed elements or array overflow
    if ( self % arrayType == NOT_ARRAY) call fatalError(Here,'Trying to add result without starting array')
    if ( self % arrayType /= RES_ARRAY) call fatalError(Here,'Arrays with mixed content are not allowed')
    if ( self % arrayTop + 2 > self % arrayLimit) then
      call fatalError(Here,'Array overflow. To many elements were provided')
    end if

    ! Print value and std
    call self % output % printNum( trim(self % num2char(val)))
    call self % output % printNum( trim(self % num2char(std)))

    ! Update state
    self % arrayTop = self % arrayTop + 2

  end subroutine addResult_scalar

  !!
  !! Add rank1 array of result to the array
  !!
  subroutine addResult_rank1(self,val,std)
    class(outputFile), intent(inout)        :: self
    real(defReal), dimension(:), intent(in) :: val
    real(defReal), dimension(:),intent(in)  :: std
    integer(shortInt)                       :: N, i
    character(100),parameter :: Here ='addResult_rank1 (outputFile_class.f90)'

    N = size(val)
    if (N /= size(std)) call fatalError(Here,'val and std have diffrent size.')

    ! Add all individual entries
    do i=1,N
      call self % addResult_scalar(val(i),std(i))
    end do

  end subroutine addResult_rank1

  !!
  !! Print complete single result entry
  !!
  subroutine printResult(self,val,std,name)
    class(outputFile), intent(inout)  :: self
    real(defReal),intent(in)          :: val
    real(defReal),intent(in)          :: std
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printResult (outputFile_class.f90)'

    !*** Check that name is unique in current block
    !*** Add name to keywords occupied in current block

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call fatalError(Here,'In block: '// self % current_block_name // &
                          ' Cannot print new entry: '//name//          &
                          ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Print single result as a 2 element array
    call self % output % startEntry(name)
    call self % output % startArray([2])
    call self % output % printNum(trim(self % num2Char(val)))
    call self % output % printNum(trim(self % num2Char(std)))
    call self % output % endArray()
    call self % output % endEntry()

  end subroutine printResult


  !!
  !! Add value of defReal to the array
  !!
  subroutine addValue_defReal_scalar(self,val)
    class(outputFile), intent(inout) :: self
    real(defReal), intent(in)        :: val
    character(100),parameter :: Here ='addValue_defReal_scalar (outputFile_class.f90)'

    ! Check is array is undefined
    if( self % arrayType == UNDEF_ARRAY) then
      ! Update state
      self % arrayType = VAL_ARRAY_REAL
      self % arrayTop  = 0
      self % arrayLimit = product([self % shapeBuffer])

      ! Inform printer about the array
      call self % output % startArray(self % shapeBuffer)
    end if

    ! Check for writing into non-existant array, array of mixed elements or array overflow
    if ( self % arrayType == NOT_ARRAY) then
      call fatalError(Here,'Trying to add result without starting array')

    else if ( self % arrayType /= VAL_ARRAY_REAL) then
      call fatalError(Here,'Arrays with mixed content are not allowed')

    else if ( self % arrayTop + 1 > self % arrayLimit) then
      call fatalError(Here,'Array overflow. To many elements were provided')

    end if

    ! Print value and std
    call self % output % printNum( trim(self % num2char(val)))

    ! Update state
    self % arrayTop = self % arrayTop + 1

  end subroutine addValue_defReal_scalar

  !!
  !! Add rank1 array of defReal to the array
  !!
  subroutine addValue_defReal_rank1(self,val)
    class(outputFile), intent(inout)        :: self
    real(defReal), dimension(:), intent(in) :: val
    integer(shortInt)                       ::  i

    ! Add all individual entries
    do i=1,size(val)
      call self % addValue_defReal_scalar(val(i))
    end do

  end subroutine addValue_defReal_rank1

  !!
  !! Print complete single result entry
  !!
  subroutine printValue_defReal(self,val,name)
    class(outputFile), intent(inout)  :: self
    real(defReal),intent(in)          :: val
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printValue_defReal (outputFile_class.f90)'

    !*** Check that name is unique in current block
    !*** Add name to keywords occupied in current block

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call fatalError(Here,'In block: '// self % current_block_name // &
                          ' Cannot print new entry: '//name//          &
                          ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Print single result as a 2 element array
    call self % output % startEntry(name)
    call self % output % printNum(trim(self % num2Char(val)))
    call self % output % endEntry()

  end subroutine printValue_defReal

  !!
  !! Add value of shortInt to the array
  !!
  subroutine addValue_shortInt_scalar(self,val)
    class(outputFile), intent(inout) :: self
    integer(shortInt), intent(in)    :: val
    integer(longInt)                 :: val_t
    character(100),parameter :: Here ='addValue_shortInt_scalar (outputFile_class.f90)'

    val_t = val
    call self % addValue_longInt_scalar(val_t)

  end subroutine addValue_shortInt_scalar

  !!
  !! Add rank1 array of shortInt to the array
  !!
  subroutine addValue_shortInt_rank1(self,val)
    class(outputFile), intent(inout)            :: self
    integer(shortInt), dimension(:), intent(in) :: val
    integer(shortInt)                           ::  i

    ! Add all individual entries
    do i=1,size(val)
      call self % addValue_shortInt_scalar(val(i))
    end do

  end subroutine addValue_shortInt_rank1

  !!
  !! Print complete single result entry
  !!
  subroutine printValue_shortInt(self,val,name)
    class(outputFile), intent(inout)  :: self
    integer(shortInt), intent(in)     :: val
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printValue_shortInt (outputFile_class.f90)'

    !*** Check that name is unique in current block
    !*** Add name to keywords occupied in current block

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call fatalError(Here,'In block: '// self % current_block_name // &
                          ' Cannot print new entry: '//name//          &
                          ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Print single result as a 2 element array
    call self % output % startEntry(name)
    call self % output % printNum(trim(self % num2Char(val)))
    call self % output % endEntry()

  end subroutine printValue_shortInt


  !!
  !! Add value of longInt to the array
  !!
  subroutine addValue_longInt_scalar(self,val)
    class(outputFile), intent(inout) :: self
    integer(longInt), intent(in)     :: val
    character(100),parameter :: Here ='addValue_longInt_scalar (outputFile_class.f90)'

    ! Check is array is undefined
    if( self % arrayType == UNDEF_ARRAY) then
      ! Update state
      self % arrayType = VAL_ARRAY_INT
      self % arrayTop  = 0
      self % arrayLimit = product([self % shapeBuffer])

      ! Inform printer about the array
      call self % output % startArray(self % shapeBuffer)
    end if

    ! Check for writing into non-existant array, array of mixed elements or array overflow
    if ( self % arrayType == NOT_ARRAY) then
      call fatalError(Here,'Trying to add result without starting array')

    else if ( self % arrayType /= VAL_ARRAY_INT) then
      call fatalError(Here,'Arrays with mixed content are not allowed')

    else if ( self % arrayTop + 1 > self % arrayLimit) then
      call fatalError(Here,'Array overflow. To many elements were provided')

    end if

    ! Print value and std
    call self % output % printNum( trim(self % num2char(val)))

    ! Update state
    self % arrayTop = self % arrayTop + 1

  end subroutine addValue_longInt_scalar

  !!
  !! Add rank1 array of shortInt to the array
  !!
  subroutine addValue_longInt_rank1(self,val)
    class(outputFile), intent(inout)            :: self
    integer(longInt), dimension(:), intent(in)  :: val
    integer(shortInt)                           ::  i

    ! Add all individual entries
    do i=1,size(val)
      call self % addValue_longInt_scalar(val(i))
    end do

  end subroutine addValue_longInt_rank1

  !!
  !! Print complete single result entry
  !!
  subroutine printValue_longInt(self,val,name)
    class(outputFile), intent(inout)  :: self
    integer(longInt), intent(in)      :: val
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printValue_longInt (outputFile_class.f90)'

    !*** Check that name is unique in current block
    !*** Add name to keywords occupied in current block

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call fatalError(Here,'In block: '// self % current_block_name // &
                          ' Cannot print new entry: '//name//          &
                          ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Print single result as a 2 element array
    call self % output % startEntry(name)
    call self % output % printNum(trim(self % num2Char(val)))
    call self % output % endEntry()

  end subroutine printValue_longInt


  !!
  !! Add value of longInt to the array
  !!
  subroutine addValue_char_scalar(self,val)
    class(outputFile), intent(inout) :: self
    character(*), intent(in)         :: val
    character(100),parameter :: Here ='addValue_char_scalar (outputFile_class.f90)'

    ! Check is array is undefined
    if( self % arrayType == UNDEF_ARRAY) then
      ! Update state
      self % arrayType = VAL_ARRAY_CHAR
      self % arrayTop  = 0
      self % arrayLimit = product([self % shapeBuffer])

      ! Inform printer about the array
      call self % output % startArray(self % shapeBuffer)
    end if

    ! Check for writing into non-existant array, array of mixed elements or array overflow
    if ( self % arrayType == NOT_ARRAY) then
      call fatalError(Here,'Trying to add result without starting array')

    else if ( self % arrayType /= VAL_ARRAY_CHAR) then
      call fatalError(Here,'Arrays with mixed content are not allowed')

    else if ( self % arrayTop + 1 > self % arrayLimit) then
      call fatalError(Here,'Array overflow. To many elements were provided')

    end if

    ! Print value and std
    call self % output % printChar(val)

    ! Update state
    self % arrayTop = self % arrayTop + 1

  end subroutine addValue_char_scalar

  !!
  !! Add rank1 array of shortInt to the array
  !!
  subroutine addValue_char_rank1(self,val)
    class(outputFile), intent(inout)            :: self
    character(*), dimension(:), intent(in)      :: val
    integer(shortInt)                           ::  i

    ! Add all individual entries
    do i=1,size(val)
      call self % addValue_char_scalar(val(i))
    end do

  end subroutine addValue_char_rank1

  !!
  !! Print complete single result entry
  !!
  subroutine printValue_char(self,val,name)
    class(outputFile), intent(inout)  :: self
    character(*), intent(in)          :: val
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printValue_char (outputFile_class.f90)'

    !*** Check that name is unique in current block
    !*** Add name to keywords occupied in current block

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call fatalError(Here,'In block: '// self % current_block_name // &
                          ' Cannot print new entry: '//name//          &
                          ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Print single result as a 2 element array
    call self % output % startEntry(name)
    call self % output % printChar(trim(val))
    call self % output % endEntry()

  end subroutine printValue_char

  !!
  !! Print defReal to character of length nameLen using format in self
  !!
  function num2char_defReal(self,r) result(c)
    class(outputFile), intent(in) :: self
    real(defReal), intent(in)     :: r
    character(nameLen)            :: c

    c = ''
    write(c,self % real_format) r
    c = adjustl(c)

  end function num2char_defReal

  !!
  !! Print shortInt to character of length nameLen using format in self
  !!
  function num2char_shortInt(self,i) result(c)
    class(outputFile), intent(in) :: self
    integer(shortInt), intent(in) :: i
    character(nameLen)            :: c

    c = ''
    write(c,self % int_format) i
    c = adjustl(c)

  end function num2char_shortInt

  !!
  !! Print shortInt to character of length nameLen using format in self
  !!
  function num2char_longInt(self,i) result(c)
    class(outputFile), intent(in) :: self
    integer(longInt), intent(in)  :: i
    character(nameLen)            :: c

    c = ''
    write(c,self % int_format) i
    c = adjustl(c)

  end function num2char_longInt
    
end module outputFile_class
