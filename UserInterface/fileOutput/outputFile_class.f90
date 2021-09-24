module outputFile_class

  use numPrecision
  use genericProcedures,       only : fatalError, numToChar
  use stack_class,             only : stackChar
  use charTape_class,          only : charTape
  use charMap_class,           only : charMap
  use asciiOutput_inter,       only : asciiOutput
  use asciiOutputFactory_func, only : new_asciiOutput
  use iso_fortran_env,         only : OUTPUT_UNIT

  implicit none
  private

  character(*),parameter :: DEF_REAL_FORMAT = '(ES12.5)'
  character(*),parameter :: DEF_INT_FORMAT  = '(I12)'

  integer(shortInt), parameter :: NOT_ARRAY      = 0, &
                                  UNDEF_ARRAY    = 1, &
                                  RES_ARRAY      = 2, &
                                  VAL_ARRAY_REAL = 3, &
                                  VAL_ARRAY_INT  = 4, &
                                  VAL_ARRAY_CHAR = 5


  integer(shortInt), parameter :: IS_PRESENT = -1, NOT_PRESENT = -2

  !!
  !! Stack to store the dictionaries with occupied entries
  !!
  !! Note that blocks in outputFile are enumarated from 0
  !! So index in stack is idx = blockLevel + 1
  !!
  !! Imlementation must be rebust to avoid segmentation errors if outut is used
  !! in incorrect sequence (e.g. closing more blocks then were open)
  !!
  !! Interface
  !!   init -> Initialise
  !!   push -> Add new level. Grow stack if needed
  !!   pop  -> Decrease stack size by 1
  !!   add  -> Add new name at a level b_lvl
  !!   isPresent -> True if name is already present at block level b_lvl
  !!   kill -> Return to uninitialised state
  !!
  type, private :: charMapStack
    private
    integer(shortInt)                        :: lvl = 0
    type(charMap), dimension(:), allocatable :: stack
  contains
    procedure :: init => init_charMapStack
    procedure :: push => push_charMapStack
    procedure :: pop  => pop_charMapStack
    procedure :: add  => add_charMapStack
    procedure :: isPresent => isPresent_charMapStack
    procedure :: kill      => kill_charMapStack
  end type charMapStack

  !!
  !! Common interface for output files
  !!
  !! Allows to support multiple output formats while providing a single interface (No need to
  !! change code to alter/add an output format).
  !!
  !! NOTE:
  !!   Is not thread safe. Avoid multiple threads writing to the same output file.
  !!
  !! Output file:
  !!  -> Consists of multiple nested block. Each block has a name with exception of the 1st, which
  !!       is called root
  !!  -> Each entry in a block (including sub-blocks) must have a unique name
  !!  -> Each entry is either a raw "Name - Value" pair or an array
  !!  -> Arrays can have unlimited rank (dimension), but values of the enteries must be given in a
  !!     sequence in a COLUMN-MAJOR ORDER
  !!  -> Arrays can store RESULTS, REALS, INTS or CHARACTERS.
  !!  -> RESULT is a pair of reals. MEAN and STANDARD DEVIATION
  !!
  !! Writing output:
  !!  Is done by the a sequence of calls, which open block, start arrays, add values etc.
  !!  OutputFile checks that the order and number of calls is correct. It also is responsible for
  !!  handling errors.
  !!
  !! Handling errors:
  !!  If `fatalError` is true (default) then on wrong sequence of calls fatalError will be raised.
  !!  Otherwise the flag `noErrors` will be set to false and error massege saved in the `errorLog`
  !!
  !! Repeated Names:
  !!  If name in a block is repeted bahviour is governed by `logNameReuse` function. If `fatalError`
  !!  is true, the nwarning message is produced. Otherwise error is logged in `errorLog` and
  !!  `noErrors` flag set to FALSE.
  !!
  !! Private Memebers:
  !!   output      -> Output printer that determines the format
  !!   type        -> Type of the output printer
  !!   fatalErrors -> True (default) makes all errors fatal
  !!   noErrors    -> True if no errors were encountered since initialisation
  !!   errorLog    -> Log of all errors
  !!   blockLevel  -> Current block level (0 is root)
  !!   arrayTop    -> Counter of the number of elements given to the array since it start
  !!   arrayLimit  -> The number of elements that needs to be given to the array
  !!   arrayType   -> Type of the array currently written
  !!   current_block_name -> Name of the current block
  !!   current_block_name -> Name of the current array
  !!   block_name_stack   -> Stack of previous block names
  !!   shapeBuffer -> Saved shape of the current array
  !!   real_format -> Format to convert real to char
  !!   int_format -> Format to convert int to char
  !!
  !! Interface:
  !!   init           -> Initialise the output file with format
  !!   writeToFile    -> Write the output to a file
  !!   writeToConsole -> Write the output to the console
  !!   reset          -> Discard already written output.
  !!   isValid        -> Return true if no errors were raised since initialisation
  !!   getErrorLog    -> Get the error log
  !!   startBlock     -> Start new block in the output
  !!   endBlock       -> End the current block
  !!   startArray     -> Start new array in the current block
  !!   endArray       -> End array in the current block
  !!   addValue       -> Add value REAL, INT or CHAR to the array
  !!   addResult      -> Add MEAN-STANDARD DEVIATION pair to the array
  !!   printValue     -> Add a NAME-VALUE pair to the current block
  !!   printResult    -> Add a NAME - RESULT pair to the current block
  !!   num2char       -> Convert a number to character
  !!
  type, public :: outputFile
    private
    class(asciiOutput),allocatable :: output
    character(nameLen)             :: type

    ! Error handling settings
    logical(defBool) :: fatalErrors = .true. ! Enable fatalErrors on wrong logic
    logical(defBool) :: noErrors    = .true. ! No errors were encountered
    type(charTape)   :: errorLog             ! error messages separated by new line

    ! State variables
    integer(shortInt)  :: blockLevel    = 0  ! Current depth in nested blocks
    integer(shortInt)  :: arrayTop      = 0  ! Number of elements written to array
    integer(shortInt)  :: arrayLimit    = 0  ! Maximum number of elements that can be written to array
    integer(shortInt)  :: arrayType     = NOT_ARRAY
    character(nameLen) :: current_block_name = ''
    character(nameLen) :: current_array_name = ''
    type(stackChar)    :: block_name_stack

    ! Buffors
    integer(shortInt),dimension(:), allocatable :: shapeBuffer
    type(charMapStack)                          :: usedNames

    ! Formats
    character(nameLen) :: real_format = DEF_REAL_FORMAT
    character(nameLen) :: int_format  = DEF_INT_FORMAT

  contains
    procedure :: init
    procedure :: writeToFile
    procedure :: writeToConsole
    procedure :: reset

    ! Error Handling
    procedure, private :: logError
    procedure, private :: logNameReuse
    procedure          :: isValid
    procedure          :: getErrorLog

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
  !! Initialise output file
  !!
  !! Args:
  !!   type [in] -> Character to specify the format of the output file
  !!   fatalErrors [in] -> Optional. Set on/off fatalErrors oninvalid use of output file.
  !!                         Default is ON.
  !!
  subroutine init(self, type, fatalErrors)
    class(outputFile), intent(inout)     :: self
    character(*),intent(in)              :: type
    logical(defBool),optional,intent(in) :: fatalErrors

    self % type = type
    allocate( self % output, source = new_asciiOutput(self % type))

    ! Initialise name stack
    call self % usedNames % init()
    call self % usedNames % push()

    ! Change fatalErrors setting
    if(present(fatalErrors)) then
      self % fatalErrors = fatalErrors
    end if

  end subroutine init

  !!
  !! Dump collected results to file
  !!
  !! Writes the output stored in memory to a file
  !!
  !! Args:
  !!   file [in] -> Path to the output file. Must be given WITHOUT extension
  !!     Approperiate extension will recived from the printer
  !!
  !! Errors:
  !!   fatalError if there is a problem opening the file.
  !!
  subroutine writeToFile(self, file)
    class(outPutFile), intent(inout) :: self
    character(*), intent(in)         :: file
    integer(shortInt)                :: unitNum
    integer(shortInt)                :: error
    character(:), allocatable        :: file_loc
    character(99)                    :: errorMsg
    character(100), parameter :: Here = 'writeToFile (outputFile_class.f90)'

    file_loc = trim(file) // "." // self % output % extension()

    ! Open file to write
    open ( newunit   = unitNum, &
           file      = file_loc, &
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
  !! Print output file to the console (default output)
  !!
  !! Args:
  !!   None
  !!
  subroutine writeToConsole(self)
    class(outputFile), intent(inout) :: self

    call self % output % writeToFile(OUTPUT_UNIT)

  end subroutine writeToConsole

  !!
  !! Discards all results and sets output file to an initialised state
  !!
  subroutine reset(self)
    class(outputFile), intent(inout) :: self

    ! Reallocate ascii output
    deallocate(self % output)
    allocate( self % output, source = new_asciiOutput(self % type))

    ! Clean error log an reset flag
    call self % errorLog % clean()
    self % noErrors = .true.

    ! Reset initial state
    self % blockLevel = 0
    self % arrayTop   = 0
    self % arrayLimit = 0
    self % arrayType  = NOT_ARRAY
    self % current_block_name = ''
    self % current_array_name = ''

    call self % usedNames % kill()
    call self % usedNames % init()
    call self % usedNames % push()

  end subroutine reset

  !!
  !! Handle an error.
  !!
  !! Raises a fatalError is fatalError is set TRUE. Else writes the error message to the log file
  !! and changes the error flag to true.
  !!
  !! Args:
  !!   where [in] -> Character. Location of the error
  !!   what [in] -> Character. Description of the error.
  !!
  subroutine logError(self, where, what)
    class(outputFile), intent(inout) :: self
    character(*), intent(in)         :: where
    character(*), intent(in)         :: what

    if(self % fatalErrors) then ! Kill ran with fatalError
      call fatalError(where, what)

    else ! Set noErrors to false and log error
      self % noErrors = .false.
      call self % errorLog % append(where // ' ; ' // what // new_line(''))

    end if

  end subroutine logError

  !!
  !! Deal with an event when entry name is resued in a block
  !!
  !! Currently is fatalError is TRUE -> Prints warning to the screen
  !! If fatalError is False markes error flagto true and appends a log
  !!
  !! Args:
  !!  name [in] -> Name that has been reused
  !!
  subroutine logNameReuse(self, name)
    class(outputFile), intent(inout) :: self
    character(nameLen), intent(in)   :: name
    character(:), allocatable        :: msg

    msg = "WARNING: Name " // trim(name) // " is used multiple times in block " // &
           trim(self % current_block_name) // " || Output may not parse correctly!"

    if (self % fatalErrors) then ! Print to screen
      print *, msg

    else ! log error
      self % noErrors = .false.
      call self % errorLog % append(msg)
    end if

  end subroutine logNameReuse

  !!
  !! Returns .true. if there were no errors. .false. otherwise
  !!
  !! Args:
  !!   None
  !!
  pure function isValid(self) result(isIt)
    class(outputFile), intent(in) :: self
    logical(defBool)              :: isIt

    isIt = self % noErrors

  end function isValid

  !!
  !! Returns error log as a allocatable character
  !!
  !! Args:
  !!   None
  !!
  function getErrorLog(self) result(errorLog)
    class(outputFile), intent(in) :: self
    character(:),allocatable      :: errorLog

    errorLog = self % errorLog % expose()

  end function getErrorLog

  !!
  !! Start new nested block
  !!
  !! Args:
  !!  name [in] -> Name of the new block (must be unique)
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called when writting an array
  !!   Calls logNameReuse if name is not unique
  !!
  subroutine startBlock(self, name)
    class(outputFile), intent(inout) :: self
    character(nameLen), intent(in)   :: name
    character(100), parameter :: Here ='startBlock (outputFile_class.f90)'

    ! Check that name is unique in current block
    if (self % usedNames % isPresent(name, self % blockLevel)) call self % logNameReuse(name)

    ! Add name to keywords occupied in current block
    call self % usedNames % push()
    call self % usedNames % add(name, self % blockLevel)

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call self % logError(Here,'In block: '// self % current_block_name // &
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
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called when writting an array
  !!     -Tries to end the root block.
  !!
  subroutine endBlock(self)
    class(outputFile), intent(inout) :: self
    character(100), parameter :: Here ='endBlock (outputFile_class.f90)'

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call self % logError(Here,'In block: '// self % current_block_name // &
                                ' Cannot exit from block: '//                &
                                ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Check that is not in root block
    if ( self % blockLevel == 0 ) call self % logError(Here,'Trying to exit from root block')

    ! Update state. Protect against blockLevel == 0
    if ( self % blockLevel > 0) then
      self % blockLevel = self % blockLevel - 1
      call self % block_name_stack % pop( self % current_block_name)
    end if

    ! End block in output
    call self % output % endBlock()

    ! Decrease name stack
    call self % usedNames % pop()

  end subroutine endBlock

  !!
  !! Prepare for writing array. Store shape information
  !!
  !! We don't know the type of content yet so we do not inform printer about shape
  !!
  !! Args:
  !!   name [in] -> Name of the array
  !!   shape [in] -> Shape of the array e.g. [N,M] or [1,3,4]
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called when writting another array
  !!     -Shape is invalid
  !!   Calls logNameReuse if name is not unique
  !!
  subroutine startArray(self, name, shape)
    class(outputFile), intent(inout) :: self
    character(nameLen), intent(in)   :: name
    integer(shortInt),dimension(:)   :: shape
    character(100), parameter :: Here ='startArray (outputFile_class.f90)'

    ! Check that name is unique in current block
    if (self % usedNames % isPresent(name, self % blockLevel)) call self % logNameReuse(name)

    ! Add name to keywords occupied in current block
    call self % usedNames % add(name, self % blockLevel)

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call self % logError(Here,'In block: '// self % current_block_name // &
                                ' Cannot start new array: '//name//          &
                                ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Check whether shape is degenerate
    if (size(shape) == 0 .or. product(shape) == 0) then
      call self % logError(Here,'In block: '// self % current_block_name // &
                                ' Cannot start new array: '//name//          &
                                ' Becouse the shape is degenerate')
    end if

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
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Array is closed before all elements are provided
  !!
  subroutine endArray(self)
    class(outputFile), intent(inout) :: self
    character(100), parameter :: Here = 'endArray ( outputFile_class.f90)'

    ! Check that all entries for current array were given
    if(self % arrayTop /= self % arrayLimit) then
      call self % logError(Here, 'Cannot close array: ' // trim(self % current_array_name) // &
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
  !! Args:
  !!   val [in] -> Mean value
  !!   std [in] -> Associated standard deviation
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addResult_scalar(self, val, std)
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
    if ( self % arrayType == NOT_ARRAY) call self % logError(Here,'Trying to add result without starting array')
    if ( self % arrayType /= RES_ARRAY) call self % logError(Here,'Arrays with mixed content are not allowed')
    if ( self % arrayTop + 2 > self % arrayLimit) then
      call self % logError(Here,'Array overflow. To many elements were provided')
    end if

    ! Print value and std
    call self % output % printNum( trim(self % num2char(val)))
    call self % output % printNum( trim(self % num2char(std)))

    ! Update state
    self % arrayTop = self % arrayTop + 2

  end subroutine addResult_scalar

  !!
  !! Add results to the array
  !!
  !! Args:
  !!   val [in] -> Array of mean valuse
  !!   std [in] -> Associated standard deviations
  !!
  !! Erorrs:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addResult_rank1(self, val, std)
    class(outputFile), intent(inout)         :: self
    real(defReal), dimension(:), intent(in)  :: val
    real(defReal), dimension(:), intent(in)  :: std
    integer(shortInt)                        :: N, i
    character(100),parameter :: Here ='addResult_rank1 (outputFile_class.f90)'

    N = size(val)
    if (N /= size(std)) call self % logError(Here,'val and std have diffrent size.')

    ! Add all individual entries
    do i = 1, N
      call self % addResult_scalar(val(i),std(i))
    end do

  end subroutine addResult_rank1

  !!
  !! Print complete single result entry
  !!
  !! Args:
  !!   val [in] -> Mean value
  !!   std [in] -> Associated standard deviation
  !!   name [in] -> Unique name of this results (within current block)
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called when an array is open
  !!   Calls logNameReuse if name is not unique
  !!
  subroutine printResult(self, val, std, name)
    class(outputFile), intent(inout)  :: self
    real(defReal),intent(in)          :: val
    real(defReal),intent(in)          :: std
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printResult (outputFile_class.f90)'

    ! Check that name is unique in current block
    if (self % usedNames % isPresent(name, self % blockLevel)) call self % logNameReuse(name)

    ! Add name to keywords occupied in current block
    call self % usedNames % add(name, self % blockLevel)

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call self % logError(Here,'In block: '// self % current_block_name // &
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
  !! Add real value to the array
  !!
  !! Args:
  !!   val [in] -> the value
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addValue_defReal_scalar(self, val)
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
      call self % logError(Here,'Trying to add result without starting array')

    else if ( self % arrayType /= VAL_ARRAY_REAL) then
      call self % logError(Here,'Arrays with mixed content are not allowed')

    else if ( self % arrayTop + 1 > self % arrayLimit) then
      call self % logError(Here,'Array overflow. To many elements were provided')

    end if

    ! Print value and std
    call self % output % printNum( trim(self % num2char(val)))

    ! Update state
    self % arrayTop = self % arrayTop + 1

  end subroutine addValue_defReal_scalar

  !!
  !! Add real values to the array
  !!
  !! Args:
  !!   val [in] -> An array of values
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addValue_defReal_rank1(self, val)
    class(outputFile), intent(inout)        :: self
    real(defReal), dimension(:), intent(in) :: val
    integer(shortInt)                       ::  i

    ! Add all individual entries
    do i=1,size(val)
      call self % addValue_defReal_scalar(val(i))
    end do

  end subroutine addValue_defReal_rank1

  !!
  !! Print real value to the block
  !!
  !! Args:
  !!  val [in] -> The value
  !!  name [in] -> Unique name of the value in the block
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called when an array is open
  !!     -Name is not unique in the block
  !!
  subroutine printValue_defReal(self, val, name)
    class(outputFile), intent(inout)  :: self
    real(defReal),intent(in)          :: val
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printValue_defReal (outputFile_class.f90)'

    ! Check that name is unique in current block
    if (self % usedNames % isPresent(name, self % blockLevel)) call self % logNameReuse(name)

    ! Add name to keywords occupied in current block
    call self % usedNames % add(name, self % blockLevel)

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call self % logError(Here,'In block: '// self % current_block_name // &
                               ' Cannot print new entry: '//name//          &
                               ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Print single result as a 2 element array
    call self % output % startEntry(name)
    call self % output % printNum(trim(self % num2Char(val)))
    call self % output % endEntry()

  end subroutine printValue_defReal

  !!
  !! Add short int value to the array
  !!
  !! Args:
  !!   val [in] -> the value
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addValue_shortInt_scalar(self, val)
    class(outputFile), intent(inout) :: self
    integer(shortInt), intent(in)    :: val
    integer(longInt)                 :: val_t
    character(100),parameter :: Here ='addValue_shortInt_scalar (outputFile_class.f90)'

    val_t = val
    call self % addValue_longInt_scalar(val_t)

  end subroutine addValue_shortInt_scalar

  !!
  !! Add short int values to the array
  !!
  !! Args:
  !!   val [in] -> An array of values
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addValue_shortInt_rank1(self, val)
    class(outputFile), intent(inout)            :: self
    integer(shortInt), dimension(:), intent(in) :: val
    integer(shortInt)                           ::  i

    ! Add all individual entries
    do i=1,size(val)
      call self % addValue_shortInt_scalar(val(i))
    end do

  end subroutine addValue_shortInt_rank1

  !!
  !! Print short int value to the block
  !!
  !! Args:
  !!  val [in] -> The value
  !!  name [in] -> Unique name of the value in the block
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called when an array is open
  !!   Calls logNameReuse if name is not unique
  !!
  subroutine printValue_shortInt(self, val, name)
    class(outputFile), intent(inout)  :: self
    integer(shortInt), intent(in)     :: val
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printValue_shortInt (outputFile_class.f90)'

    ! Check that name is unique in current block
    if (self % usedNames % isPresent(name, self % blockLevel)) call self % logNameReuse(name)

    ! Add name to keywords occupied in current block
    call self % usedNames % add(name, self % blockLevel)

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call self % logError(Here,'In block: '// self % current_block_name // &
                                ' Cannot print new entry: '//name//          &
                                ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Print single result as a 2 element array
    call self % output % startEntry(name)
    call self % output % printNum(trim(self % num2Char(val)))
    call self % output % endEntry()

  end subroutine printValue_shortInt

  !!
  !! Add long int value to the array
  !!
  !! Args:
  !!   val [in] -> the value
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addValue_longInt_scalar(self, val)
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
      call self % logError(Here,'Trying to add result without starting array')

    else if ( self % arrayType /= VAL_ARRAY_INT) then
      call self % logError(Here,'Arrays with mixed content are not allowed')

    else if ( self % arrayTop + 1 > self % arrayLimit) then
      call self % logError(Here,'Array overflow. To many elements were provided')

    end if

    ! Print value and std
    call self % output % printNum( trim(self % num2char(val)))

    ! Update state
    self % arrayTop = self % arrayTop + 1

  end subroutine addValue_longInt_scalar

  !!
  !! Add long int values to the array
  !!
  !! Args:
  !!   val [in] -> An array of values
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addValue_longInt_rank1(self, val)
    class(outputFile), intent(inout)            :: self
    integer(longInt), dimension(:), intent(in)  :: val
    integer(shortInt)                           ::  i

    ! Add all individual entries
    do i=1,size(val)
      call self % addValue_longInt_scalar(val(i))
    end do

  end subroutine addValue_longInt_rank1

  !!
  !! Print long int value to the block
  !!
  !! Args:
  !!  val [in] -> The value
  !!  name [in] -> Unique name of the value in the block
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called when an array is open
  !!   Calls logNameReuse if name is not unique
  !!
  subroutine printValue_longInt(self, val, name)
    class(outputFile), intent(inout)  :: self
    integer(longInt), intent(in)      :: val
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printValue_longInt (outputFile_class.f90)'

    ! Check that name is unique in current block
    if (self % usedNames % isPresent(name, self % blockLevel)) call self % logNameReuse(name)

    ! Add name to keywords occupied in current block
    call self % usedNames % add(name, self % blockLevel)

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call self % logError(Here,'In block: '// self % current_block_name // &
                                ' Cannot print new entry: '//name//          &
                                ' Becouse is writing array: '// self % current_array_name)
    end if

    ! Print single result as a 2 element array
    call self % output % startEntry(name)
    call self % output % printNum(trim(self % num2Char(val)))
    call self % output % endEntry()

  end subroutine printValue_longInt

  !!
  !! Add character value to the array
  !!
  !! Args:
  !!   val [in] -> the value
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addValue_char_scalar(self, val)
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
      call self % logError(Here,'Trying to add result without starting array')

    else if ( self % arrayType /= VAL_ARRAY_CHAR) then
      call self % logError(Here,'Arrays with mixed content are not allowed')

    else if ( self % arrayTop + 1 > self % arrayLimit) then
      call self % logError(Here,'Array overflow. To many elements were provided')

    end if

    ! Print value and std
    call self % output % printChar(trim(val))

    ! Update state
    self % arrayTop = self % arrayTop + 1

  end subroutine addValue_char_scalar

  !!
  !! Add character values to the array
  !!
  !! Args:
  !!   val [in] -> An array of values
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called before an array is opened
  !!     -Currently opened array has different type
  !!     -To many entries are fed into the opened array
  !!
  subroutine addValue_char_rank1(self, val)
    class(outputFile), intent(inout)            :: self
    character(*), dimension(:), intent(in)      :: val
    integer(shortInt)                           ::  i

    ! Add all individual entries
    do i = 1, size(val)
      call self % addValue_char_scalar(val(i))
    end do

  end subroutine addValue_char_rank1

  !!
  !! Print character value to the block
  !!
  !! Args:
  !!  val [in] -> The value
  !!  name [in] -> Unique name of the value in the block
  !!
  !! Errors:
  !!   Calls logError if:
  !!     -Called when an array is open
  !!   Calls logNameReuse if name is not unique
  !!
  subroutine printValue_char(self, val, name)
    class(outputFile), intent(inout)  :: self
    character(*), intent(in)          :: val
    character(nameLen), intent(in)    :: name
    character(100), parameter :: Here = 'printValue_char (outputFile_class.f90)'

    ! Check that name is unique in current block
    if (self % usedNames % isPresent(name, self % blockLevel)) call self % logNameReuse(name)

    ! Add name to keywords occupied in current block
    call self % usedNames % add(name, self % blockLevel)

    ! Check that currently is not writing array
    if ( self % arrayType /= NOT_ARRAY) then
      call self % logError(Here,'In block: '// self % current_block_name // &
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
  !! Args:
  !!   r [in] -> Real value
  !!
  function num2char_defReal(self, r) result(c)
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
  !! Args:
  !!   i [in] -> Value of the integer
  !!
  function num2char_shortInt(self, i) result(c)
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
  !! Args:
  !!   i [in] ->Value of the integer
  !!
  function num2char_longInt(self, i) result(c)
    class(outputFile), intent(in) :: self
    integer(longInt), intent(in)  :: i
    character(nameLen)            :: c

    c = ''
    write(c,self % int_format) i
    c = adjustl(c)

  end function num2char_longInt

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! charMapStack Procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Initialise charMapStack to a initial size
  !!
  !! Args:
  !!   None
  !!
  subroutine init_charMapStack(self)
    class(charMapStack), intent(inout) :: self

    ! In case of double initialisation
    call self % kill()

    allocate(self % stack(5))
    self % lvl = 0

  end subroutine init_charMapStack

  !!
  !! Add extra level to the stack
  !!
  !! Args:
  !!  None
  !!
  subroutine push_charMapStack(self)
    class(charMapStack), intent(inout)       :: self
    type(charMap), dimension(:), allocatable :: temp
    integer(shortInt)                        :: i

    ! Increment level
    self % lvl = self % lvl + 1

    ! If run out of space grow the stack
    if (self % lvl > size(self % stack)) then
      allocate(temp(size(self % stack) * 2))
      do i = 1,size(self % stack)
        temp(i) = self % stack(i)
      end do
      call move_alloc(temp, self % stack)
    end if

  end subroutine push_charMapStack

  !!
  !! Remove top level from the stack
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   Poping empty stack does nothing!
  !!
  subroutine pop_charMapStack(self)
    class(charMapStack), intent(inout) :: self

    ! Clean dictionary
    if (self % lvl > 0) then
      call self % stack(self % lvl) % kill()
      ! Decrement level
      self % lvl = self % lvl - 1
    end if

  end subroutine pop_charMapStack

  !!
  !! Add new entery at the given block level
  !!
  !! Args:
  !!   name [in]  -> Entry name
  !!   b_lvl [in] -> Block level (idx = b_lvl + 1)
  !!
  !! Errors:
  !!   Does nothing if b_lvl is out of range
  !!
  subroutine add_charMapStack(self, name, b_lvl)
    class(charMapStack), intent(inout) :: self
    character(nameLen), intent(in)     :: name
    integer(shortInt), intent(in)      :: b_lvl

    if (b_lvl >= 0 .and. b_lvl < self % lvl) then
      call self % stack(b_lvl + 1) % add(name, IS_PRESENT)
    end if

  end subroutine add_charMapStack

  !!
  !! Check if the name is present at the given block lvl
  !!
  !! Args:
  !!   name [in]  -> Entry name
  !!   b_lvl [in] -> Block level (idx = b_lvl + 1)
  !!
  !! Result:
  !!  True if it is. False otherwsie
  !!
  !! Errors:
  !!   Returns FALSE is b_lvl is out of range
  !!
  function isPresent_charMapStack(self, name, b_lvl) result(isIt)
    class(charMapStack), intent(inout) :: self
    character(nameLen), intent(in)     :: name
    integer(shortInt), intent(in)      :: b_lvl
    logical(defBool)                   :: isIt

    isIt = b_lvl >= 0 .and. b_lvl < self % lvl

    if (isIt) then
      isIt = IS_PRESENT == self % stack(b_lvl + 1) % getOrDefault(name, NOT_PRESENT)
    end if

  end function isPresent_charMapStack

  !!
  !! Return to uninitialised state
  !!
  subroutine kill_charMapStack(self)
    class(charMapStack), intent(inout) :: self

    if(allocated(self % stack)) deallocate(self % stack)
    self % lvl = 0

  end subroutine kill_charMapStack

end module outputFile_class
