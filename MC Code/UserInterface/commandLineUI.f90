module commandLineUI

  use numPrecision
  use genericProcedures,  only : fatalError, numToChar, targetNotFound
  use hashFunctions_func, only : FNV_1

  implicit none
  private

  !!
  !! This module is used to define and parse command line options
  !! Command line is parsed during first execution of getInputFile/clOptionPresent/getFromCL
  !! so all command line options need to be defined with addClOption before the first enquiry.
  !!
  !! There is a special command line input for an Input file, which has no keyword.
  !! The first command line argument that does not match any defined keyword is interpreted to be
  !! an input file!
  !!
  !! If the program is called with one of the special help keywords from HELP_KEYWORDS parameter
  !! anywhere in the command line input, parser will stop execution and display all defined options.
  !!
  !! Following procedures are avalible to read command line options
  !!  -> getInputFile(char) - Subroutine which loads input file into allocatable char
  !!  -> getFromCl(var,keyword,opt) - Subtroutine that loads opt-th argument of option under keyword
  !!                                  into var. Var must be shortInt, defReal or allocatable char
  !!  -> clOptionIsPresent(keyword) - Returns .true. if the option uder keyword was invoked
  !!
  public :: addClOption
  public :: getInputFile
  public :: clOptionIsPresent
  public :: getFromCL

  interface getFromCL
    module procedure getFromCL_shortInt
    module procedure getFromCL_defReal
    module procedure getFromCl_char
  end interface

  !!
  !! Module parameters
  !!
  integer(shortInt), parameter :: UNDEFINED_TYPE = 0, &
                                  INT_TYPE       = 1, &
                                  REAL_TYPE      = 2, &
                                  CHAR_TYPE      = 3

  character(*), dimension(2),parameter :: HELP_KEYWORDS = ['-help','-h   ']

  !!
  !! Type to hold a single option argument (type + character with value)
  !!
  type, private :: optionArg
    integer(shortInt)        :: type = UNDEFINED_TYPE
    character(:),allocatable :: val
  contains
    procedure :: typeChar => typeChar_optionArg
    procedure :: verify   => verify_optionArg
  end type optionArg

  !!
  !! Type that holds definition and contents of an option
  !! NOTE: Options are unique. Only one option of a given type is allowed
  !!
  type, private :: optionDescriptor
    logical(defBool)                         :: isLoaded
    integer(shortInt)                        :: keywordHash
    type(optionArg),dimension(:),allocatable :: optArgs
    character(:),allocatable                 :: keyword
    character(:),allocatable                 :: help
  contains
    procedure :: init         => init_optionDescriptor
    procedure :: read         => read_optionDescriptor
    generic   :: operator(==) => eq_optionDescriptor
    procedure :: display      => display_optionDescriptor

    procedure,private :: eq_optionDescriptor
  end type optionDescriptor

  !!
  !! Module components
  !!
  logical(defBool)                                  :: parsed = .false.
  character(:), allocatable                         :: input
  type(optionDescriptor), dimension(:), allocatable :: options

contains

  !!
  !! NOTE: Ideally this should be preprocessor macro resolved at compile-time
  !!       For now it is left as a subroutine revolved at runtime
  !!
  !! Subroutine which declares a new command line option:
  !!  keyword  -> Option keyword, should begin with '-' or '--' but doesn't have to
  !!  Narg     -> Number of arguments that follows keyword
  !!  argTypes -> type of argument ('int', 'real', 'char'). Must have size of Narg
  !!  help     -> help string displayed when help is invoked
  !!
  subroutine addClOption(keyword,Narg,argTypes,help)
    character(*), intent(in)                        :: keyword
    integer(shortInt), intent(in)                   :: Narg
    character(*),dimension(:),intent(in)            :: argTypes
    character(*),intent(in)                         :: help
    type(optionDescriptor),dimension(:),allocatable :: temp
    integer(shortInt)                               :: N
    character(100), parameter :: Here = 'addClOption (commandLinieUI.f90)'

    ! Extend or allocate options array
    if (allocated(options)) then
      ! Find new size
      N = size(options)
      N = N + 1

      ! Create new extended array
      allocate(temp(N))
      temp(1:N-1) = options

      ! Move allocation
      call move_alloc(temp, options)
    else
      ! Allocate memory
      N = 1
      allocate(options(N))
    end if

    ! Initialise new entery
    call options(N) % init(keyword,Narg,argTypes,help)

    ! Check that new entry is unique
    if (any( options(1:N-1) == options(N))) then
      call parseError('Option with keyword: '//trim(options(N) % keyword)//' is already present')
    end if

  end subroutine addClOption

  !!
  !! Returns true if a given option is present
  !!
  function clOptionIsPresent(keyword) result(isIt)
    character(*),intent(in) :: keyword
    logical(defBool)        :: isIt
    integer(shortInt)       :: idx

    if(.not.parsed) call parseCL()

    idx = findOptionIdx(trim(adjustl(keyword)))

    isIt = idx /= targetNotFound

    if(isIt) isIt = options(idx) % isLoaded

  end function clOptionIsPresent

  !!
  !! Read path to input file into string of fixed length
  !! Returns error if input path is not present
  !!
  subroutine getInputFile(string)
    character(:),allocatable, intent(out) :: string
    character(100), parameter :: Here= 'getInputFile (commandLineUI.f90)'

    if(.not.parsed) call parseCL()

    if(.not.allocated(input)) call fatalError(Here, 'No input file was provided')

    ! Fill with blanks and load input
    string =input

  end subroutine getInputFile

  !!
  !! Parse command line arguments
  !!
  subroutine parseCL()
    integer(shortInt)         :: argCount
    integer(shortInt)         :: i, keywords
    integer(shortInt)         :: argLen
    integer(shortInt)         :: idx
    character(:), allocatable :: string
    character(100), parameter :: Here = 'parseCL (commandLinieUI.f90)'

    argCount = command_argument_count()

    ! There are no command arguments
    if( argCount == 0) then
      parsed = .true.
      return
    end if

    ! Read arguments
    i = 1
    do while (i <= argCount)
      ! Read length of argument
      call get_command_argument(i, length = argLen)

      ! Allocate string for a command line argument
      if(allocated(string)) deallocate(string)
      allocate(character(argLen) :: string)

      ! Get command argument
      call get_command_argument(i,string)
      string = trim(adjustl(string))

      ! Special case "-help" -> display options and stop execution
      if( any(HELP_KEYWORDS == string)) then
        call displayCloptions
        stop
      end if

      ! Look for command
      idx = findOptionIdx(string)

      if (idx == targetNotFound) then
        ! If keyword was not found there are two options
        ! Either it is input (if input wasn't already read) or unrecognised keyword
        if( .not.allocated(input)) then
          input = string
        else
          call parseError("Option: "//trim(string)//" is not defined")
        end if

      else
        call options(idx) % read(i,argCount)

      end if
      ! Increment position
      i = i +1
    end do
    ! Set to parsed
    parsed = .true.

  end subroutine parseCL

  !!
  !! Search for the keyword among defined options
  !!  Returns targetNotFound if there is no option under keyword
  !!
  function findOptionIdx(keyword) result(idx)
    character(*),intent(in) :: keyword
    integer(shortInt)       :: idx
    integer(shortInt)       :: hash, i
    logical(defBool)        :: same

    ! Hash the keyword
    call FNV_1(trim(adjustl(keyword)),hash)

    ! Loop over all options
    do i=1,size(options)
      ! Compare hashes
      same = hash == options(i) % keywordHash

      ! Compare keywords
      if(same) same = trim(adjustl(keyword)) == options(i) % keyword

      if(same) then
        idx = i
        return
      end if
    end do

    ! Keyword was not found
    idx = targetNotFound

  end function findOptionIdx

  !!
  !! Obtain integer argument from option under "keyword"
  !! Take opt-th argument
  !!
  subroutine getFromCL_shortInt(var,keyword,opt)
    integer(shortInt), intent(out) :: var
    character(*),intent(in)        :: keyword
    integer(shortInt),intent(in)   :: opt
    integer(shortInt)              :: idx
    character(100), parameter :: Here = 'getFromCL_shortInt (commandLineUI.f90)'

    if(.not.parsed) call parseCL()

    ! Find index of the option
    idx = findOptionIdx(keyword)

    ! Check if option is defined
    if( idx == targetNotFound) then
      call fatalError(Here,'Option with a keyword: '//trim(adjustl(keyword))//' was not defined.')
    end if

    ! Check if option was invoked (has values)
    if (.not.options(idx) % isLoaded) then
      call fatalError(Here,'Option with a keyword: '//trim(adjustl(keyword))//' was not invoked.')
    end if

    ! Check that value of opt checks out
    if ( opt < 1) call fatalError(Here,'Option argument index must be < 1')
    if ( opt > size(options(idx) % optArgs) ) then
      call fatalError(Here,'Option '//trim(adjustl(keyword))// 'has only'// &
                            numToChar(size(options(idx) % optArgs))//' input arguments')
    end if

    ! Check that type of option argument matches
    if (options(idx) % optArgs(opt) % type /= INT_TYPE) then
      call fatalError(Here,numToChar(opt)//'-th argument of option '//trim(adjustl(keyword))//&
                           'is an integer')
    end if

    ! Read the value
    read(options(idx) % optArgs(opt) % val, '(I30)' ) var

  end subroutine getFromCL_shortInt

  !!
  !! Obtain real argument from option under "keyword"
  !! Take opt-th argument
  !!
  subroutine getFromCL_defReal(var,keyword,opt)
    real(defReal), intent(out)     :: var
    character(*),intent(in)        :: keyword
    integer(shortInt),intent(in)   :: opt
    integer(shortInt)              :: idx
    character(100), parameter :: Here = 'getFromCL_defReal (commandLineUI.f90)'

    if(.not.parsed) call parseCL()

    ! Find index of the option
    idx = findOptionIdx(keyword)

    ! Check if option is defined
    if( idx == targetNotFound) then
      call fatalError(Here,'Option with a keyword: '//trim(adjustl(keyword))//' was not defined.')
    end if

    ! Check if option was invoked (has values)
    if (.not.options(idx) % isLoaded) then
      call fatalError(Here,'Option with a keyword: '//trim(adjustl(keyword))//' was not invoked.')
    end if

    ! Check that value of opt checks out
    if ( opt < 1) call fatalError(Here,'Option argument index must be < 1')
    if ( opt > size(options(idx) % optArgs) ) then
      call fatalError(Here,'Option '//trim(adjustl(keyword))// 'has only'// &
                            numToChar(size(options(idx) % optArgs))//' input arguments')
    end if

    ! Check that type of option argument matches
    if (options(idx) % optArgs(opt) % type /= REAL_TYPE) then
      call fatalError(Here,numToChar(opt)//'-th argument of option '//trim(adjustl(keyword))//&
                           'is a real number')
    end if

    ! Read the value
    read(options(idx) % optArgs(opt) % val, '(G30.30)' ) var

  end subroutine getFromCL_defReal

  !!
  !! Obtain real argument from option under "keyword"
  !! Take opt-th argument
  !!
  subroutine getFromCL_char(var,keyword,opt)
    character(:),allocatable, intent(out)  :: var
    character(*),intent(in)                :: keyword
    integer(shortInt),intent(in)           :: opt
    integer(shortInt)                      :: idx
    character(100), parameter :: Here = 'getFromCL_char (commandLineUI.f90)'

    if(.not.parsed) call parseCL()

    ! Find index of the option
    idx = findOptionIdx(keyword)

    ! Check if option is defined
    if( idx == targetNotFound) then
      call fatalError(Here,'Option with a keyword: '//trim(adjustl(keyword))//' was not defined.')
    end if

    ! Check if option was invoked (has values)
    if (.not.options(idx) % isLoaded) then
      call fatalError(Here,'Option with a keyword: '//trim(adjustl(keyword))//' was not invoked.')
    end if

    ! Check that value of opt checks out
    if ( opt < 1) call fatalError(Here,'Option argument index must be < 1')
    if ( opt > size(options(idx) % optArgs) ) then
      call fatalError(Here,'Option '//trim(adjustl(keyword))// 'has only'// &
                            numToChar(size(options(idx) % optArgs))//' input arguments')
    end if

    ! Check that type of option argument matches
    if (options(idx) % optArgs(opt) % type /= CHAR_TYPE) then
      call fatalError(Here,numToChar(opt)//'-th argument of option '//trim(adjustl(keyword))//&
                           'is a character')
    end if

    ! Read the value
    var = options(idx) % optArgs(opt) % val

  end subroutine getFromCL_char

  !!
  !! Display parser error and stop execution
  !!
  subroutine parseError(Error)
    character(*), intent(in) :: Error
    print *, "Command line parse error"
    print *, Error
    print *
    call displayClOptions()
    stop
  end subroutine parseError

  !!
  !! Print description of all command line options
  !!
  subroutine displayClOptions()
    integer(shortInt)        :: i

    print *, "Usage: (Executable) <input File>  [Options]"
    print *
    print *, "Following options are available: "

    do i=1,size(options)
      call options(i) % display()
    end do

  end subroutine displayClOptions

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Option descriptor procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Initialise option Descriptor
  !!
  subroutine init_optionDescriptor(self,keyword,Narg,argTypes,help)
    class(optionDescriptor), intent(inout)     :: self
    character(*), intent(in)                   :: keyword
    integer(shortInt), intent(in)              :: Narg
    character(*),dimension(:),intent(in)       :: argTypes
    character(*),intent(in)                    :: help
    integer(shortInt)                          :: i
    character(100),parameter :: Here ='init_optionDescriptor (commandLinieUI.f90)'

    ! Load arguments
    self % keyword = trim(adjustl(keyword))
    self % help = help

    if(Narg < 0) call parseError('Number of arguments must be non -ve')

    allocate( self % optArgs(Narg))

    ! Check if keyword is zero length
    if (len_trim(keyword) == 0) call parseError('Keyword cannot consist of blanks')

    ! Error if argTypes is wong size. Ignore argTypes if Narg == 0
    if( size(argTypes) /= Narg .and. Narg /= 0 ) then
       call parseError('Argument type needs to be given for all arguments')
    end if

    do i=1,Narg
      select case(trim(argTypes(i)))
        case('int')
          self % optArgs(i) % type = INT_TYPE

        case('real')
          self % optArgs(i) % type = REAL_TYPE

        case('char')
          self % optArgs(i) % type = CHAR_TYPE

        case default
          call parseError("Unrecognised argument type: "//trim(argTypes(i))//" Must be &
                              & 'int', 'real' or 'char' ") 
      end select
    end do

    ! Hash keyword
    call FNV_1(self % keyword, self % keywordHash)

    ! Set to unloaded
    self % isLoaded = .false.

  end subroutine init_optionDescriptor

  !!
  !! Read arguments of an option from command line,
  !! given position(pos) of the keyword and number of CL arguments(argCount)
  !!   Increment pos by the number of arguments at the end
  !!
  subroutine read_optionDescriptor(self,pos,argCount)
    class(optionDescriptor), intent(inout) :: self
    integer(shortInt), intent(inout)       :: pos
    integer(shortInt), intent(in)          :: argCount
    character(:),allocatable               :: string
    integer(shortInt)                      :: i, argLen
    integer(shortInt)                      :: Narg
    character(100), parameter :: Here = 'read_optionDescriptor (commandLinieUI.f90)'

    ! Change status to loaded
    self % isLoaded = .true.

    ! Get number of option arguments
    Narg = size(self % optArgs)

    ! Check for overflow of command arguments
    if (pos + Narg > argCount) then
      call parseError("There is not enough CL arguments to read option: "//self % keyword)
    end if

    do i=1,Narg
      ! Get length of command argument
      call get_command_argument(pos + i, length = argLen)

      ! Allocate string for a command line argument
      if(allocated(string)) deallocate(string)
      allocate(character(argLen) :: string)

      ! Get command argument
      call get_command_argument(pos + i, string)

      string = trim(adjustl(string))

      ! Load string
      self % optArgs(i) % val = string

      ! Verify string
      call self % optArgs(i) % verify(self % keyword)
    end do

    ! Increment pos
    pos = pos + Narg

  end subroutine read_optionDescriptor

  !!
  !! Returns .true. if options have the same keyword
  !!
  elemental function eq_optionDescriptor(LHS,RHS) result(equal)
    class(optionDescriptor), intent(in) :: LHS
    class(optionDescriptor), intent(in) :: RHS
    logical(defBool)                    :: equal

    ! Compare hashes
    equal = LHS % keywordHash == RHS % keywordHash

    ! If hases are the same compare keywords
    if (equal) equal = LHS % keyword == RHS % keyword

  end function eq_optionDescriptor

  !!
  !! Print definition of the option to the screen
  !!
  subroutine display_optionDescriptor(self)
    class(optionDescriptor), intent(in) :: self
    character(:),allocatable            :: keyword
    character(nameLen)                  :: format2, format
    integer(shortInt)                   :: start,end, lines, i

    ! Append keyword with argument types
    keyword = self % keyword
    do i=1,size(self % optArgs)
      keyword = keyword // ' '// self % optArgs(i) % typeChar()
    end do

    ! Determine how meany lines need to be printed
    lines = ceiling( len(self % help)/50.0)

    ! Print first line
    start = 1
    end = min(start+50-1,len(self % help))

    format  = '(TR3,A'//numToChar(len(keyword))//',T25,A'//numToChar(end-start+1)//')'
    print format, keyword, self % help(start:end)

    ! Print subsequent lines
    do i=2,lines
      start =  1 + (i-1) * 50
      end = min(start+50-1,len(self % help))

      format2 = '(T25, A'//numToChar(end-start+1)//')'
      print format2, self % help(start:end)
    end do

  end subroutine display_optionDescriptor

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Option Argument procedures
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
  !!
  !! Return 6 character string describing required argument type
  !!
  elemental function typeChar_optionArg(self) result(typeChar)
    class(optionArg), intent(in) :: self
    character(6)                 :: typeChar

    select case(self % type)
      case(INT_TYPE)
        typeChar = '<int> '

      case(REAL_TYPE)
        typeChar = '<real>'

      case(CHAR_TYPE)
        typeChar = '<char>'

      case default
        typeChar = '<udef>'

    end select
  end function typeChar_optionArg

  !!
  !! Returns error if loaded string does not match type
  !!
  subroutine verify_optionArg(self, option)
    class(optionArg), intent(in) :: self
    character(*),intent(in)      :: option
    integer(shortInt)            :: readErr
    integer(shortInt)            :: dummyInt
    real(defReal)                :: dummyReal
    character(100), parameter :: Here = 'verify_optionArg (commandLinieUI.f90)'

    readErr = 0
    select case(self % type)
      case(INT_TYPE)
        read(self % val, '(I30)' , iostat = readErr) dummyInt
        if ( readErr /= 0) call parseError(option//': Wrong input given to integer value')

      case(REAL_TYPE)
        ! Check if the given value is an integer. dummyInt must be present
        read(self % val, '(I30)' , iostat = readErr) dummyInt
        if( readErr == 0) call parseError(option//': Integer given in place of real value')

        ! Check if the givenvalue is a real
        read(self % val, '(G30.30)' , iostat = readErr) dummyReal
        if ( readErr /= 0) call parseError(option//': Wrong input given to real value')
    end select
  end subroutine verify_optionArg

end module commandLineUI
