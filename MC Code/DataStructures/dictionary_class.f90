module dictionary_class

  use numPrecision
  use genericProcedures, only: linFind, searchError, fatalError, targetNotFound

  implicit none
  private
  ! *** Below is not up to date !
  ! Due to bug in verision 4.8 of gfortran it is necessary to compile this code with newer compiler.
  ! This code was tested after compilation with gfortran 5.2
  !
  ! It is necessary to predefine size of stored char. Gfortran 4.8 does not support allocatable
  ! character with deffered length.
  !
  ! WARNING: ***************************************************************************************
  ! Becouse the current default assigment(=) for dictionary is deep copy using
  ! dictionary valued functions will result in memory leaks. Following statment needs to be avoided
  !
  !  <dictionary variable> = functionDict()
  !
  !  function functionDict()
  !     type(dictionary), intent(out)  :: functionDict
  !                ....
  ! END OF WARNING**********************************************************************************
  !
  ! For now the dictionary is limited to following enteries:
  ! ->   scalar defReal                                       - real(defReal)      :: a
  ! ->   1D array of defReal                                  - real(defReal)      :: a(:)
  ! ->   scalar shortInt                                      - integer(shortInt)  :: i
  ! ->   1D array of shortInt                                 - integer(shortInt)  :: i(:)
  ! ->   scalar character of length "charLen" defined below   - character(charLen) :: c
  ! ->   array of charcters of length "charLen" defined below - character(charLen) :: c(:)
  ! ->   another dictionary                                   - type(dictionary)   :: dict
  !
  ! To add additional structure <type> to store it is necessary to :
  ! 0) Define new type parameter for <type>
  ! 1) Add put_<type>   subroutine to dictDontent class
  ! 2) Add another entery for <type> to select type in deepCopy_dictCont subroutine
  ! 3) Add store_<type> subroutine to dictionary class
  ! 4) Add get<type> function to dictionary class
  ! 5) Add keys<type> function
  !
  ! To increase maximum rank would be messy. Don't do it! If you are determined you can probably
  ! hack existing code by storing higher rank arrays in rank1 array. Their structure should then
  ! be stored seperatly in integers in dictContent class.
  !
  ! When using getReal or getRealArray it is possible to provide keyword associated with an integer
  ! to atutomaticly conver integer to real e.g:
  ! For following dictionary dict:
  ! integerKey 3;
  ! realKey    3.1;
  !
  ! Function: dict % getReal('integerKey') will return 3.0;
  !
  ! Output of getDict is a shallow copy of the dictionary. Therefore it is possible to provide it
  ! as an argument to a subroutine without worring about a memory leak
  !

  integer(shortInt),parameter,public :: charLen  = max(nameLen,pathLen)
  integer(shortInt),parameter        :: empty    = 0
  integer(shortInt),parameter        :: numInt   = 1
  integer(shortInt),parameter        :: numReal  = 2
  integer(shortInt),parameter        :: word     = 3
  integer(shortInt),parameter        :: nestDict = 4
  integer(shortInt),parameter        :: arrInt   = 5
  integer(shortInt),parameter        :: arrReal  = 6
  integer(shortInt),parameter        :: arrWord  = 7

  integer(shortInt),parameter        :: defStride = 20

  !!
  !! Type to store a single entery in a dictionary
  !! Uses a single allocatable variable for diffrent type of contant
  !!
  type,public :: dictContent
    ! Allocatable space for all content types
    integer(shortInt)                           :: int0_alloc
    integer(shortInt),dimension(:),allocatable  :: int1_alloc
    real(defReal)                               :: real0_alloc
    real(defReal),dimension(:),allocatable      :: real1_alloc
    character(charLen)                          :: char0_alloc
    character(charLen),dimension(:),allocatable :: char1_alloc
    ! *** Note that dictionary is defined as ponter not allocatable
    ! *** This is becouse gfortran > 7.0 supports circular derived types with
    ! *** allocatable keyword. This line may change in a future
    class(dictionary), pointer                  :: dict0_alloc => null()

    ! dictContent type ID
    integer(shortInt)                           :: type = empty
  contains
    procedure :: kill          => kill_dictCont
    procedure :: copy          => copy_dictCont
    procedure :: getType       => getType_dictContent

  end type dictContent

  !!
  !! Dictionary type
  !!
  type, public :: dictionary
    private
    ! Dictionary storage array
    character(nameLen),dimension(:),allocatable  :: keywords
    type(dictContent),dimension(:), allocatable  :: entries
    ! Dictionary state information
    integer(shortInt)                            :: maxSize = 0         ! Maximum size of a dictionary
    integer(shortInt)                            :: dictLen = 0         ! Current size of the dictionary
    integer(shortInt)                            :: stride  = defStride ! Extension size when reached maxSize
  contains

    generic    :: store => store_real      ,&
                           store_realArray ,&
                           store_int       ,&
                           store_intArray  ,&
                           store_char      ,&
                           store_charArray ,&
                           store_dict
    procedure  :: isPresent
    procedure  :: getReal
    procedure  :: getRealArray
    procedure  :: getInt
    procedure  :: getIntArray
    procedure  :: getChar
    procedure  :: getCharArray
    procedure  :: getDict

    procedure  :: keysReal
    procedure  :: keysRealArray
    procedure  :: keysInt
    procedure  :: keysIntArray
    procedure  :: keysChar
    procedure  :: keysCharArray
    procedure  :: keysDict
    !procedure  :: keysDict_type
    procedure  :: keys

    generic    :: assignment(=) => copy

    procedure  :: copy    => copy_dictionary
    procedure  :: init
    procedure  :: kill => kill_dictionary
    final      :: final_dictionary

    procedure, private :: extendBy
    procedure, private :: getEmptyIdx

    procedure, private :: store_real
    procedure, private :: store_realArray
    procedure, private :: store_int
    procedure, private :: store_intArray
    procedure, private :: store_char
    procedure, private :: store_charArray
    procedure, private :: store_dict

  end type dictionary

contains

  !!
  !! Helper function to get next empty index in a dictionary storage array
  !! If storage arrays are full it extends dictionary by stride
  !!
  function getEmptyIdx(self,keyword) result (idx)
    class(dictionary), intent(inout)   :: self
    character(nameLen), intent(in)     :: keyword
    integer(shortInt)                  :: idx
    integer(shortInt)                  :: i
    character(100),parameter           :: Here='getEmptyIdx (dictionary_class.f90)'

    ! Extend dictionary if required
    if (self % dictLen >= self % maxSize) then
      call self % extendBy(self % stride)

    end if

    ! Check whether keywods contains characters

    if (adjustl(keyword) == adjustl('')) then
      call fatalError(Here,"Keyword contains only blanks : '' ")
    end if

    ! Find if keyword is alrady present
    i = linFind(self % keywords,keyword)
    if (i > 0 ) call fatalError(Here, 'Keyword: ' // keyword // ' is already present in dictionary')

    ! Increase counter and return avalible index
    self % dictLen = self % dictLen + 1
    idx = self % dictLen

  end function getEmptyIdx

  !!
  !! Extends size of dictionary storage arrays by variable stride
  !!
  subroutine extendBy(self,stride)
    class(dictionary), intent(inout)             :: self
    integer(shortInt)                            :: stride
    character(nameLen),dimension(:),allocatable  :: keywords
    type(dictContent),dimension(:), allocatable  :: entries
    integer(shortInt)                            :: newSize, oldSize, i
    character(100),parameter                     :: Here='extendBy (dictionary_class.f90)'

    if(.not.(allocated(self % keywords).and.allocated(self % entries))) then
      call fatalError(Here,'An attempt was made to extend uninitialised dictionary')

    end if

    oldSize = self % maxSize
    newSize = oldSize + stride

    allocate(keywords(newSize))
    allocate(entries(newSize) )

    keywords(1:oldSize) = self % keywords

    do i=1,oldSize
     entries(i) = self % entries(i)

    end do

    deallocate( self % keywords)
    deallocate( self % entries)

    call move_alloc(keywords, self % keywords)
    call move_alloc(entries , self % entries )

    self % maxSize = newSize

  end subroutine extendBy


  !!
  !! Initialise a dictionary
  !! Choose size of storage array and allocate them
  !!
  subroutine init(self,maxSize,stride)
    class(dictionary), intent(inout)         :: self
    integer(shortInt), intent(in)            :: maxSize
    integer(shortInt), intent(in), optional  :: stride
    character(100),parameter                 :: Here='init (dictionary_class.f90)'

    if(allocated(self % keywords).or.allocated(self % entries)) then
      call fatalError(Here,'Attempting to reinitialise a dictionary is forbidden')
    end if

    self % maxSize = maxSize

    allocate(self % keywords(maxSize))
    allocate(self % entries(maxSize))


    ! Keywords can (perhaps?) allocate with some garbage inside. Make shure all enteries are blank.
    self % keywords = ''

    if (present(stride)) self % stride = stride

  end subroutine init

  !!
  !! Dealloacte storage space
  !! Returns dictionary to uninitialised state (including stride)
  !!
  subroutine kill_dictionary(self)
    class(dictionary), intent(inout) :: self
    integer(shortInt)                :: i
    logical(defBool)                 :: keysAllocated
    logical(defBool)                 :: entAllocated
    character(100),parameter         :: Here='kill (dictionary_class.f90)'

    keysAllocated = allocated(self % keywords)
    entAllocated  = allocated(self % entries)

    if(keysAllocated .neqv. entAllocated) then
      call fatalError(Here,'Immposible state. Keywors or entries is allocated without the other')

    elseif (keysAllocated .and. entAllocated) then

      ! Expression below could extend only to self % dictLen but lets make double shure we
      ! kill all data. There is no need to optimise this and it is more robust this way.
      do i=1,self % maxSize
        call self % entries(i) % kill()
      end do
      deallocate(self % keywords)
      deallocate(self % entries)

      self % maxSize = 0
      self % dictLen = 0
      self % stride = defStride

    end if

  end subroutine kill_dictionary

  !!
  !! Finalisation Subroutine
  !!
  subroutine final_dictionary(self)
    type(dictionary),intent(inout) :: self

    call self % kill()
  end subroutine


  !!
  !! Copy one dictionary to another
  !! Overwrites LHS
  !!
  subroutine copy_dictionary(LHS,RHS)
    class(dictionary), intent(inout) :: LHS
    class(dictionary), intent(in)    :: RHS
    integer(shortInt)                :: rhsSize, stride
    integer(shortInt)                :: i

    ! Clean LHS dictionary
    call LHS % kill()

    ! Reainitialise LHS dictionary
    rhsSize = size( RHS % keywords)
    stride = RHS % stride

    call LHS % init(rhsSize,stride)

    ! Copy Keywords and entries
    LHS % keywords = RHS % keywords

    do i=1,RHS % maxSize
      call LHS % entries(i) % copy(RHS % entries(i) )

    end do
    ! Copy state settings
    LHS % dictLen = RHS % dictLen
    LHS % maxSize = RHS % maxSize
    LHS % stride  = RHS % stride

  end subroutine copy_dictionary

  !!
  !! Checks if given keyword is present in a dictionary
  !!
  function isPresent(self,keyword) result(isIt)
    class(dictionary), intent(in) :: self
    character(*), intent(in)      :: keyword
    logical(defBool)              :: isIt
    integer(shortInt)             :: idx

    idx  = linFind(self % keywords, keyword)
    isIt = .not.(idx == targetNotFound)

  end function isPresent

  !!
  !! Reads a real rank 0 entery from a dictionary
  !! If keyword is associated with an integer it converts it to real
  !!
  function getReal(self,keyword) result(value)
    class(dictionary), intent(in)  :: self
    character(*),intent(in)        :: keyword
    real(defReal)                  :: value
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getReal (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    select case (self % entries(idx) % getType())
      case(numReal)
        value = self % entries(idx) % real0_alloc

      case(numInt)
        value = real(self % entries(idx) % int0_alloc, defReal)

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a real or int')

    end select

  end function getReal

  !!
  !! Reads a real rank 1 from a dictionary
  !! If keyword is associated with an integer it converts it to real
  !!
  function getRealArray(self,keyword) result(value)
    class(dictionary), intent(in)          :: self
    character(*),intent(in)                :: keyword
    real(defReal),dimension(:),allocatable :: value
    integer(shortInt)                      :: idx
    character(100),parameter               :: Here='getRealArray (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    select case (self % entries(idx) % getType())
      case(arrReal)
        value = self % entries(idx) % real1_alloc

      case(arrInt)
        value = real(self % entries(idx) % int1_alloc, defReal)

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a real array or int array')

    end select

  end function getRealArray

  !!
  !! Reads a integer rank 0 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !!
  function getInt(self,keyword) result(value)
    class(dictionary), intent(in)  :: self
    character(*),intent(in)        :: keyword
    integer(shortInt)              :: value
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getInt (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    select case (self % entries(idx) % getType())
      case(numInt)
        value = self % entries(idx) % int0_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not an integer')

    end select

  end function getInt

  !!
  !! Reads a integer rank 1 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !!
  function getIntArray(self,keyword) result(value)
    class(dictionary), intent(in)              :: self
    character(*),intent(in)                    :: keyword
    integer(shortInt),dimension(:),allocatable :: value
    integer(shortInt)                          :: idx
    character(100),parameter                   :: Here='getIntArray (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    select case (self % entries(idx) % getType())
      case(arrInt)
        value = self % entries(idx) % int1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not integer array')

    end select

  end function getIntArray

  !!
  !! Reads a character rank 0 from a dictionary
  !!
  function getChar(self,keyword) result(value)
    class(dictionary), intent(in)  :: self
    character(*),intent(in)        :: keyword
    character(charLen)             :: value
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getChar (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    select case (self % entries(idx) % getType())
      case(word)
        value = self % entries(idx) % char0_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a character')

    end select


  end function getChar

  !!
  !! Reads a character rank 1 from a dictionary
  !!
  function getCharArray(self,keyword) result(value)
    class(dictionary), intent(in)               :: self
    character(*),intent(in)                     :: keyword
    character(charLen),dimension(:),allocatable :: value
    integer(shortInt)                           :: idx
    character(100),parameter                    :: Here='getCharArray (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    select case (self % entries(idx) % getType())
      case(arrWord)
        value = self % entries(idx) % char1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a character array')

    end select

  end function getCharArray

  !!
  !! Reads a dictionary rank 0 from a dictionary
  !!
  function getDict(self,keyword) result(value)
    class(dictionary), intent(in)  :: self
    character(*),intent(in)        :: keyword
    type(dictionary)               :: value
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getChar (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    select case (self % entries(idx) % getType())
      case(nestDict)
        value = self % entries(idx) % dict0_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a dictionary')

    end select

  end function getDict

  !!
  !! Returns an array of all keywords associated with a real rank 0
  !!
  function keysReal(self) result(keys)
    class(dictionary), intent(in)                :: self
    character(nameLen),dimension(:), allocatable :: keys
    logical(defBool),dimension(:),allocatable    :: mask
    integer(shortInt)                            :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == numReal)

    keys = pack(self % keywords(1:L), mask)

  end function keysReal

  !!
  !! Returns an array of all keywords associated with a real rank 1
  !!
  function keysRealArray(self) result(keys)
    class(dictionary), intent(in)                :: self
    character(nameLen),dimension(:), allocatable :: keys
    logical(defBool),dimension(:),allocatable    :: mask
    integer(shortInt)                            :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == arrReal)

    keys = pack(self % keywords(1:L), mask)

  end function keysRealArray

  !!
  !! Returns an array of all keywords associated with an integer rank 0
  !!
  function keysInt(self) result(keys)
    class(dictionary), intent(in)                :: self
    character(nameLen),dimension(:), allocatable :: keys
    logical(defBool),dimension(:),allocatable    :: mask
    integer(shortInt)                            :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == numInt)

    keys = pack(self % keywords(1:L), mask)

  end function keysInt

  !!
  !! Returns an array of all keywords associated with an integer rank 0
  !!
  function keysIntArray(self) result(keys)
    class(dictionary), intent(in)                :: self
    character(nameLen),dimension(:), allocatable :: keys
    logical(defBool),dimension(:),allocatable    :: mask
    integer(shortInt)                            :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == arrInt)

    keys = pack(self % keywords(1:L), mask)

  end function keysIntArray

  !!
  !! Returns an array of all keywords associated with an character rank 0
  !!
  function keysChar(self) result(keys)
    class(dictionary), intent(in)                :: self
    character(nameLen),dimension(:), allocatable :: keys
    logical(defBool),dimension(:),allocatable    :: mask
    integer(shortInt)                            :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == word)

    keys = pack(self % keywords(1:L), mask)

  end function keysChar

  !!
  !! Returns an array of all keywords associated with an character rank 1
  !!
  function keysCharArray(self) result(keys)
    class(dictionary), intent(in)                :: self
    character(nameLen),dimension(:), allocatable :: keys
    logical(defBool),dimension(:),allocatable    :: mask
    integer(shortInt)                            :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == arrWord)

    keys = pack(self % keywords(1:L), mask)

  end function keysCharArray

  !!
  !! Returns an array of all keywords associated with a dictionary rank 0
  !!
  function keysDict(self) result(keys)
    class(dictionary), intent(in)                :: self
    character(nameLen),dimension(:), allocatable :: keys
    logical(defBool),dimension(:),allocatable    :: mask
    integer(shortInt)                            :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == nestDict)

    keys = pack(self % keywords(1:L), mask)

  end function keysDict

!  !!
!  !! Return all dictionarys with a specific type keword
!  !!
!  function keysDict_type(self,dictType) result(keys)
!    class(dictionary), intent(in)               :: self
!    character(*), intent(in)                    :: dictType
!    character(nameLen),dimension(:),allocatable :: keys
!    character(nameLen),dimension(:),allocatable :: allDict
!    logical(defBool), dimension(:), allocatable :: mask
!    integer(shortInt)                           :: L,i
!    character(charLen)                          :: typeTemp
!    type(dictionary)                            :: tempDict
!
!
!    ! Due to compiler bugs copy keysDict_all
!    L = self % dictLen
!    allocate( mask(L) )
!
!    mask = (self % entries(1:L) % getType() == nestDict)
!
!    allDict = pack(self % keywords(1:L), mask)
!    deallocate(mask)
!
!    ! Find all indexes in keys that match requested type
!    L = size(allDict)
!    allocate(mask(L))
!
!    do i=1,L
!      tempDict = self % getDict(allDict(i))
!
!      if (tempDict % isPresent('type')) then
!        typeTemp = tempDict % getChar('type')
!        mask(i) = (trim(adjustl(dictType)) == trim(adjustl(typeTemp)))
!
!      else
!        mask(i) = .false.
!
!      end if
!    end do
!
!    keys = pack(allDict,mask)
!
!  end function keysDict_type


  !!
  !! Returns an array of all keywords
  !!
  function keys(self)
    class(dictionary), intent(in)                :: self
    character(nameLen),dimension(:), allocatable :: keys
    integer(shortInt)                            :: L

    L    = self % dictLen
    keys = self % keywords(1:L)

  end function keys

  !!
  !! Stores a real rank 0 in dictionary
  !!
  subroutine store_real(self,keywordArgument,entry)
    class(dictionary), intent(inout)  :: self
    character(*), intent(in)          :: keywordArgument
    real(defReal), intent(in)         :: entry
    character(nameLen)                :: keyword
    integer(shortInt)                 :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    ! Load into dictionary content
    self % entries(idx) % real0_alloc = entry
    self % entries(idx) % type = numReal

  end subroutine store_real

  !!
  !! Stores a real rank 1 in dictionary
  !!
  subroutine store_realArray(self,keywordArgument,entry)
    class(dictionary), intent(inout)        :: self
    character(*), intent(in)                :: keywordArgument
    real(defReal), dimension(:), intent(in) :: entry
    character(nameLen)                      :: keyword
    integer(shortInt)                       :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    ! Load into dictionary content
    self % entries(idx) % real1_alloc = entry
    self % entries(idx) % type = arrReal

  end subroutine store_realArray

  !!
  !! Stores a integer rank 0 in dictionary
  !!
  subroutine store_int(self,keywordArgument,entry)
    class(dictionary), intent(inout)  :: self
    character(*), intent(in)          :: keywordArgument
    integer(shortInt), intent(in)     :: entry
    character(nameLen)                :: keyword
    integer(shortInt)                 :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    ! Load into dictionary content
    self % entries(idx) % int0_alloc = entry
    self % entries(idx) % type = numInt

  end subroutine store_int

  !!
  !! Stores a integer rank 1 in dictionary
  !!
  subroutine store_intArray(self,keywordArgument,entry)
    class(dictionary), intent(inout)            :: self
    character(*), intent(in)                    :: keywordArgument
    integer(shortInt), dimension(:), intent(in) :: entry
    character(nameLen)                          :: keyword
    integer(shortInt)                           :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    ! Load into dictionary content
    self % entries(idx) % int1_alloc = entry
    self % entries(idx) % type = arrInt

  end subroutine store_intArray

  !!
  !! Stores a character rank 0 in dictionary
  !!
  subroutine store_char(self,keywordArgument,entry)
    class(dictionary), intent(inout)  :: self
    character(*), intent(in)          :: keywordArgument
    character(*), intent(in)          :: entry
    character(nameLen)                :: keyword
    integer(shortInt)                 :: idx

    keyword   = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    ! Load into dictionary content
    self % entries(idx) % char0_alloc = entry
    self % entries(idx) % type = word

  end subroutine store_char

  !!
  !! Stores a character rank 1 in dictionary
  !!
  subroutine store_charArray(self,keywordArgument,entry)
    class(dictionary), intent(inout)            :: self
    character(*), intent(in)                    :: keywordArgument
    character(*), dimension(:), intent(in)      :: entry
    character(nameLen)                          :: keyword
    integer(shortInt)                           :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    ! Load into dictionary content
    self % entries(idx) % char1_alloc = entry
    self % entries(idx) % type = arrWord

  end subroutine store_charArray


  !!
  !! Stores a dictionary rank 0 in dictionary
  !!
  subroutine store_dict(self,keywordArgument,entry)
    class(dictionary), intent(inout)            :: self
    character(*), intent(in)                    :: keywordArgument
    type(dictionary), intent(in)                :: entry
    character(nameLen)                          :: keyword
    integer(shortInt)                           :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword

    ! Load into dictionary content
    allocate(self % entries(idx) % dict0_alloc)

    self % entries(idx) % dict0_alloc = entry
    self % entries(idx) % type = nestDict

  end subroutine store_dict

  !!
  !! Copies dictContent
  !!
  subroutine copy_dictCont(LHS,RHS)
    class(dictContent), intent(inout) :: LHS
    type(dictContent), intent(in)     :: RHS

    ! Remove contents of a target
    call LHS % kill()

    ! Copy dictContent type
    LHS % type = RHS % type

    ! Copy individual enteries
    LHS % int0_alloc = RHS % int0_alloc
    if (allocated( RHS % int1_alloc)) LHS % int1_alloc = RHS % int1_alloc


    LHS % real0_alloc = RHS % real0_alloc
    if (allocated( RHS % real1_alloc)) LHS % real1_alloc = RHS % real1_alloc


    LHS % char0_alloc = RHS % char0_alloc
    if (allocated( RHS % char1_alloc)) LHS % char1_alloc = RHS % char1_alloc

    if (associated( RHS % dict0_alloc)) then
      allocate(LHS % dict0_alloc)
      LHS % dict0_alloc = RHS % dict0_alloc
    end if

  end subroutine copy_dictCont


  !!
  !! Deallocates dictContent
  !!
  subroutine kill_dictCont(self)
    class(dictContent), intent(inout) :: self

    ! Deallocate allocatable components
    if(allocated(self % int1_alloc)) deallocate (self % int1_alloc)

    if(allocated(self % real1_alloc)) deallocate (self % real1_alloc)

    if(allocated(self % char1_alloc)) deallocate (self % char1_alloc)

    ! Clean nested dictionaries. Kill before deallocation to avoid memory leaks
    if(associated(self % dict0_alloc)) then
      call self % dict0_alloc % kill()
      deallocate (self % dict0_alloc)
    end if

    ! Change type of content to empty
    self % type = empty

  end subroutine kill_dictCont

  !!
  !! Returns type of a given dictContent
  !!
  elemental function getType_dictContent(self) result(type)
    class(dictContent), intent(in)     :: self
    integer(shortInt)                  :: type

    type = self % type

  end function getType_dictContent

end module dictionary_class
