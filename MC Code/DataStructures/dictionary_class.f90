module dictionary_class

  use numPrecision
  use genericProcedures, only: linFind, searchError, fatalError, targetNotFound

  implicit none
  private
  ! Implementation of the dictionary uses finalisation procedures that are implemented into gfortran
  ! version 4.9 and higher. There is a bug in gfortran untill version 7.0 which causes false warning
  ! from -Wsurprising to appear with final procedure. This warning mey be suppressed with
  ! -Wno-surprising flag.
  !
  ! Becouse of the finalisation local dictionary defined inside procedure will be properly
  ! deallocated. Thus:
  !
  !   subroutine memLeak(dictIn)
  !     class(dictionary),intent(in) :: dictIN
  !     type(dictionary)             :: locDict
  !     locDict = dictIN
  !   end subroutine
  !
  ! will not cause a memory leak despite lack of "call locDict % kill()" before "end subroutine"
  !
  ! Size of the stored char is predefined to be maximum of nameLen and pathLen.
  ! Maximum size of keyword is equal to nameLen
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
  ! 0) Define new type parameter for <type> type
  ! 1) Add additional allocatable space for the <type> in dictContent
  ! 3) Add another entery for <type> to copy_dictContent
  ! 4) Add another entery for <type> to kill_dictContent
  ! 5) Add store_<type> subroutine to dictionary class. Add it to generic store
  ! 6) Add get<type> subroutine to dictionary class. Add it to generic store
  ! 7) Add get<type>orDefault subroutine to dictionary class. Ad it to generic storeOrDefault
  ! 8) Add keys<type> function
  !
  ! There are two generic subroutines to extract data from a dictionary.
  ! 1) get(<type>,keyword) subroutine
  !    - error is returned if keyword is not present in dictionary
  !    - error is returned if <type> does not match type of content under keyword
  !    - exception to above is for <type>=real which can read integer content
  !    - hard copy of the content is made including subdictionaries
  !    - if <type> is array it needs to be allocatable or pointer. It will be deallocated and
  !      reallocated inside the subroutine
  !    - if <type> is scalar character and its length is too short to fit the trimed content under
  !      keyword error is returned
  !    - if <type> is array character and its length is too short to fit the trimmed content inder
  !      keyword error is returned
  !
  ! 2) getOrDefault(<type>,keyword,default)
  !    - getOrDefault does not support <type> = dictionary
  !    - default needs to be same type and rank as <type>
  !    - if keyword is not find in dictionary default is returned
  !    - if <type> is scalar or array character error is returned if len(<type>) < len(default)
  !      or if <type> is too short to fit trimmed content under keyword
  !
  !
  ! When retreaving real or real Array it is possible to provide keyword associated with an integer
  ! For following dictionary dict:
  ! integerKey 3;
  ! realKey    3.1;
  !
  ! Subroutine: dict % getReal(r,'integerKey') will put 3.0 into r
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
    ! *** Note that dictionary is defined as pointer not allocatable
    ! *** This is becouse gfortran < 7.0 does not supports circular derived types with
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
    !*** Obsolete functional access interface will be removed
    procedure  :: getReal      => getReal_old
    procedure  :: getRealArray => getRealArray_old
    procedure  :: getInt       => getInt_old
    procedure  :: getIntArray  => getIntArray_old
    procedure  :: getChar      => getChar_old
    procedure  :: getCharArray => getCharArray_old
    procedure  :: getDict      => getDict_old
    !***********************

    generic    :: get => getReal_new,&
                         getRealArray_alloc_new,&
                         getRealArray_ptr_new,&
                         getInt_new,&
                         getIntArray_alloc_new,&
                         getIntArray_ptr_new,&
                         getChar_new,&
                         getCharArray_alloc_new,&
                         getCharArray_ptr_new,&
                         getDict_new

    procedure,private :: getReal_new
    procedure,private :: getRealArray_alloc_new
    procedure,private :: getRealArray_ptr_new
    procedure,private :: getInt_new
    procedure,private :: getIntArray_alloc_new
    procedure,private :: getIntArray_ptr_new
    procedure,private :: getChar_new
    procedure,private :: getCharArray_alloc_new
    procedure,private :: getCharArray_ptr_new
    procedure,private :: getDict_new

    generic :: getOrDefault => getOrDefault_real ,&
                               getOrDefault_realArray_alloc ,&
                               getOrDefault_realArray_ptr ,&
                               getOrDefault_int ,&
                               getOrDefault_intArray_alloc ,&
                               getOrDefault_intArray_ptr ,&
                               getOrDefault_char ,&
                               getOrDefault_charArray_alloc ,&
                               getOrDefault_charArray_ptr

    procedure,private :: getOrDefault_real
    procedure,private :: getOrDefault_realArray_alloc
    procedure,private :: getOrDefault_realArray_ptr
    procedure,private :: getOrDefault_int
    procedure,private :: getOrDefault_intArray_alloc
    procedure,private :: getOrDefault_intArray_ptr
    procedure,private :: getOrDefault_char
    procedure,private :: getOrDefault_charArray_alloc
    procedure,private :: getOrDefault_charArray_ptr

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
    i = linFind(self % keywords(1:self % dictLen), keyword)
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
  !! Loads a real rank 0 from a dictionary into provided variable
  !! If keyword is associated with an integer it converts it to real
  !!
  subroutine getReal_new(self,value,keyword)
    class(dictionary), intent(in)  :: self
    real(defReal),intent(inout)    :: value
    character(*),intent(in)        :: keyword
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

  end subroutine getReal_new

  !!
  !! Loads a real rank 1 from dictionary into provided variable
  !! If keyword is associated with an integer it converts it to real
  !! Variable needs to be allocatable. It will be deallocated before assignment
  !!
  subroutine getRealArray_alloc_new(self,value,keyword)
    class(dictionary), intent(in)                        :: self
    real(defReal),dimension(:),allocatable,intent(inout) :: value
    character(*),intent(in)                              :: keyword
    integer(shortInt)                                    :: idx
    character(100),parameter                             :: Here='getRealArray_alloc (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    if(allocated(value)) deallocate(value)

    select case (self % entries(idx) % getType())
      case(arrReal)
        value = self % entries(idx) % real1_alloc

      case(arrInt)
        value = real(self % entries(idx) % int1_alloc, defReal)

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a real array or int array')

    end select

  end subroutine getRealArray_alloc_new

  !!
  !! Loads a real rank 1 from dictionary into provided variable
  !! If keyword is associated with an integer it converts it to real
  !! Variable needs to be pointer. It will be deallocated before assignment
  !!
  subroutine getRealArray_ptr_new(self,value,keyword)
    class(dictionary), intent(in)                      :: self
    real(defReal),dimension(:),pointer,intent(inout)   :: value
    character(*),intent(in)                            :: keyword
    integer(shortInt)                                  :: idx, N
    character(100),parameter                           :: Here='getRealArray_ptr (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    if(associated(value)) deallocate(value)

    select case (self % entries(idx) % getType())
      case(arrReal)
        N = size (self % entries(idx) % real1_alloc)
        allocate(value(N))
        value = self % entries(idx) % real1_alloc

      case(arrInt)
        value = real(self % entries(idx) % int1_alloc, defReal)

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a real array or int array')

    end select

  end subroutine getRealArray_ptr_new

  !!
  !! Loads a integer rank 0 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !!
  subroutine getInt_new(self,value,keyword)
    class(dictionary), intent(in)  :: self
    integer(shortInt),intent(inout):: value
    character(*),intent(in)        :: keyword
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

  end subroutine getInt_new

  !!
  !! Loads a integer rank 1 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !! Variable needs to be allocatable. It will be deallocated before assignment
  !!
  subroutine getIntArray_alloc_new(self,value,keyword)
    class(dictionary), intent(in)                            :: self
    integer(shortInt),dimension(:),allocatable,intent(inout) :: value
    character(*),intent(in)                                  :: keyword
    integer(shortInt)                                        :: idx
    character(100),parameter                   :: Here='getIntArray_alloc (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    if(allocated(value)) deallocate(value)

    select case (self % entries(idx) % getType())
      case(arrInt)
        value = self % entries(idx) % int1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not integer array')

    end select

  end subroutine getIntArray_alloc_new


  !!
  !! Loads a integer rank 1 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !! Variable needs to be pointer. It will be deallocated before assignment
  !!
  subroutine getIntArray_ptr_new(self,value,keyword)
    class(dictionary), intent(in)                            :: self
    integer(shortInt),dimension(:),pointer,intent(inout)     :: value
    character(*),intent(in)                                  :: keyword
    integer(shortInt)                                        :: idx,N
    character(100),parameter                   :: Here='getIntArray_ptr (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    if(associated(value)) deallocate(value)

    select case (self % entries(idx) % getType())
      case(arrInt)
        N = size(self % entries(idx) % int1_alloc)
        allocate(value(N))
        value = self % entries(idx) % int1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not integer array')

    end select

  end subroutine getIntArray_ptr_new

  !!
  !! Reads a character rank 0 from a dictionary
  !!
  subroutine getChar_new(self,value,keyword)
    class(dictionary), intent(in)  :: self
    character(*),intent(inout)     :: value
    character(*),intent(in)        :: keyword
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getChar (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    select case (self % entries(idx) % getType())
      case(word)
        ! Check if the content character fits into value
        if( len(value) < len_trim(self % entries(idx) % char0_alloc)) then
          call fatalError(Here,'value character is to short to store content. Increase its length')
        end if

        value = self % entries(idx) % char0_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a character')

    end select

  end subroutine getChar_new

  !!
  !! Reads a character rank 1 from a dictionary
  !! Variable needs to be allocatable. It will be deallocated before assignment
  !!
  subroutine getCharArray_alloc_new(self,value,keyword)
    class(dictionary), intent(in)                             :: self
    character(*),dimension(:),allocatable,intent(inout)       :: value
    character(*),intent(in)                                   :: keyword
    integer(shortInt)                                         :: idx
    character(100),parameter                    :: Here='getCharArray_alloc (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    if(allocated(value)) deallocate(value)

    select case (self % entries(idx) % getType())
      case(arrWord)
        ! Check if the content character fits into value. Any is required becouse len_trim returns
        ! an array
        if( any(len(value) < len_trim(self % entries(idx) % char1_alloc))) then
          call fatalError(Here,'value character is to short to store content. Increase its length')
        end if

        value = self % entries(idx) % char1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a character array')

    end select

  end subroutine getCharArray_alloc_new

  !!
  !! Reads a character rank 1 from a dictionary
  !! Variable needs to be pointer. It will be deallocated before assignment
  !!
  subroutine getCharArray_ptr_new(self,value,keyword)
    class(dictionary), intent(in)                             :: self
    character(*),dimension(:),pointer,intent(inout)           :: value
    character(*),intent(in)                                   :: keyword
    integer(shortInt)                                         :: idx
    character(100),parameter                       :: Here='getCharArray_ptr (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    if(associated(value)) deallocate(value)

    select case (self % entries(idx) % getType())
      case(arrWord)
        ! Check if the content character fits into value. Any is required becouse len_trim returns
        ! an array
        if( any( len(value) < len_trim(self % entries(idx) % char1_alloc)) ) then
          call fatalError(Here,'value character is to short to store content. Increase its length')
        end if

        ! Use mold to approperiatly allocate the pointer
        allocate(value ( size(self % entries(idx) % char1_alloc) ))
        value = self % entries(idx) % char1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a character array')

    end select

  end subroutine getCharArray_ptr_new

  !!
  !! Reads a dictionary rank 0 from a dictionary
  !!
  subroutine getDict_new(self,value,keyword)
    class(dictionary), intent(in)   :: self
    class(dictionary),intent(inout) :: value
    character(*),intent(in)         :: keyword
    integer(shortInt)               :: idx
    character(100),parameter        :: Here='getDict (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    select case (self % entries(idx) % getType())
      case(nestDict)
        value = self % entries(idx) % dict0_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a dictionary')

    end select

  end subroutine getDict_new

  !!
  !! Loads a real rank 0 from a dictionary into provided variable
  !! If keyword is associated with an integer it converts it to real
  !!
  subroutine getOrDefault_real(self,value,keyword,default)
    class(dictionary), intent(in)  :: self
    real(defReal),intent(inout)    :: value
    character(*),intent(in)        :: keyword
    real(defReal),intent(in)       :: default
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getOrDefault_real (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)

    if (idx == targetNotFound) then
      value = default
      return
    end if

    select case (self % entries(idx) % getType())
      case(numReal)
        value = self % entries(idx) % real0_alloc

      case(numInt)
        value = real(self % entries(idx) % int0_alloc, defReal)

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a real or int')

    end select

  end subroutine getOrDefault_real

  !!
  !! Loads a real rank 1 from dictionary into provided variable
  !! If keyword is associated with an integer it converts it to real
  !! Variable needs to be allocatable. It will be deallocated before assignment
  !!
  subroutine getOrDefault_realArray_alloc(self,value,keyword,default)
    class(dictionary), intent(in)                        :: self
    real(defReal),dimension(:),allocatable,intent(inout) :: value
    character(*),intent(in)                              :: keyword
    real(defReal),dimension(:),intent(in)                :: default
    integer(shortInt)                                    :: idx
    character(100),parameter         :: Here='getOrDefault_realArray_allocc (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    if(allocated(value)) deallocate(value)

    if (idx == targetNotFound) then
      value = default
      return
    end if

    select case (self % entries(idx) % getType())
      case(arrReal)
        value = self % entries(idx) % real1_alloc

      case(arrInt)
        value = real(self % entries(idx) % int1_alloc, defReal)

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a real array or int array')

    end select

  end subroutine getOrDefault_realArray_alloc

  !!
  !! Loads a real rank 1 from dictionary into provided variable
  !! If keyword is associated with an integer it converts it to real
  !! Variable needs to be pointer. It will be deallocated before assignment
  !!
  subroutine getOrDefault_realArray_ptr(self,value,keyword,default)
    class(dictionary), intent(in)                      :: self
    real(defReal),dimension(:),pointer,intent(inout)   :: value
    character(*),intent(in)                            :: keyword
    real(defReal),dimension(:),intent(in)              :: default
    integer(shortInt)                                  :: idx, N
    character(100),parameter            :: Here='getOrDefault_realArray_ptr (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)

    if(associated(value)) deallocate(value)

    if (idx == targetNotFound) then
      allocate(value( size(default) ))
      value = default
      return

    end if

    select case (self % entries(idx) % getType())
      case(arrReal)
        N = size (self % entries(idx) % real1_alloc)
        allocate(value(N))
        value = self % entries(idx) % real1_alloc

      case(arrInt)
        value = real(self % entries(idx) % int1_alloc, defReal)

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a real array or int array')

    end select

  end subroutine getOrDefault_realArray_ptr

  !!
  !! Loads a integer rank 0 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !!
  subroutine getOrDefault_int(self,value,keyword,default)
    class(dictionary), intent(in)  :: self
    integer(shortInt),intent(inout):: value
    character(*),intent(in)        :: keyword
    integer(shortInt), intent(in)  :: default
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getOrDefault_int (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)

    if (idx == targetNotFound) then
      value = default
      return

    end if

    select case (self % entries(idx) % getType())
      case(numInt)
        value = self % entries(idx) % int0_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not an integer')

    end select

  end subroutine getOrDefault_int

  !!
  !! Loads a integer rank 1 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !! Variable needs to be allocatable. It will be deallocated before assignment
  !!
  subroutine getOrDefault_intArray_alloc(self,value,keyword,default)
    class(dictionary), intent(in)                            :: self
    integer(shortInt),dimension(:),allocatable,intent(inout) :: value
    character(*),intent(in)                                  :: keyword
    integer(shortInt),dimension(:),intent(in)                :: default
    integer(shortInt)                                        :: idx
    character(100),parameter           :: Here='getOrDefault_intArray_alloc (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)

    if(allocated(value)) deallocate(value)

    if (idx == targetNotFound) then
      value = default
      return
    end if

    select case (self % entries(idx) % getType())
      case(arrInt)
        value = self % entries(idx) % int1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not integer array')

    end select

  end subroutine getOrDefault_intArray_alloc


  !!
  !! Loads a integer rank 1 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !! Variable needs to be pointer. It will be deallocated before assignment
  !!
  subroutine getOrDefault_intArray_ptr(self,value,keyword,default)
    class(dictionary), intent(in)                            :: self
    integer(shortInt),dimension(:),pointer,intent(inout)     :: value
    character(*),intent(in)                                  :: keyword
    integer(shortInt),dimension(:),intent(in)                :: default
    integer(shortInt)                                        :: idx,N
    character(100),parameter             :: Here='getOrDefault_intArray_ptr (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)

    if(associated(value)) deallocate(value)

    if (idx == targetNotFound) then
      allocate(value( size(default) ))
      value = default
      return

    end if


    select case (self % entries(idx) % getType())
      case(arrInt)
        N = size(self % entries(idx) % int1_alloc)
        allocate(value(N))
        value = self % entries(idx) % int1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not integer array')

    end select

  end subroutine getOrDefault_intArray_ptr

  !!
  !! Reads a character rank 0 from a dictionary
  !!
  subroutine getOrDefault_char(self,value,keyword,default)
    class(dictionary), intent(in)  :: self
    character(*),intent(inout)     :: value
    character(*),intent(in)        :: keyword
    character(*),intent(in)        :: default
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getOrDefault_char (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)

    ! Check if default exceeds length of a value charcter
    if (len(default) > len(value) ) then
      call fatalError(Here,'Default character does not fit into value')
    end if

    if (idx == targetNotFound) then
      value = default
      return

    end if

    select case (self % entries(idx) % getType())
      case(word)
        ! Check if the content character fits into value
        if( len(value) < len_trim(self % entries(idx) % char0_alloc)) then
          call fatalError(Here,'value character is to short to store content. Increase its length')
        end if

        value = self % entries(idx) % char0_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a character')

    end select

  end subroutine getOrDefault_char

  !!
  !! Reads a character rank 1 from a dictionary
  !! Variable needs to be allocatable. It will be deallocated before assignment
  !!
  subroutine getOrDefault_charArray_alloc(self,value,keyword,default)
    class(dictionary), intent(in)                             :: self
    character(*),dimension(:),allocatable,intent(inout)       :: value
    character(*),intent(in)                                   :: keyword
    character(*),dimension(:),intent(in)                      :: default
    character(charLen),dimension(size(default))               :: loc_Char
    integer(shortInt)                                         :: idx
    character(100),parameter           :: Here='getOrDefault_charArray_alloc(dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)

    if(allocated(value)) deallocate(value)

    ! Check if provided default character array exceeeds maximum character length
    if (len(default) > len(value) ) then
      call fatalError(Here,'Default character does not fit into value')
    end if

    if (idx == targetNotFound) then
      loc_Char = default
      value = loc_Char
      return
    end if

    select case (self % entries(idx) % getType())
      case(arrWord)
        ! Check if the content character fits into value. Any is required becouse len_trim returns
        ! an array
        if( any( len(value) < len_trim(self % entries(idx) % char1_alloc)) ) then
          call fatalError(Here,'value character is to short to store content. Increase its length')
        end if

        value = self % entries(idx) % char1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a character array')

    end select

  end subroutine getOrDefault_charArray_alloc

  !!
  !! Reads a character rank 1 from a dictionary
  !! Variable needs to be pointer. It will be deallocated before assignment
  !!
  subroutine getOrDefault_charArray_ptr(self,value,keyword,default)
    class(dictionary), intent(in)                             :: self
    character(*),dimension(:),pointer,intent(inout)           :: value
    character(*),intent(in)                                   :: keyword
    character(*),dimension(:),intent(in)                      :: default
    integer(shortInt)                                         :: idx
    character(100),parameter            :: Here='getOrDefault_charArray_ptr (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)

    if(associated(value)) deallocate(value)

    ! Check if provided default character array exceeeds maximum character length
    if (len(default) > len(value) ) then
      call fatalError(Here,'Default character does not fit into value')
    end if

    if (idx == targetNotFound) then
      !loc_char = default
      allocate(value(size(default)))
      value = default
      return

    end if

    select case (self % entries(idx) % getType())
      case(arrWord)
        ! Check if the content character fits into value. Any is required becouse len_trim returns
        ! an array
        if( any( len(value) < len_trim(self % entries(idx) % char1_alloc)) ) then
          call fatalError(Here,'value character is to short to store content. Increase its length')
        end if

        ! Use mold to approperiatly allocate the pointer
        allocate(value ( size(self % entries(idx) % char1_alloc) ))
        value = self % entries(idx) % char1_alloc

      case default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a character array')

    end select

  end subroutine getOrDefault_charArray_ptr


  !!
  !! Reads a real rank 0 entery from a dictionary
  !! If keyword is associated with an integer it converts it to real
  !!
  function getReal_old(self,keyword) result(value)
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

  end function getReal_old

  !!
  !! Reads a real rank 1 from a dictionary
  !! If keyword is associated with an integer it converts it to real
  !!
  function getRealArray_old(self,keyword) result(value)
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

  end function getRealArray_old

  !!
  !! Reads a integer rank 0 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !!
  function getInt_old(self,keyword) result(value)
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

  end function getInt_old

  !!
  !! Reads a integer rank 1 from a dictionary
  !! If keyword is associated with real integer (i.e. 1.0) it returns an error
  !!
  function getIntArray_old(self,keyword) result(value)
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

  end function getIntArray_old

  !!
  !! Reads a character rank 0 from a dictionary
  !!
  function getChar_old(self,keyword) result(value)
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


  end function getChar_old

  !!
  !! Reads a character rank 1 from a dictionary
  !!
  function getCharArray_old(self,keyword) result(value)
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

  end function getCharArray_old

  !!
  !! Reads a dictionary rank 0 from a dictionary
  !!
  function getDict_old(self,keyword) result(value)
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

  end function getDict_old

  !!
  !! Returns an array of all keywords associated with a real rank 0
  !!
  subroutine keysReal(self,keys)
    class(dictionary), intent(in)                            :: self
    character(nameLen),dimension(:), allocatable,intent(out) :: keys
    logical(defBool),dimension(:),allocatable                :: mask
    integer(shortInt)                                        :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == numReal)

    keys = pack(self % keywords(1:L), mask)

  end subroutine keysReal

  !!
  !! Returns an array of all keywords associated with a real rank 1
  !!
  subroutine keysRealArray(self,keys)
    class(dictionary), intent(in)                            :: self
    character(nameLen),dimension(:), allocatable,intent(out) :: keys
    logical(defBool),dimension(:),allocatable                :: mask
    integer(shortInt)                                        :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == arrReal)

    keys = pack(self % keywords(1:L), mask)

  end subroutine keysRealArray

  !!
  !! Returns an array of all keywords associated with an integer rank 0
  !!
  subroutine keysInt(self,keys)
    class(dictionary), intent(in)                            :: self
    character(nameLen),dimension(:), allocatable,intent(out) :: keys
    logical(defBool),dimension(:),allocatable                :: mask
    integer(shortInt)                                        :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == numInt)

    keys = pack(self % keywords(1:L), mask)

  end subroutine keysInt

  !!
  !! Returns an array of all keywords associated with an integer rank 0
  !!
  subroutine keysIntArray(self,keys)
    class(dictionary), intent(in)                            :: self
    character(nameLen),dimension(:), allocatable,intent(out) :: keys
    logical(defBool),dimension(:),allocatable                :: mask
    integer(shortInt)                                        :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == arrInt)

    keys = pack(self % keywords(1:L), mask)

  end subroutine keysIntArray

  !!
  !! Returns an array of all keywords associated with an character rank 0
  !!
  subroutine keysChar(self,keys)
    class(dictionary), intent(in)                            :: self
    character(nameLen),dimension(:), allocatable,intent(out) :: keys
    logical(defBool),dimension(:),allocatable                :: mask
    integer(shortInt)                                        :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == word)

    keys = pack(self % keywords(1:L), mask)

  end subroutine keysChar

  !!
  !! Returns an array of all keywords associated with an character rank 1
  !!
  subroutine keysCharArray(self,keys)
    class(dictionary), intent(in)                            :: self
    character(nameLen),dimension(:), allocatable,intent(out) :: keys
    logical(defBool),dimension(:),allocatable                :: mask
    integer(shortInt)                                        :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == arrWord)

    keys = pack(self % keywords(1:L), mask)

  end subroutine keysCharArray

  !!
  !! Returns an array of all keywords associated with a dictionary rank 0
  !!
  subroutine keysDict(self,keys)
    class(dictionary), intent(in)                            :: self
    character(nameLen),dimension(:), allocatable,intent(out) :: keys
    logical(defBool),dimension(:),allocatable                :: mask
    integer(shortInt)                                        :: L

    L = self % dictLen
    allocate( mask(L) )

    mask = (self % entries(1:L) % getType() == nestDict)

    keys = pack(self % keywords(1:L), mask)

  end subroutine keysDict

  !!
  !! Returns an array of all keywords
  !!
  subroutine keys(self,keysArr)
    class(dictionary), intent(in)                            :: self
    character(nameLen),dimension(:), allocatable,intent(out) :: keysArr
    integer(shortInt)                                        :: L

    L       = self % dictLen
    keysArr = self % keywords(1:L)

  end subroutine keys

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
    class(dictionary), intent(in)               :: entry
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
      allocate(LHS % dict0_alloc, source = RHS % dict0_alloc )
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
