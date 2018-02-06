module dictionary_class

  use numPrecision
  use genericProcedures, only: linFind, searchError, fatalError

  implicit none
  private
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
  ! 1) Add put_<type>   subroutine to dictDontent class
  ! 2) Add another entery for <type> to select type in deepCopy_dictCont subroutine
  ! 3) Add store_<type> subroutine to dictionary class
  ! 4) Add get<type> function to dictionary class
  !
  ! To increase maximum rank would be messy. Don't do it! If you are determined you can probably
  ! hack existing code by storing higher rank arrays in rank1 array. Their structure should then
  ! be stored seperatly in integers in dictContent class.
  !
  ! At the moment there are no functions to return all keywords associated with a given type inside .
  ! the dictionary. It can be easly implemented if needed by adding an integer variable to
  ! dictContent and a function to retrive it. Then integer values for all enteries could be pulled
  ! and use to mask keywords.
  !

  integer(shortInt),parameter,public :: charLen  = nameLen
  integer(shortInt),parameter        :: empty    = 0
  integer(shortInt),parameter        :: numInt   = 1
  integer(shortInt),parameter        :: numReal  = 2
  integer(shortInt),parameter        :: word     = 3
  integer(shortInt),parameter        :: nestDict = 4
  integer(shortInt),parameter        :: arrInt   = 5
  integer(shortInt),parameter        :: arrReal  = 6
  integer(shortInt),parameter        :: arrWord  = 7

  type,public :: dictContent
    private
    class(*),pointer               :: rank0_ptr => null()
    class(*),dimension(:),pointer  :: rank1_ptr => null()
    integer(shortInt)              :: type      = empty
  contains
    generic   :: put  => put_real      ,&
                         put_realArray ,&
                         put_int       ,&
                         put_intArray  ,&
                         put_char      ,&
                         put_charArray ,&
                         put_dict
    procedure :: pull_rank0
    procedure :: pull_rank1
    procedure :: kill => kill_dictCont

    generic   :: assignment(=) => shallowCopy
    procedure :: deepCopy      => deepCopy_dictCont
    procedure :: shallowCopy   => shallowCopy_dictCont
    procedure :: getType       => getType_dictContent

    procedure,private :: put_real
    procedure,private :: put_int
    procedure,private :: put_char
    procedure,private :: put_realArray
    procedure,private :: put_intArray
    procedure,private :: put_charArray
    procedure,private :: put_dict

  end type dictContent



  type, public :: dictionary
    private
    character(nameLen),dimension(:),allocatable  :: keywords
    type(dictContent),dimension(:), allocatable  :: entries
    integer(shortInt)                            :: maxSize = 0
    integer(shortInt)                            :: dictLen = 0
    integer(shortInt)                            :: stride  = 20
  contains

    generic    :: store => store_real      ,&
                           store_realArray ,&
                           store_int       ,&
                           store_intArray  ,&
                           store_char      ,&
                           store_charArray ,&
                           store_dict
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

    generic    :: assignment(=) => deepCopy
    procedure  :: extendBy
    procedure  :: getEmptyIdx
    procedure  :: deepCopy => deepCopy_dictionary
    procedure  :: init
    procedure  :: kill => kill_dictionary

    procedure, private :: store_real
    procedure, private :: store_realArray
    procedure, private :: store_int
    procedure, private :: store_intArray
    procedure, private :: store_char
    procedure, private :: store_charArray
    procedure, private :: store_dict

  end type dictionary

contains

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
     ! call shallowCopy(entries(i), self % entries(i))
    end do

    deallocate( self % keywords)
    deallocate( self % entries)

    call move_alloc(keywords, self % keywords)
    call move_alloc(entries , self % entries )

    self % maxSize = newSize

    !*** Debug ***!
   ! print *, ' Extending'
    ! ******************* !
  end subroutine extendBy



!  subroutine deepCopy_dict(LHS,RHS)
!    class(dictionary),intent(out) :: LHS
!    type(dictionary), intent(in)  :: RHS
!    integer(shortInt)             :: i
!    character(nameLen)            :: locKeyword
!    class(*),pointer              :: temp_rank0
!    class(*),pointer              :: temp_rank1
!
!    LHS % init(RHS % maxSize, RHS % stride)
!
!    do i=1,size(RHS % keywords)
!      locKeyword = RHS % keywords(i)
!
!    end do
!
!
!  end subroutine deepCopy_dict

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

    ! Keywords can allocate with some garbage inside. Make shure all enteries are blank.
    self % keywords = ''

    print *, self % keywords

    if (present(stride)) self % stride = stride

  end subroutine init

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

      ! Loop below could extend only to self % dictLen but lets make double shure we kill all data.
      ! There is no need to optimise this and it is more robust this way
      do i=1,size(self % keywords)
        call self % entries(i) % kill()
      end do

      deallocate(self % keywords)
      deallocate(self % entries)

      self % maxSize = 0
      self % dictLen = 0
    end if

  end subroutine kill_dictionary



  subroutine deepCopy_dictionary(LHS,RHS)
    class(dictionary), intent(inout) :: LHS
    class(dictionary), intent(in)  :: RHS
    integer(shortInt)              :: rhsSize, stride
    integer(shortInt)              :: i

    call LHS % kill()

    rhsSize = size( RHS % keywords)
    stride = RHS % stride

    call LHS % init(rhsSize,stride)

    LHS % keywords = RHS % keywords

    do i=1,RHS % dictLen
      call LHS % entries(i) % deepCopy( RHS % entries(i))
    end do

    LHS % dictLen = RHS % dictLen
    LHS % maxSize = RHS % maxSize
    LHS % stride  = RHS % stride

  end subroutine deepCopy_dictionary


  function getReal(self,keyword) result(value)
    class(dictionary), intent(in)  :: self
    character(*),intent(in)        :: keyword
    real(defReal)                  :: value
    class(*),pointer               :: temp_ptr
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getReal (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    call self % entries(idx) % pull_rank0(temp_ptr)

    select type(temp_ptr)
      type is (real(defReal))
        value = temp_ptr

      class default
        call fatalError(Here,'Value under keyword ' // keyword // ' is not a real')

    end select

  end function getReal



  function getRealArray(self,keyword) result(value)
    class(dictionary), intent(in)          :: self
    character(*),intent(in)                :: keyword
    real(defReal),dimension(:),allocatable :: value
    class(*),dimension(:),pointer          :: temp_ptr
    integer(shortInt)                      :: idx
    character(100),parameter               :: Here='getRealArray (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    call self % entries(idx) % pull_rank1(temp_ptr)

    select type(temp_ptr)
      type is (real(defReal))
        value = temp_ptr

      class default
        call fatalError(Here,'Value under keyword ' // keyword // ' is not a real')

    end select

  end function getRealArray



  function getInt(self,keyword) result(value)
    class(dictionary), intent(in)  :: self
    character(*),intent(in)        :: keyword
    integer(shortInt)              :: value
    class(*),pointer               :: temp_ptr
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getInt (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    call self % entries(idx) % pull_rank0(temp_ptr)

    select type(temp_ptr)
      type is (integer(shortInt))
        value = temp_ptr

      class default
        call fatalError(Here,'Value under keyword ' // keyword // ' is not an integer')

    end select

  end function getInt



  function getIntArray(self,keyword) result(value)
    class(dictionary), intent(in)              :: self
    character(*),intent(in)                    :: keyword
    integer(shortInt),dimension(:),allocatable :: value
    class(*),dimension(:),pointer              :: temp_ptr
    integer(shortInt)                          :: idx
    character(100),parameter                   :: Here='getIntArray (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    call self % entries(idx) % pull_rank1(temp_ptr)

    select type(temp_ptr)
      type is (integer(shortInt))
        value = temp_ptr

      class default
        call fatalError(Here,'Value under keyword ' // keyword // ' is not an integer')

    end select

  end function getIntArray



  function getChar(self,keyword) result(value)
    class(dictionary), intent(in)  :: self
    character(*),intent(in)        :: keyword
    character(charLen)             :: value
    class(*),pointer               :: temp_ptr
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getChar (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    call self % entries(idx) % pull_rank0(temp_ptr)

    select type(temp_ptr)
      type is (character(*))
        value = temp_ptr

      class default
        call fatalError(Here,'Value under keyword ' // keyword // ' is not a character')

    end select


  end function getChar



  function getCharArray(self,keyword) result(value)
    class(dictionary), intent(in)               :: self
    character(*),intent(in)                     :: keyword
    character(charLen),dimension(:),allocatable :: value
    character(charLen),dimension(:),pointer     :: localPointer
    class(*),dimension(:),pointer               :: temp_ptr
    integer(shortInt)                           :: idx
    character(100),parameter                    :: Here='getCharArray (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    call self % entries(idx) % pull_rank1(temp_ptr)

    select type(temp_ptr)
      type is (character(*))
        !**** Warning ****
        ! For some reason "localPointer" and this syntax is required to associate local pointer
        ! with class(*) pointer from container. It seems as if pointer to character array was
        ! like an array of pointers.
        ! Rank Remapping plays a role here -> investigate further
        localPointer(1:size(temp_ptr)) => temp_ptr(1:size(temp_ptr))

      class default
        call fatalError(Here,'Array under keyword ' // keyword // ' is not a character array')

    end select

    value = localPointer

  end function getCharArray


  function getDict(self,keyword) result(value)
    class(dictionary), intent(in)  :: self
    character(*),intent(in)        :: keyword
    type(dictionary)               :: value
    class(*),pointer               :: temp_ptr
    integer(shortInt)              :: idx
    character(100),parameter       :: Here='getChar (dictionary_class.f90)'

    idx = linFind(self % keywords, keyword)
    call searchError(idx,Here)

    call self % entries(idx) % pull_rank0(temp_ptr)

    select type(temp_ptr)
      type is (dictionary)
        value = temp_ptr

      class default
        call fatalError(Here,'Entery under keyword ' // keyword // ' is not a dictionary')

    end select

  end function getDict


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



  subroutine store_real(self,keywordArgument,entry)
    class(dictionary), intent(inout)  :: self
    character(*), intent(in)          :: keywordArgument
    real(defReal), intent(in)         :: entry
    character(nameLen)                :: keyword
    integer(shortInt)                 :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    call self % entries(idx) % put(entry)

  end subroutine store_real



  subroutine store_realArray(self,keywordArgument,entry)
    class(dictionary), intent(inout)        :: self
    character(*), intent(in)                :: keywordArgument
    real(defReal), dimension(:), intent(in) :: entry
    character(nameLen)                      :: keyword
    integer(shortInt)                       :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    call self % entries(idx) % put(entry)

  end subroutine store_realArray



  subroutine store_int(self,keywordArgument,entry)
    class(dictionary), intent(inout)  :: self
    character(*), intent(in)          :: keywordArgument
    integer(shortInt), intent(in)     :: entry
    character(nameLen)                :: keyword
    integer(shortInt)                 :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    call self % entries(idx) % put(entry)

  end subroutine store_int



  subroutine store_intArray(self,keywordArgument,entry)
    class(dictionary), intent(inout)            :: self
    character(*), intent(in)                    :: keywordArgument
    integer(shortInt), dimension(:), intent(in) :: entry
    character(nameLen)                          :: keyword
    integer(shortInt)                           :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    call self % entries(idx) % put(entry)

  end subroutine store_intArray



  subroutine store_char(self,keywordArgument,entry)
    class(dictionary), intent(inout)  :: self
    character(*), intent(in)          :: keywordArgument
    character(*), intent(in)          :: entry
    character(nameLen)                :: keyword
    integer(shortInt)                 :: idx

    keyword   = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    call self % entries(idx) % put(entry)

  end subroutine store_char



  subroutine store_charArray(self,keywordArgument,entry)
    class(dictionary), intent(inout)            :: self
    character(*), intent(in)                    :: keywordArgument
    character(*), dimension(:), intent(in)      :: entry
    character(nameLen)                          :: keyword
    integer(shortInt)                           :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    call self % entries(idx) % put(entry)

  end subroutine store_charArray



  subroutine store_dict(self,keywordArgument,entry)
    class(dictionary), intent(inout)            :: self
    character(*), intent(in)                    :: keywordArgument
    type(dictionary), intent(in)                :: entry
    character(nameLen)                          :: keyword
    integer(shortInt)                           :: idx

    keyword = keywordArgument

    idx = self % getEmptyIdx(keyword)

    self % keywords(idx)  = keyword
    call self % entries(idx) % put(entry)

  end subroutine store_dict



  subroutine shallowCopy_dictCont(LHS,RHS)
    class(dictContent), intent(out)  :: LHS
    type(dictContent), intent(in)    :: RHS

    LHS % rank0_ptr => RHS % rank0_ptr
    LHS % rank1_ptr => RHS % rank1_ptr

  end subroutine shallowCopy_dictCont

  subroutine deepCopy_dictCont(LHS,RHS)
    class(dictContent), intent(out)  :: LHS
    type(dictContent), intent(in)    :: RHS
    class(*), pointer                :: rank0         ! Local pointer to rank0 object
    class(*), dimension(:),pointer   :: rank1         ! Local pointer to rank1 object
    real(defReal),pointer            :: real_ptr
    integer(shortInt),pointer        :: int_ptr
    character(charLen),pointer       :: char_ptr
    type(dictionary),pointer         :: dict_ptr

    real(defReal),dimension(:),pointer :: realArray_ptr
    integer(shortInt),dimension(:),pointer :: intArray_ptr
    character(charLen),dimension(:),pointer :: charArray_ptr
    character(charLen),dimension(:),allocatable :: charArray_alloc

    ! Copy rank 0 class(*) polymorphic pointer
    if(associated(RHS % rank0_ptr)) then
      rank0 => RHS % rank0_ptr

      select type (rank0)
        type is( real(defReal))
          allocate(real_ptr)
          real_ptr = rank0
          LHS % rank0_ptr => real_ptr

        type is( integer(shortInt))
          allocate(int_ptr)
          int_ptr = rank0
          LHS % rank0_ptr => int_ptr

        type is( character(*))
          allocate(char_ptr)
          char_ptr = rank0
          LHS % rank0_ptr => char_ptr

        type is( dictionary)
          allocate(dict_ptr)
          call dict_ptr % deepCopy(rank0)
          LHS % rank0_ptr => dict_ptr
        ! Extra entry for a dictionary

      end select

    end if

    ! Copy rank 1 class(*) polymorphic pointer
    if(associated(RHS % rank1_ptr)) then
      rank1 => RHS % rank1_ptr

      select type (rank1)
        type is( real(defReal))
          allocate( realArray_ptr(size(rank1)) )
          realArray_ptr = rank1
          LHS % rank1_ptr => realArray_ptr

        type is( integer(shortInt))
          allocate( intArray_ptr(size(rank1)) )
          intArray_ptr = rank1
          LHS % rank1_ptr => intArray_ptr

        type is( character(*))
          ! Again slightly weird code. I am using allocatable character array to temporary
          ! store the data. I am quite confused about how character arrays (effectivly array of
          ! arrays in Fortran! [I think...]) interact with pointers. The code below works OK
          ! nevertheless.

          charArray_ptr(1:size(rank1)) => rank1(1:size(rank1)) ! Obtain ptr to character array
          charArray_alloc = charArray_ptr                      ! Store character array in local copy

          allocate( charArray_ptr(size(rank1)))      ! Allocate new space under the ptr
          charArray_ptr = charArray_alloc            ! Copy values from local copy to final memory
          LHS % rank1_ptr => charArray_ptr           ! Attach pointer to the new memory

      end select
    end if

  end subroutine deepCopy_dictCont



  subroutine kill_dictCont(self)
    class(dictContent), intent(inout) :: self

    if(associated(self % rank0_ptr)) deallocate (self % rank0_ptr)
    if(associated(self % rank1_ptr)) deallocate (self % rank1_ptr)

  end subroutine kill_dictCont



  subroutine pull_rank0(self,polymorph_ptr)
    class(dictContent), intent(in)   :: self
    class(*),pointer, intent(inout)  :: polymorph_ptr
    character(100),parameter         :: Here='pull_rank0 (dictionary_class.f90)'

    if (associated(self % rank0_ptr)) then
      polymorph_ptr => self % rank0_ptr

    elseif (associated(self % rank1_ptr)) then
      call fatalError(Here,'Content under requested keyword is of rank 1')

    else
      call fatalError(Here,'Attempting to retrive data from empty dictContent')

    end if

  end subroutine pull_rank0



  subroutine pull_rank1(self,polymorph_ptr)
    class(dictContent), intent(in)                 :: self
    class(*),pointer, dimension(:), intent(inout)  :: polymorph_ptr
    character(100),parameter                       :: Here='pull_rank1 (dictionary_class.f90)'

    if (associated(self % rank0_ptr)) then
      call fatalError(Here,'Content under requested keyword is of rank 0')

    elseif (associated(self % rank1_ptr)) then
      polymorph_ptr => self % rank1_ptr

    else
      call fatalError(Here,'Attempting to retrive data from empty dictContent')

    end if


  end subroutine pull_rank1



  subroutine put_real(self,input)
    class(dictContent), intent(inout) :: self
    real(defReal), intent(in)         :: input
    real(defReal), pointer            :: newMemory
!    character(100),parameter          :: Here='put_real (dictionary_class.f90)'

    if(associated(self % rank0_ptr)) deallocate (self % rank0_ptr)
    if(associated(self % rank1_ptr)) deallocate (self % rank1_ptr)

    allocate(newMemory)
    newMemory = input

    self % rank0_ptr => newMemory
    self % type  =  numReal

  end subroutine put_real



  subroutine put_realArray(self,input)
    class(dictContent), intent(inout)      :: self
    real(defReal),dimension(:), intent(in) :: input
    real(defReal),dimension(:),pointer     :: newMemory
    integer(shortInt)                      :: inputSize

    if(associated(self % rank0_ptr)) deallocate (self % rank0_ptr)
    if(associated(self % rank1_ptr)) deallocate (self % rank1_ptr)

    inputSize = size(input)
    allocate( newMemory(inputSize) )
    newMemory = input

    self % rank1_ptr => newMemory
    self % type  =  arrReal

  end subroutine put_realArray



  subroutine put_int(self,input)
    class(dictContent), intent(inout) :: self
    integer(shortInt), intent(in)     :: input
    integer(shortInt), pointer        :: newMemory

    if(associated(self % rank0_ptr)) deallocate (self % rank0_ptr)
    if(associated(self % rank1_ptr)) deallocate (self % rank1_ptr)

    allocate(newMemory)
    newMemory = input

    self % rank0_ptr => newMemory
    self % type = numInt

  end subroutine put_int



  subroutine put_intArray(self,input)
    class(dictContent), intent(inout)          :: self
    integer(shortInt),dimension(:), intent(in) :: input
    integer(shortInt),dimension(:),pointer     :: newMemory
    integer(shortInt)                          :: inputSize

    if(associated(self % rank0_ptr)) deallocate (self % rank0_ptr)
    if(associated(self % rank1_ptr)) deallocate (self % rank1_ptr)

    inputSize = size(input)
    allocate( newMemory(inputSize) )
    newMemory = input

    self % rank1_ptr => newMemory
    self % type = arrInt

  end subroutine put_intArray


  subroutine put_char(self,input)
    class(dictContent), intent(inout) :: self
    character(*), intent(in)          :: input
    character(charLen), pointer       :: newMemory
    !integer(shortInt)                 :: inputLen

    if(associated(self % rank0_ptr)) deallocate (self % rank0_ptr)
    if(associated(self % rank1_ptr)) deallocate (self % rank1_ptr)

    !inputLen  = len(input)

    allocate(newMemory)
    newMemory = input

    self % rank0_ptr => newMemory
    self % type = word

  end subroutine put_char

  subroutine put_charArray(self,input)
    class(dictContent), intent(inout)        :: self
    character(*), dimension(:),intent(in)    :: input
    character(charLen),dimension(:), pointer :: newMemory
    integer(shortInt)                        :: inputSize
    integer(shortInt)                        :: inputLen

    if(associated(self % rank0_ptr)) deallocate (self % rank0_ptr)
    if(associated(self % rank1_ptr)) deallocate (self % rank1_ptr)

    inputSize = size(input)
    !inputLen = len(input)

    allocate(newMemory(inputSize) )

    newMemory = input

    self % rank1_ptr => newMemory
    self % type = arrWord

  end subroutine put_charArray



  subroutine put_dict(self,input)
    class(dictContent), intent(inout)    :: self
    type(dictionary), intent(in)         :: input
    type(dictionary), pointer            :: newMemory

    if(associated(self % rank0_ptr)) deallocate (self % rank0_ptr)
    if(associated(self % rank1_ptr)) deallocate (self % rank1_ptr)

    allocate(newMemory)
    newMemory = input     ! Hard Copy the dictionary

    self % rank0_ptr => newMemory
    self % type = nestDict

  end subroutine put_dict


  elemental function getType_dictContent(self) result(type)
    class(dictContent), intent(in)     :: self
    integer(shortInt)                  :: type

    type = self % type

  end function getType_dictContent

end module dictionary_class
