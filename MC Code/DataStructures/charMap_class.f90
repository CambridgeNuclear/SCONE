module charMap_class

  use numPrecision
  use genericProcedures,  only : fatalError
  use hashFunctions_func, only : FNV_1, knuthHash

  implicit none
  private

  !! Local parameters
  integer(shortInt),parameter :: EMPTY    = 0
  integer(shortInt),parameter :: TAKEN    = 1
  integer(shortInt),parameter :: DELETED  = 2
  real(defReal),parameter     :: MAX_LOAD = 0.6

  !!
  !! Helper structure to group key with its hash and value
  !! Also contains status of the entry
  !!
  type, private :: content
    integer(shortInt)  :: status = EMPTY
    character(nameLen) :: key
    integer(shortInt)  :: hash
    integer(shortInt)  :: val
  end type


  !!
  !! Maps characters 'nameLen' long to shortInts
  !! Dictionary with character keys and integer content
  !! Implemented as hash table with open adressing
  !! Implementation is based on intMap
  !!
  type, public :: charMap
    private
    integer(shortInt)                        :: Nexp      ! Nexp = log2(N)
    integer(shortInt)                        :: N     = 0 ! Current size of the map is a power of 2
    integer(shortInt)                        :: Load  = 0 ! Number of occupied entries
    type(content), dimension(:), allocatable :: map
  contains
    procedure :: init
    procedure :: add
    procedure :: get
    procedure :: getOrDefault
    procedure :: kill
    procedure :: length

    ! Private Procedures
    procedure, private :: grow
  end type charMap

contains

  !!
  !! Initialise charMap with a desired size
  !! It will be allocated to size equal to smallest larger power of 2 > N
  !! Size will never be smaller then 8
  !!
  subroutine init(self,N)
    class(charMap), intent(inout) :: self
    integer(shortInt), intent(in) :: N
    integer(shortInt)             :: N_bar
    integer(shortInt)             :: nextPow2
    character(100), parameter :: Here = 'init (charMap_class.f90)'

    if( N <= 0) call fatalError(Here,'Size needs to be +ve')

    ! Clean map
    call self % kill()

    ! Find minumim size required for N entries
    N_bar = int(N/MAX_LOAD)

    ! Find next power of 2
    ! We are acounting for zero by using N-1 in leadz
    nextPow2 = bit_size(N_bar) - leadz(N_bar-1)

    ! Assign size to 8 if small value is provided
    if (nextPow2 < 3 ) nextPow2 = 3

    ! Allocate storage space
    self % Load = 0
    self % Nexp = nextPow2
    self % N    = 2**nextPow2
    allocate(self % map( self % N))

    ! SET map keys to EMTPY
    self % map % status = EMPTY

  end subroutine init

  !!
  !! Kill charMap
  !!
  subroutine kill(self)
    class(charMap), intent(inout) :: self

    ! Deallocate space
    if(allocated(self % map)) deallocate(self % map)

    ! Set default values of parameters
    self % N    = 0
    self % Load = 0

  end subroutine kill

  !!
  !! Returns current length of the map
  !!
  pure function length(self) result(L)
    class(charMap), intent(in) :: self
    integer(shortInt)          :: L

    L = self % Load

  end function length

  !!
  !! Add a new entry to map or overwrite the existing one
  !!
  subroutine add(self,key,val)
    class(charMap), intent(inout)  :: self
    character(nameLen), intent(in) :: key
    integer(shortInt), intent(in)  :: val
    integer(shortInt)              :: hash
    integer(shortInt)              :: idx
    logical(defBool)               :: sameEntry

    ! Check for initialisation or growth
    if (self % N == 0) then ! Initialise map
      call self % init(1)

    else if (real(self % Load +1) / self % N > MAX_LOAD) then ! Double storage space
      call self % grow()

    end if

    ! Hash character key to integer
    call FNV_1(key,hash)

    ! Hash Hashed key to map to index
    idx = knuthHash(hash, self % Nexp) + 1

    ! Find next non empty place
    sameEntry = .false.
    do while (self % map(idx) % status /= EMPTY)
      ! Exit if the entry with the same key was found
      ! Compare hash before character for speed.
      ! Sorry for nested ifs
      if( self % map(idx) % hash == hash) then
        if(self % map(idx) % key == key) then
          sameEntry = .true.
          exit

        end if
      end if

      ! Increment position
      idx = idx + 1
      ! Go to the beggining of table if overflow
      if(idx > self % N) idx = 1

    end do

    ! Increment load if not overwriting existing entry
    if (.not.sameEntry) self % Load = self % Load +1

    ! Load key and value
    self % map(idx) % key    = key
    self % map(idx) % hash   = hash
    self % map(idx) % val    = val
    self % map(idx) % status = TAKEN

  end subroutine add

  !!
  !! Retrive entry in the map
  !! Returns error if key is not present
  !!
  function get(self,key) result(val)
    class(charMap), intent(in)     :: self
    character(nameLen), intent(in) :: key
    integer(shortInt)              :: val
    integer(shortInt)              :: hash
    integer(shortInt)              :: idx
    character(100), parameter :: Here = 'get (charMap_class.f90)'

    ! Give error if map is uninitialised (empty)
    if( .not.allocated(self % map)) then
      call fatalError(Here,'Target key: '// key // ' was not found. Map is empty')
    end if

    ! Hash character key to integer
    call FNV_1(key,hash)

    ! Hash Hashed key to map to index
    idx = knuthHash(hash, self % Nexp) + 1

    ! Look for the entry
    do while (self % map(idx) % status /= EMPTY)
      ! Exit if the entry with the same key was found
      ! Compare hash before character for speed.
      ! Sorry for nested ifs
      if( self % map(idx) % hash == hash) then
        if( self % map(idx) % key == key) then
          val = self % map(idx) % val
          return
        end if
      end if

      ! Increment position
      idx = idx + 1
      ! Go to the beggining of table if overflow
      if(idx > self % N) idx = 1

    end do

    call fatalError(Here,'Target key: '// trim(key) // ' was not found')

    ! Avoid compiler warning
    val = 0

  end function get

  !!
  !! Retrive entry in the map
  !! If key is not present return val == default
  !!
  pure function getOrDefault(self,key,default) result(val)
    class(charMap), intent(in)     :: self
    character(nameLen), intent(in) :: key
    integer(shortInt), intent(in)  :: default
    integer(shortInt)              :: val
    integer(shortInt)              :: hash
    integer(shortInt)              :: idx
    character(100), parameter :: Here = 'get (charMap_class.f90)'

    ! Give default if map is uninitialised (empty)
    if( .not.allocated(self % map)) then
      val = default
      return
    end if

    ! Hash character key to integer
    call FNV_1(key,hash)

    ! Hash Hashed key to map to index
    idx = knuthHash(hash, self % Nexp) + 1

    ! Look for the entry
    do while (self % map(idx) % status /= EMPTY)
      ! Exit if the entry with the same key was found
      ! Compare hash before character for speed.
      ! Sorry for nested ifs
      if( self % map(idx) % hash == hash) then
        if( self % map(idx) % key == key) then
          val = self % map(idx) % val
          return
        end if
      end if

      ! Increment position
      idx = idx + 1
      ! Go to the beggining of table if overflow
      if(idx > self % N) idx = 1

    end do

    ! Entry was not found -> return default
    val = default

  end function getOrDefault

  !!
  !! Increase size of the map by factor of 2
  !!
  subroutine grow(self)
    class(charMap), intent(inout) :: self
    type(charMap)                 :: tempMap
    integer(shortInt)             :: i

    ! When growing rehasing is required. We will just create new map and copy all values in
    call tempMap % init( self % Load * 2)

    ! Loop throuth current table and rehash non-empty entries
    do i=1,self % N
      associate (entry => self % map(i) )
        if(entry % status == TAKEN) then
          call tempMap % add( entry % key, entry % val)

        end if
      end associate
    end do

    ! Deallocate current storage and load new storage
    self % Nexp = tempMap % Nexp
    self % N    = tempMap % N
    call move_alloc(tempMap % map, self % map)

  end subroutine grow
    
end module charMap_class
