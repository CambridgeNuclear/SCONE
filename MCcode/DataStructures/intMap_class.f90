module intMap_class

  use numPrecision
  use genericProcedures , only : fatalError, numToChar
  use hashFunctions_func, only : knuthHash

  implicit none
  private

  !! Local parameters
  integer(shortInt),parameter :: EMPTY     = 0
  integer(shortInt),parameter :: TAKEN     = 1
  integer(shortInt),parameter :: DELETED   = 2
  integer(shortInt),parameter :: END_TOKEN = -789
  real(defReal),parameter     :: MAX_LOAD  = 0.6

  !!
  !! Helper type to wrap key and value
  !! Also contains status of the entry
  !!
  type,private :: content
    integer(shortInt) :: status = EMPTY
    integer(shortInt) :: key = 0
    integer(shortInt) :: val
  end type

  !!
  !! Maps integers to shortInt to shortInt
  !!
  !! Dictionary with shortInt key and shortInt content
  !! Implemented as a hash table with open adressing
  !!
  !! NOTE: Following structure can be used to loop over entire map
  !! i = map % begin()
  !! do while (i == map % end())
  !!   ! Access value with: map % atVal(i)
  !!   ! Access key with: map % atKey(i)
  !!   i = map % next(i)
  !! end do
  !!
  !! Private members:
  !!   Nexp -> Exponent of 2 for a size of the map = log2(N)
  !!   N    -> Current size of the space avalible
  !!   Load -> Number of entries in the map
  !!   map  -> Array of content
  !!
  !! Interface:
  !!   init         -> Sets initial size. Use to pre-allocate memory
  !!   add          -> Add new entry to the map
  !!   get          -> Get entry from map. Error if not present
  !!   getOrDefault -> Get entry from map. Return default if not present
  !!   del          -> Remove an entry from the map
  !!   length       -> Return number of entries in the map
  !!   begin        -> Return index in array to first occupied element
  !!   atVal        -> Return value under index in array
  !!   atKey        -> Return key under index in array
  !!   next         -> Return index in array of the next occupied element
  !!   end          -> Return index in array for the last occupied element
  !!   kill         -> Return map to uninitialised state
  !!
  type, public :: intMap
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
    procedure :: del
    procedure :: length
    procedure :: begin
    procedure :: atVal
    procedure :: atKey
    procedure :: next
    procedure :: end
    procedure :: kill

    ! Private procedures
    procedure, private :: grow
  end type intMap


contains

  !!
  !! Initialise intMap with a desired size
  !!
  !! It will be allocated to size equal to smallest larger power of 2 > N
  !! Size will never be smaller then 8
  !!
  !! Args:
  !!   N [in] -> Desired size of the map
  !!
  !! Errors:
  !!   FatalError if N is 0 or -ve
  !!
  subroutine init(self, N)
    class(intMap), intent(inout) :: self
    integer(shortInt), intent(in):: N
    integer(shortInt)            :: N_bar
    integer(shortInt)            :: nextPow2
    character(100), parameter :: Here = 'init (intMap_class.f90)'

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
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(intMap), intent(inout) :: self

    ! Deallocate space
    if(allocated(self % map)) deallocate(self % map)

    ! Set default values of parameters
    self % N    = 0
    self % Load = 0

  end subroutine kill

  !!
  !! Returns number of elements in the map
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of elements ion the map (Load)
  !!
  !! Errors:
  !!   None
  !!
  pure function length(self) result(L)
    class(intMap), intent(in)  :: self
    integer(shortInt)          :: L

    L = self % Load

  end function length

  !!
  !! Add a new entry to map
  !!
  !! If key is not present -> add new key and value
  !! If key is present -> change value under key to val
  !!
  !! Args:
  !!   key [in] -> Integer key
  !!   val [in] -> Integer value to be set under key
  !!
  !! Errors:
  !!   None
  !!
  subroutine add(self, key, val)
    class(intMap), intent(inout)  :: self
    integer(shortInt), intent(in) :: key
    integer(shortInt), intent(in) :: val
    integer(shortInt)             :: hash

    ! Check for initialisation or growth
    if (self % N == 0) then ! Initialise map
      call self % init(1)

    else if (real(self % Load +1) / self % N > MAX_LOAD) then ! Double storage space
      call self % grow()

    end if

    ! Calculate Hash
    hash = knuthHash(key, self % Nexp) + 1

    ! Find next non empty place
    do while (self % map(hash) % status /= EMPTY)
      ! Exit if the entry with the same key was found
      if( self % map(hash) % key == key) exit

      ! Increment position
      hash = hash + 1
      ! Go to the beggining of table if overflow
      if(hash > self % N) hash = 1

    end do

    ! Increment load if not overwriting existing entry
    if (self % map(hash) % status /= TAKEN) self % Load = self % Load +1

    ! Load key and value
    self % map(hash) % key    = key
    self % map(hash) % val    = val
    self % map(hash) % status = TAKEN

  end subroutine add

  !!
  !! Retrive entry in the map
  !!
  !! Arg:
  !!   key [in] -> Integer key to access value
  !!
  !! Result:
  !!   Value stored under key
  !!
  !! Errors:
  !!   FatalError if key is not present
  !!
  function get(self, key) result(val)
    class(intMap), intent(in)     :: self
    integer(shortInt), intent(in) :: key
    integer(shortInt)             :: val
    integer(shortInt)             :: hash
    character(100), parameter :: Here = 'get (intMap_class.f90)'

    ! Give error if map is uninitialised (empty)
    if( self % Load == 0) then
      call fatalError(Here,'Target key: '// numToChar(key) // ' was not found. Map is empty')
    end if

    ! Calculate Hash
    hash = knuthHash(key, self % Nexp) + 1

    ! Look for the entry
    do while (self % map(hash) % status /= EMPTY)
      ! Exit if the entry with the same kay was found
      if( self % map(hash) % status == TAKEN .and. self % map(hash) % key == key) then
        val = self % map(hash) % val
        return

      end if
      ! Increment position
      hash = hash + 1
      ! Go to the beggining of table if overflow
      if(hash > self % N) hash = 1

    end do

    call fatalError(Here,'Target key: '// numToChar(key) // ' was not found')

    ! Avoid compier warning
    val = huge(val)

  end function get

  !!
  !! Retrive entry in the map with default value
  !!
  !! Arg:
  !!   key [in] -> Integer key to access value
  !!
  !! Result:
  !!   If key is present -> returns value under key
  !!   If key is not present -> returns default
  !!
  !! Errors:
  !!   None
  !!
  pure function getOrDefault(self,key,default) result(val)
    class(intMap), intent(in)     :: self
    integer(shortInt), intent(in) :: key
    integer(shortInt), intent(in) :: default
    integer(shortInt)             :: val
    integer(shortInt)             :: hash

    ! Give default if map is uninitialised (empty)
    if( .not.allocated(self % map)) then
      val = default
      return
    end if

    ! Calculate Hash
    hash = knuthHash(key, self % Nexp) + 1

    ! Look for the entry
    do while (self % map(hash) % status /= EMPTY)
      ! Exit if the entry with the same kay was found
      if( self % map(hash) % status == TAKEN .and. self % map(hash) % key == key) then
        val = self % map(hash) % val
        return

      end if
      ! Increment position
      hash = hash + 1
      ! Go to the beggining of table if overflow
      if(hash > self % N) hash = 1

    end do

    ! Entry was not found -> return default
    val = default

  end function getOrDefault

  !!
  !! Delete element from the map
  !!
  !! If key is not present in the dictionary nothing happens
  !!
  !! Args:
  !!   key [in] -> Key to be deleted
  !!
  !! Errors:
  !!   None
  !!
  subroutine del(self, key)
    class(intMap), intent(inout)  :: self
    integer(shortInt), intent(in) :: key
    integer(shortInt)             :: hash

    ! Quit if map is empty
    if(self % load == 0) return

    ! Calculate Hash
    hash = knuthHash(key, self % Nexp) + 1

    ! Look for the entry
    do while (self % map(hash) % status == TAKEN)
      ! Exit if the entry with the same kay was found
      if( self % map(hash) % key == key) then
        self % map(hash) % status = DELETED
        return

      end if
      ! Increment position
      hash = hash + 1
      ! Go to the beggining of table if overflow
      if(hash > self % N) hash = 1

    end do

  end subroutine del

  !!
  !! Return index in array of the first Occupied element
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Returns index END_TOKEN if map is empty
  !!
  !! Error:
  !!   None
  !!
  pure function begin(self) result(idx)
    class(intMap), intent(in) :: self
    integer(shortInt)         :: idx

    ! Return 0 if map is empty
    if (self % load == 0) then
      idx = END_TOKEN
      return
    end if

    ! First non-empty element
    do idx = 1, self % N
      if( self % map(idx) % status == TAKEN) return
    end do

    ! Should never be executed
    idx = END_TOKEN

  end function begin

  !!
  !! Return VALUE under an index in array
  !!
  !! Args:
  !!   idx [in] -> Index to access data in the array
  !!
  !! Result:
  !!   Value stored in the array under index idx
  !!
  !! Errors:
  !!   fatalError if element under idx is not TAKEN (EMPTY, DELETED or out of bounds)
  !!
  function atVal(self, idx) result (val)
    class(intMap), intent(in)     :: self
    integer(shortInt), intent(in) :: idx
    integer(shortInt)             :: val
    character(100), parameter :: Here = 'atVal (intMap_class.f90)'

    ! Check bounds
    if (idx <= 0 .or. idx > self % N) then
      call fatalError(Here, "Index is outside of bounds or map is uninitialised:" // numToChar(idx))

    else if ( self % map(idx) % status /= TAKEN) then
      call fatalError(Here, "Index refers to unoccupied entry:" // numToChar(idx))
    end if

    ! Return value
    val = self % map(idx) % val

  end function atVal

  !!
  !! Return KEY under an index in the array
  !!
  !! Args:
  !!   idx [in] -> Index to access data in the array
  !!
  !! Result:
  !!   KEY of entry at index idx in the array
  !!
  !! Errors:
  !!   fatalError if element under idx is not TAKEN (EMPTY, DELETED or out of bounds)
  !!
  function atKey(self, idx) result(key)
    class(intMap), intent(in)     :: self
    integer(shortInt), intent(in) :: idx
    integer(shortInt)             :: key
    character(100), parameter :: Here = 'atKey (intMap_class.f90)'

    ! Check bounds and status
    if (idx <= 0 .or. idx > self % N) then
      call fatalError(Here, "Index is outside of bounds or map is uninitialised:" // numToChar(idx))

    else if ( self % map(idx) % status /= TAKEN) then
      call fatalError(Here, "Index refers to unoccupied entry:" // numToChar(idx))
    end if

    ! Return key
    key = self % map(idx) % key

  end function atKey

  !!
  !! Find next TAKEN element following index idx
  !!
  !! Args:
  !!   idx [in] -> Index in array to start search
  !!
  !! Result:
  !!   Index of the next TAKEN element in the array
  !!   If there are no more TAKEN elements before the end of array
  !!     return END_TOKEN
  !!
  !! Errors:
  !!   fatalError if idx is out-of-bounds
  !!
  function next(self, idx) result(next_idx)
    class(intMap), intent(in)     :: self
    integer(shortInt), intent(in) :: idx
    integer(shortInt)             :: next_idx
    character(100), parameter :: Here = 'next (intMap_class.f90)'

    ! Check bounds and status
    if (idx <= 0 .or. idx > self % N) then
      call fatalError(Here, "Index is outside of bounds or map is uninitialised:" // numToChar(idx))
    end if

    ! Loop until the next element
    do next_idx = idx +1, self % N
      if(self % map(next_idx) % status == TAKEN) return
    end do

    ! Reached the end of array
    next_idx = END_TOKEN

  end function next

  !!
  !! Return END_TOKEN
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   END_TOKEN parameter
  !!
  !! Errors:
  !!   None
  !!
  pure function end(self) result(idx)
    class(intMap), intent(in) :: self
    integer(shortInt)         :: idx

    idx = END_TOKEN

  end function end

  !!
  !! Increase size of the map by factor of 2
  !!
  !! Doubles the size of the map
  !! Utility function used by the map
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  subroutine grow(self)
    class(intMap), intent(inout) :: self
    type(intMap)                 :: tempMap
    integer(shortInt)            :: i

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

    
end module intMap_class
