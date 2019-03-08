module intMap_class

  use numPrecision
  use genericProcedures , only : fatalError, numToChar
  use hashFunctions_func, only : knuthHash

  implicit none
  private

  !! Local parameters
  integer(shortInt),parameter :: EMPTY    = 0
  integer(shortInt),parameter :: TAKEN    = 1
  integer(shortInt),parameter :: DELETED  = 2
  real(defReal),parameter     :: MAX_LOAD = 0.6

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
  !! Dictionary with shortInt key and shortInt content
  !! Implemented as a hash table with open adressing
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
    procedure :: kill
    procedure :: length

    ! Private procedures
    procedure, private :: grow
  end type intMap


contains

  !!
  !! Initialise intMap with a desired size
  !! It will be allocated to size equal to smallest larger power of 2 > N
  !! Size will never be smaller then 8
  !!
  subroutine init(self,N)
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
  !! Kill intMap
  !!
  subroutine kill(self)
    class(intMap), intent(inout) :: self

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
    class(intMap), intent(in)  :: self
    integer(shortInt)          :: L

    L = self % Load

  end function length

  !!
  !! Add a new entry to map or overwrite the existing one
  !!
  subroutine add(self,key,val)
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
  !! Returns error if key is not present
  !!
  function get(self,key) result(val)
    class(intMap), intent(in)     :: self
    integer(shortInt), intent(in) :: key
    integer(shortInt)             :: val
    integer(shortInt)             :: hash
    character(100), parameter :: Here = 'get (intMap_class.f90)'

    ! Give error if map is uninitialised (empty)
    if( .not.allocated(self % map)) then
      call fatalError(Here,'Target key: '// numToChar(key) // ' was not found. Map is empty')
    end if

    ! Calculate Hash
    hash = knuthHash(key, self % Nexp) + 1

    ! Look for the entry
    do while (self % map(hash) % status /= EMPTY)
      ! Exit if the entry with the same kay was found
      if( self % map(hash) % key == key) then
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
  !! Returns entry under key.
  !! If key is not present return val == default
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
      if( self % map(hash) % key == key) then
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
  !! Increase size of the map by factor of 2
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
