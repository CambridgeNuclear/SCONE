module maps_class

  use numPrecision
  use genericProcedures , only : fatalError, numToChar
  use hashFunctions_func, only : knuthHash

  implicit none
  private

  !! Local parameters
  integer(shortInt) :: EMPTY    = 0
  real(defReal)     :: MAX_LOAD = 0.6

  !!
  !! Helper type to wrap key and value
  !!
  type, private :: content
    integer(shortInt) :: key
    integer(shortInt) :: val
  end type

  !!
  !! Maps integers to shortInt to shortInt
  !! Dictionary with shortInt key and shortInt content
  !! Implemented as a hash table with open adressing
  !! Does not support 0 as a key
  type, public :: intMap
    private
    integer(shortInt)                        :: Nexp
    integer(shortInt)                        :: N
    integer(shortInt)                        :: Load
    type(content), dimension(:), allocatable :: map
  contains
    procedure :: init
    procedure :: add
    procedure :: get
    !procedure :: remove
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
    character(100), parameter :: Here = 'init (maps_class.f90)'

    if( N <= 0) call fatalError(Here,'Size needs to be +ve')

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

  end subroutine init

  !!
  !! Add a new entry to map or overwrite the existing one
  !!
  subroutine add(self,key,val)
    class(intMap), intent(inout)  :: self
    integer(shortInt), intent(in) :: key
    integer(shortInt), intent(in) :: val
    integer(shortInt)             :: hash
    character(100), parameter :: Here =' add (maps_class.f90)'

    if (key == 0) call fatalError(Here,'Key cannot be 0')

    ! Increase size if maximum load factor is exceeded
    if (real(self % Load +1) / self % N > MAX_LOAD) then
      call self % grow()
    end if

    ! Calculate Hash
    hash = knuthHash(key, self % Nexp) + 1

    ! Find next non empty place
    do while (self % map(hash) % key /= EMPTY)
      ! Exit if the entry with the same key was found
      if( self % map(hash) % key == key) exit

      ! Increment position
      hash = hash + 1
      ! Go to the beggining of table if overflow
      if(hash > self % N) hash = 1

    end do

    ! Increment load if not overwriting existing entry
    if (self % map(hash) % key /= key) self % Load = self % Load +1

    ! Load key and value
    self % map(hash) % key = key
    self % map(hash) % val = val

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
    character(100), parameter :: Here = 'get (maps_class.f90)'

    ! Calculate Hash
    hash = knuthHash(key, self % Nexp) + 1

    ! Look for the entry
    do while (self % map(hash) % key /= EMPTY)
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

  end function get

  !!
  !! Increase size of the map by factor of 2
  !!
  subroutine grow(self)
    class(intMap), intent(inout) :: self
    type(intMap)                 :: tempMap
    integer(shortInt)            :: i

    ! When growing rehasing is required. We will just create new map and copy all values in
    call tempMap % init( self % N * 2)

    ! Loop throuth current table and rehash non-empty entries
    do i=1,self % N
      associate (entry => self % map(i) )
        if(entry % key /= EMPTY) then
          call tempMap % add( entry % key, entry % val)

        end if
      end associate
    end do

    ! Deallocate current storage and load new storage
    self % Nexp = tempMap % Nexp
    self % N    = tempMap % N
    call move_alloc(tempMap % map, self % map)

  end subroutine grow

    
end module maps_class
