module surface_inter
  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, dotProduct, linFind
  use vector_class,       only : vector
  use hashFunctions_func, only : FNV_1

  implicit none
  private
  !
  ! Unfortunatly this module will be somewhat breaking rules
  ! In order to avoid circular dependencies surface interface, surface slot and surface shelf
  ! need to be placed in the same module. This is why this one is quite monstrous.
  !


  ! Local Paramethers
  character(*),parameter         :: defTol     = '6'
  character(*),parameter         :: expSize    = '2'
  character(*),parameter         :: width      = '13'  ! defTol + expSize + 5
  integer(shortInt), parameter   :: int_width  = 13   ! Needs to be the same as width but as int
  integer(shortInt),parameter    :: UNHASHED   = 0
  character(*),parameter,public  :: realForm   = 'ES'//width//'.'//defTol//'E'//expSize !
  integer(shortInt),parameter    :: DEF_STRIDE = 20

  ! Public functions
  public :: printSurfDef  !! Print surface definition string


  !!
  !! Abstract interface for a surface
  !!
  type, abstract, public :: surface
    private
    integer(shortInt)       :: hash   = UNHASHED
    integer(shortInt)       :: surfID = -1
    integer(shortInt)       :: idx    = 0

  contains
    ! Operators and assignments
    generic                                      :: operator(.same.) => same_surf

    ! Initialisation & Indentification procedures
    procedure                                    :: id
    procedure                                    :: setId
    procedure(type),deferred                     :: type
    procedure(getDef),deferred                   :: getDef
    procedure                                    :: cannotBeBoundary
    procedure                                    :: setBoundaryConditions
    procedure                                    :: hashSurfDef
    procedure, private                           :: same_surf

    ! Runtime procedures
    procedure                                    :: myIdx
    procedure                                    :: halfspace
    procedure                                    :: bounce
    procedure(evaluate), deferred                :: evaluate
    procedure(distance), deferred                :: distance
    procedure(normalVector), deferred            :: normalVector
    procedure                                    :: boundaryTransform

  end type surface

  !!
  !! Slot to store diffrent polymorphic surfaces in the single array
  !!
  type,public,extends(surface) :: surfaceSlot
    private
    class(surface),allocatable :: slot
  contains
    ! Loading into slot
    generic           :: assignment(=) => copy_slot
    procedure         :: moveAllocFrom => moveAllocFrom_slot
    procedure         :: load          => load_slot
    procedure,private :: copy_slot

    ! Initialisation & Identification procedures
    procedure :: id                    => id_slot
    procedure :: setId                 => setId_slot
    procedure :: type                  => type_slot
    procedure :: getDef                => getDef_slot
    procedure :: cannotBeBoundary      => cannotBeBoundary_slot
    procedure :: setBoundaryConditions => setBoundaryConditions_slot
    procedure :: hashSurfDef           => hashSurfDef_slot

    ! Run time procedures
    procedure :: evaluate          => evaluate_slot
    procedure :: bounce            => bounce_slot
    procedure :: distance          => distance_slot
    procedure :: normalVector      => normalVector_slot
    procedure :: boundaryTransform => boundaryTransform_slot

  end type surfaceSlot

  !!
  !! Stores number of polymorphic surfaces
  !! Public interface:
  !!  shelf     -> array of surface slots
  !!  addUnique -> adds a surface with unique definition and id to Shelf
  !!  getOrAdd  -> returns ID & Idx of provided surface, If it isn't present it adds it to the shelf
  !!  getIdx    -> returns Idx of surface with ID. returns "targetNotFound" if ID is not present
  !!  trimSize  -> Resizes array to contain no empty entries
  !!  init      -> initialises with maximum size and stride
  !!  kill      -> returns to uninitialised state
  !!
  type,public :: surfaceShelf
    !**private
    type(surfaceSlot),dimension(:),allocatable, public :: shelf
    integer(shortInt)                                  :: Nmax   = 0
    integer(shortInt)                                  :: N      = 0
    integer(shortInt)                                  :: stride = DEF_STRIDE
    integer(shortInt)                                  :: nextID = 1

  contains
    procedure :: init       => init_shelf
    procedure :: kill       => kill_shelf
    procedure :: trimSize   => trimSize_shelf

    procedure :: addUnique  => addUnique_shelf
    procedure :: getOrAdd   => getOrAdd_shelf
    procedure :: getIdx     => getIdx_shelf

    ! Private procedures
    procedure, private :: grow   => grow_shelf
    procedure, private :: resize => resize_shelf
    procedure, private :: freeID => freeId_shelf

  end type surfaceShelf


  abstract interface

    !!
    !! Returns surface hashed definition
    !!
    elemental function hash(self)
      import :: shortInt,&
                surface
      class(surface), intent(in) :: self
      integer(shortInt)          :: hash
    end function hash

    !!
    !! Return character with name of surface type
    !!
    elemental function type(self)
      import :: nameLen,&
                surface
      class(surface), intent(in) :: self
      character(nameLen)         :: type
    end function type

    !!
    !! Returns string containing surface def into allocatable string
    !! Structure of surf def can be diffrent for each surf type
    !!
    pure subroutine getDef(self,string)
      import :: surface
      class(surface), intent(in)             :: self
      character(:),allocatable,intent(inout) :: string
    end subroutine getDef


    !!
    !! Return a value of the surface expression
    !!
    elemental subroutine evaluate(self,res, r)
      import :: surface, &
                defReal, &
                vector, &
                surfaceShelf
      class(surface), intent(in)              :: self
      real(defReal), intent(out)              :: res
      type(vector), intent(in)                :: r
    end subroutine evaluate

    !!
    !! Return +ve distance to surface from point r along direction u
    !! Return INFINITY if there is no crossing
    !! Also return index of the surface being X-ed
    !!
    elemental subroutine distance(self, dist, idx, r, u)
      import :: surface, &
                defReal,&
                vector, &
                shortInt, &
                surfaceShelf
      class(surface), intent(in)              :: self
      real(defReal), intent(out)              :: dist
      integer(shortInt), intent(out)          :: idx
      type(vector), intent(in)                :: r
      type(vector), intent(in)                :: u
    end subroutine distance

    !!
    !! Return vector normal to the surface for a point r on the surface
    !! Vector is pointing into +ve halfspace
    !! No check if r lies on the surface is performed
    !!
    elemental function normalVector(self, r) result(normal)
      import :: surface, &
                vector, &
                surfaceShelf
      class(surface), intent(in)      :: self
      type(vector), intent(in)        :: r
      type(vector)                    :: normal
    end function normalVector

  end interface

contains

!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!!
!! surface procedures
!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!!
  !!
  !! Returns surface ID
  !!
  elemental function id(self)
    class(surface), intent(in) :: self
    integer(shortInt)          :: id

    id = self % surfID

  end function id

  !!
  !! Sets surface ID
  !!
  elemental subroutine setId(self,ID)
    class(surface), intent(inout) :: self
    integer(shortInt), intent(in) :: ID

    self % surfID = ID

  end subroutine setId

  !!
  !! Determine whether a point occupies the positive or negative halfspace of a surface
  !! Point can also be located on a surface - must include direction to determine halfspace
  !! Returns .true. for +ve Halfspace
  !!
  elemental function halfspace(self, r, u) result(position)
    class(surface), intent(in)              :: self
    type(vector),intent(in)                 :: r
    type(vector),intent(in)                 :: u
    real(defReal)                           :: res
    logical(defBool)                        :: position

    call self % evaluate(res, r)

    ! Set halfspace based on sign of res
    if( res > ZERO) then
      position = infront
    else
      position = behind
    end if

    ! If res is within surface tolerance choose halfspace based on direction
    if(abs(res) < surface_tol) then
      position = (self % normalVector(r) .dot. u) > ZERO
    end if

  end function halfspace

  !!
  !! Reflect a direction of a point r at the surface and nudge it along new direction
  !!
  elemental subroutine bounce(self, r, u)
    class(surface), intent(in)                 :: self
    type(vector),intent(inout)                 :: r
    type(vector),intent(inout)                 :: u
    type(vector)                               :: normal
    real(defReal)                              :: magSquared

    normal     = self % normalVector(r)
    magSquared = normal .dot. normal

    u = u - TWO * (u .dot. normal) * normal / magSquared
    r = r + NUDGE * u

  end subroutine bounce

  !!
  !! Always returns true.
  !! Surfaces that support boundary conditions need to overwrite this procedure
  !!
  function cannotBeBoundary(self) result(itCant)
    class(surface), intent(in) :: self
    logical(defBool)           :: itCant

    itCant = .true.

  end function cannotBeBoundary


  !!
  !! Provide Boundary Condition as array of integers
  !! Base class returns error. Surfaces that support BC need to override this subroutine
  !!
  subroutine setBoundaryConditions(self, BC)
    class(surface), intent(inout)               :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here='setBoundaryConditions (surface_inter.f90)'

    call fatalError(Here,'Surface: ' // self % type() // ' does not accept BCs')

  end subroutine setBoundaryConditions


  !!
  !! Perform BC coordinate transformation of a point r with direction u
  !! By default surfaces can have no BC.
  !! Thus by default they need to give error if the boundaryTransform is called
  !! If surface supports BC it needs to override this subroutine
  !!
  subroutine boundaryTransform(self, r, u)
    class(surface), intent(in)     :: self
    type(vector), intent(inout)    :: r
    type(vector), intent(inout)    :: u
    character(100),parameter :: Here = 'boundaryTransform (surface_inter.f90)'

    call fatalError(Here,"Surface: " // self % type() // " doesn't support BC." )

  end subroutine boundaryTransform

  !!
  !! Returns surface index in the shelf
  !!
  elemental function myIdx(self) result(idx)
    class(surface), intent(in) :: self
    integer(shortInt)          :: idx

    idx = self % idx

  end function myIdx

  !!
  !! Compare definition of surfaces
  !! Returns .true. if surfaces are the same (in space)
  !! Surfaces can differ in ID and name
  !!
  elemental function same_surf(LHS,RHS) result(same)
    class(surface),intent(in) :: LHS
    class(surface),intent(in) :: RHS
    logical(defBool)          :: same
    character(:),allocatable  :: LHS_def
    character(:),allocatable  :: RHS_def
    logical(defBool)          :: notHashed, sameHashes

    notHashed = (LHS % hash == UNHASHED) .or. (RHS % hash == UNHASHED)
    sameHashes = LHS % hash == RHS % hash

    if (sameHashes .or. notHashed) then
      ! Obtain definition strings
      call LHS % getDef(LHS_def)
      call RHS % getDef(RHS_def)

      ! Compare
      same = LHS_def == RHS_def

    else
      same = .false.

    end if

  end function same_surf

  !!
  !! Hashes surface definitions and saves hash for surface comparisons
  !!
  subroutine hashSurfDef(self)
    class(surface), intent(inout) :: self
    character(:),allocatable      :: defString

    ! Obtain definition string
    call self % getDef(defString)

    ! Hash and store hashed definition string
    call FNV_1(defString, self % hash)

  end subroutine hashSurfDef

!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!!
!! surfaceSlot procedures
!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!!

  !!
  !! Load allocatable surface on RHS into LHS slot
  !! Be carefull about loading slots into slots.
  !! It will work but the chain will hurt performance
  !!
  subroutine load_slot(LHS,RHS)
    class(surfaceSlot), intent(inout)         :: LHS
    class(surface),allocatable, intent(inout) :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)
    call move_alloc(RHS , LHS % slot)

    LHS % hash   = LHS % slot % hash
    LHS % surfID = LHS % slot % id()
    LHS % idx    = LHS % slot % idx

  end subroutine load_slot

  !!
  !! Copy RHS into LHS
  !! Be carefull about loading slots into slots.
  !! It will work but the chain will hurt performance
  !!
  elemental subroutine copy_slot(LHS,RHS)
    class(surfaceSlot),intent(inout) :: LHS
    class(surface), intent(in)       :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)
    allocate(LHS % slot, source = RHS)

    LHS % hash   = LHS % slot % hash
    LHS % surfID = LHS % slot % id()
    LHS % idx    = LHS % slot % idx

  end subroutine copy_slot

  !!
  !! Move allocation from RHS slot into LHS slot
  !!
  elemental subroutine moveAllocFrom_slot(LHS,RHS)
    class(surfaceSlot), intent(inout) :: LHS
    class(surfaceSlot), intent(inout) :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)
    call move_alloc(RHS % slot, LHS % slot)

    LHS % hash   = LHS % slot % hash
    LHS % surfID = LHS % slot % id()
    LHS % idx    = LHS % slot % idx

  end subroutine moveAllocFrom_slot

  !!
  !! Get id from the slot
  !!
  elemental function id_slot(self) result(id)
    class(surfaceSlot), intent(in) :: self
    integer(shortInt)              :: id

    id = self % slot % id()
    !id = self % surfId

  end function id_slot

  !!
  !! Set ID inside the slot
  !!
  elemental subroutine setId_slot(self,Id)
    class(surfaceSlot), intent(inout) :: self
    integer(shortInt), intent(in)     :: Id

    call self % slot % setId(Id)
    self % surfId = Id

  end subroutine setId_slot

  !!
  !! Get type from the slot
  !!
  elemental function type_slot(self) result(type)
    class(surfaceSlot), intent(in) :: self
    character(nameLen)             :: type

    type = self % slot % type()

  end function type_slot

  !!
  !! Get surface definition string from the slot
  !!
  pure subroutine getDef_slot(self, string)
    class(surfaceSlot), intent(in)          :: self
    character(:),allocatable, intent(inout) :: string

    call self % slot % getDef(string)

  end subroutine getDef_slot

  !!
  !! Determine if surface in slot can be boundary
  !!
  function cannotBeBoundary_slot(self) result(itCant)
    class(surfaceSlot), intent(in) :: self
    logical(defBool)               :: itCant

    itCant = self % slot % cannotBeBoundary()

  end function cannotBeBoundary_slot

  !!
  !! Set boundary conditions on the surface inside the slot
  !!
  subroutine setBoundaryConditions_slot(self,BC)
    class(surfaceSlot), intent(inout)          :: self
    integer(shortInt), dimension(:),intent(in) :: BC

    call self % slot % setBoundaryConditions(BC)

  end subroutine setBoundaryConditions_slot

  !!
  !! Hash surface definition inside the slot
  !!
  subroutine hashSurfDef_slot(self)
    class(surfaceSlot), intent(inout) :: self

    ! Hash surface inside the slot
    call self % slot % hashSurfDef()

    ! Load new hash into slot itself
    self % hash = self % slot % hash

  end subroutine hashSurfDef_slot

  !!
  !! Evaluate surface expression inside slot
  !!
  elemental subroutine evaluate_slot(self, res, r)
    class(surfaceSlot), intent(in) :: self
    real(defReal), intent(out)     :: res
    type(vector), intent(in)       :: r

    call self % slot % evaluate(res, r)

  end subroutine evaluate_slot

  !!
  !! Reflect from surface inside the slot
  !!
  elemental subroutine bounce_slot(self, r, u)
    class(surfaceSlot), intent(in)     :: self
    type(vector), intent(inout)        :: r
    type(vector), intent(inout)        :: u

    call self % slot % bounce(r, u)

  end subroutine bounce_slot

  !!
  !! Evaluate distance from the surface in the slot
  !!
  elemental subroutine distance_slot(self, dist, idx, r , u)
    class(surfaceSlot), intent(in) :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt),intent(out)  :: idx
    type(vector), intent(in)       :: r
    type(vector), intent(in)       :: u

    call self % slot % distance(dist, idx, r, u)

  end subroutine distance_slot

  !!
  !! Normal vector into +ve halfspace of the surface inside the slot
  !!
  elemental function normalVector_slot(self, r) result(normal)
    class(surfaceSlot), intent(in) :: self
    type(vector), intent(in)       :: r
    type(vector)                   :: normal

    normal = self % slot % normalVector(r)

  end function normalVector_slot

  !!
  !! Boundary transform by the surface inside the slot
  !!
  subroutine boundaryTransform_slot(self, r, u)
    class(surfaceSlot), intent(in)  :: self
    type(vector), intent(inout)     :: r
    type(vector), intent(inout)     :: u

    call self % slot % boundaryTransform(r, u)

  end subroutine boundaryTransform_slot

!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!!
!! surfaceShelf procedures
!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!!

  !!
  !! Initialise shelf
  !! Give initial size and optionaly stride to override default
  !!
  subroutine init_shelf(self,N,stride)
    class(surfaceShelf), intent(inout)     :: self
    integer(shortInt), intent(in)          :: N
    integer(shortInt), optional,intent(in) :: stride

    ! Return to uninitialised state
    call self % kill()

    ! Allocate storage space
    allocate( self % shelf(N))

    ! Assign constants
    self % Nmax = N
    if(present(stride)) self % stride = stride

  end subroutine init_shelf

  !!
  !! Adds allocatable surface to the shelf
  !! Surf needs to be allocated upon entry
  !! Deallocates surface in the process
  !! Throws error if surface with the same definition or id is already present
  !! Returns surfIdx in the surfaceShelf
  !!
  subroutine addUnique_shelf(self,surf,surfIdx)
    class(surfaceShelf), intent(inout)       :: self
    class(surface),allocatable,intent(inout) :: surf
    integer(shortInt), intent(out)           :: surfIdx
    integer(shortInt)                        :: N
    character(100),parameter :: Here = 'addUnique_shelf (surface_inter.f90)'

    if(.not.allocated(surf)) call fatalError(Here, 'surf argument is not allocated')

    ! Load surface counter for convenience
    N = self % N

    if( N == 0) then ! First surface is  being loaded
      surf % idx = 1
      call self % shelf(1) % load(surf)
      surfIdx = 1

    else
      ! There are already surfaces present (need to check uniqueness)
      associate ( storedSurf => self % shelf(1:N) )

        if ( any( storedSurf % id() == surf % id() )) then        ! ID is already present
          call fatalError(Here,'Id is already present in the surface shelf')

        else if ( any( storedSurf .same. surf )) then            ! Definition is not unique
          call fatalError(Here,'Surface with the same definition is already present')

        end if
      end associate

      ! Load surface
      surf % idx = N + 1
      call self % shelf(N+1) % load(surf)
      surfIdx = N +1

    end if

    ! Increment surface counter
    self % N = N + 1

    ! Grow storage if needed
    call self % grow()

  end subroutine addUnique_shelf

  !!
  !! Returns surfId and surfIdx of the surface with the same definition
  !! Surf needs to be allocated upon entry
  !! If there is no surface with the same definition adds it the shelf
  !! It completly ignores id of surf dummy argument
  !! surf if deallocated upon exit
  !!
  subroutine getOrAdd_shelf(self,surf,surfId,surfIdx)
    class(surfaceShelf), intent(inout)        :: self
    class(surface),allocatable, intent(inout) :: surf
    integer(shortInt), intent(out)            :: surfId
    integer(shortInt), intent(out)            :: surfIdx
    integer(shortInt)                         :: i
    character(100),parameter :: Here = 'getOrAdd_shelf (surface_inter.f90)'

    if(.not.allocated(surf)) call fatalError(Here, 'surf argument is not allocated')

    ! Find surface with the same definition
    do i=1,self % N
      if(surf .same. self % shelf(i)) then ! Return Idx & ID of existing surface
        surfId  = self % shelf(i) % id()
        surfIdx = i
        deallocate(surf)
        return

      end if
    end do

    ! Surface was not found set a next free Id and add to shelf
    surfId = self % freeId()
    call surf % setId(surfId)
    call self % addUnique(surf,surfIdx)

  end subroutine getOrAdd_shelf

  !!
  !! Returns smallest +ve free id in the shelf
  !!
  function freeId_shelf(self) result(id)
    class(surfaceShelf), intent(inout) :: self
    integer(shortInt)                  :: id

    associate( storedSurf => self % shelf(1:self % N) )

      do
        id = self % nextId

        ! Check if the proposed id is unique
        if ( all(storedSurf % id() /= id )) return

        ! Increment test id if current nextId is occupied
        self % nextId = self % nextId + 1

      end do
    end associate

  end function freeId_shelf

  !!
  !! Return surface idx given surf ID
  !! Returns "targetNotFound" if id is not present
  !!
  elemental function getIdx_shelf(self,id) result(idx)
    class(surfaceShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: id
    integer(shortInt)               :: idx
    integer(shortInt)               :: i

    associate( storedSurf => self % shelf(1:self % N))
      idx = linFind(storedSurf % id(), id)

    end associate
  end function getIdx_shelf

  !!
  !! Increases the size of the shelf by stride if N >= Nmax
  !!
  subroutine grow_shelf(self)
    class(surfaceShelf), intent(inout)         :: self
    integer(shortInt)                          :: Nmax

    ! Return if shelf does not need to grow
    if ( self % N < self % Nmax) then
      return
    end if

    ! Calculate new maximum size
    Nmax = self % Nmax + self % stride
    call self % resize(Nmax)


  end subroutine grow_shelf

  !!
  !! Set shelf size to the number of stored entries
  !!
  subroutine trimSize_shelf(self)
    class(surfaceShelf), intent(inout) :: self

    call self % resize(self % N)

  end subroutine trimSize_shelf

  !!
  !! Set size of the shelf to Nmax >= self % N
  !!
  subroutine resize_shelf(self,Nmax)
    class(surfaceShelf), intent(inout)         :: self
    integer(shortInt), intent(in)              :: Nmax
    type(surfaceSlot),dimension(:),allocatable :: tempShelf
    character(100), parameter :: Here ='resize (surface_inter.f90)'

    if( Nmax < self % N) call fatalError(Here,'Resizeing to size smaller then number of entries')

    ! Allocate new storage space
    allocate(tempShelf(Nmax) )

    ! Copy Old enteries without reallocation
    associate ( new => tempShelf(1:self % N), &
                old => self % shelf(1:self % N) )

      call new % moveAllocFrom( old)
    end associate

    ! Replace shelf and update Nmax
    call move_alloc(tempShelf, self % shelf)
    self % Nmax = Nmax

  end subroutine resize_shelf

  !!
  !! Deallocate storage space and returns shelf to uninitialised state
  !!
  subroutine kill_shelf(self)
    class(surfaceShelf), intent(inout) :: self

    if(allocated(self % shelf)) deallocate(self % shelf)
    self % Nmax   = 0
    self % N      = 0
    self % stride = DEF_STRIDE
    self % nextID = 1

  end subroutine kill_shelf


!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!!
!! Not type-bound procedures
!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!!

  !!
  !! Write parameter using local format specification defined for surfaces in this module
  !!
  pure function real2Char_surf(num) result(str)
    real(defReal), intent(in) :: num
    character(int_width)      :: str

    write(str,'('//realForm//')') num

  end function real2Char_surf

  !!
  !! Template function for writing definition strings
  !!
  pure function printSurfDef(type,coeff) result (str)
    character(nameLen), intent(in)         :: type
    real(defReal),dimension(:), intent(in) :: coeff
    character(:), allocatable              :: str
    character(size(coeff)*int_width)       :: tail
    integer(shortInt)                      :: i

    tail = ''
    do i=1,size(coeff)
      tail = trim(tail) // real2Char_surf(coeff(i))
    end do

    str = 'Type: '// trim(type)// '; Coeff:' // tail

  end function printSurfDef

end module surface_inter
