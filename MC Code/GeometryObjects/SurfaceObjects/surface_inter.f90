module surface_inter
  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, dotProduct
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
  character(*),parameter         :: defTol   = '6'
  character(*),parameter         :: expSize  = '2'
  character(*),parameter         :: width    = '13'  ! defTol + expSize + 5
  integer(shortInt), parameter   :: int_width = 13   ! Needs to be the same as width but as int
  integer(shortInt),parameter    :: UNHASHED = 0
  character(*),parameter,public  :: realForm = 'ES'//width//'.'//defTol//'E'//expSize !
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
    procedure                                    :: halfspace
    procedure                                    :: reflect
    procedure(evaluate), deferred                :: evaluate
    procedure(distance), deferred                :: distance
    procedure(normalVector), deferred            :: normalVector
    procedure(boundaryTransform), deferred       :: boundaryTransform

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
    procedure :: reflect           => reflect_slot
    procedure :: distance          => distance_slot
    procedure :: normalVector      => normalVector_slot
    procedure :: boundaryTransform => boundaryTransform_slot

  end type surfaceSlot

  !!
  !!
  !!
  type,public :: surfaceShelf
    private
    type(surfaceSlot),dimension(:),allocatable, public :: shelf
    integer(shortInt)                                  :: Nmax   = 0
    integer(shortInt)                                  :: N      = 0
    integer(shortInt)                                  :: stride = DEF_STRIDE
    integer(shortInt)                                  :: nextID = 0

  contains
    procedure :: init       => init_shelf
    procedure :: kill       => kill_shelf
    !procedure :: addUnique  => addUnique_shelf
    !procedure :: getOrAdd   => getOrAdd_shelf
    !procedure :: freeID     => freeId_shelf
    !procedure :: getIdx     => getIdx_shelf
    !procedure :: resize     => resize_shelf

    !procedure, private :: add_shelf
    !procedure, private :: grow_shelf
    !procedure, private :: containsID_shelf
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
    elemental subroutine evaluate(self,res, r, shelf)
      import :: surface, &
                defReal, &
                vector, &
                surfaceShelf
      class(surface), intent(in)              :: self
      real(defReal), intent(out)              :: res
      type(vector), intent(in)                :: r
      type(surfaceShelf), intent(in)          :: shelf
    end subroutine evaluate

    !!
    !! Return +ve distance to surface from point r along direction u
    !! Return INFINITY if there is no crossing
    !! Also return index of the surface being X-ed
    !!
    elemental subroutine distance(self, dist, idx, r, u, shelf)
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
      type(surfaceShelf), intent(in)          :: shelf
    end subroutine distance

    !!
    !! Return vector normal to the surface for a point r on the surface
    !! Vector is pointing into +ve halfspace
    !! No check if r lies on the surface is performed
    !!
    elemental function normalVector(self, r, shelf) result(normal)
      import :: surface, &
                vector, &
                surfaceShelf
      class(surface), intent(in)      :: self
      type(vector), intent(in)        :: r
      type(surfaceShelf), intent(in)  :: shelf
      type(vector)                    :: normal
    end function normalVector

    !!
    !! Perform coordinate transformation of a point r with direction u by the surface
    !! Set isVacuum to true if vacuum BC was X-ed
    !!
    subroutine boundaryTransform(self, r, u, isVacuum, shelf)
      import :: surface, &
                vector, &
                defBool, &
                surfaceShelf
      class(surface), intent(in)                 :: self
      type(vector), intent(inout)                :: r
      type(vector), intent(inout)                :: u
      logical(defBool), intent(out)              :: isVacuum
      type(surfaceShelf), intent(in)             :: shelf
    end subroutine boundaryTransform

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
  !!
  elemental function halfspace(self, r, u, shelf) result(position)
    class(surface), intent(in)              :: self
    type(vector),intent(in)                 :: r
    type(vector),intent(in)                 :: u
    type(surfaceShelf), intent(in)          :: shelf
    real(defReal)                           :: res
    logical(defBool)                        :: position

    call self % evaluate(res, r, shelf)

    ! Set halfspace based on sign of res
    if( res > ZERO) then
      position = infront
    else
      position = behind
    end if

    ! If res is within surface tolerance choose halfspace based on direction
    if(abs(res) < surface_tol) then
      position = (self % normalVector(r,shelf) .dot. u) > ZERO
    end if

  end function halfspace

  !!
  !! Reflect a particle incident on a surface and nudge it away from the surface
  !!
  elemental subroutine reflect(self, r, u, shelf)
    class(surface), intent(in)                 :: self
    type(vector),intent(inout)                 :: r
    type(vector),intent(inout)                 :: u
    type(surfaceShelf), intent(in)             :: shelf
    type(vector)                               :: normal
    real(defReal)                              :: magSquared

    normal     = self % normalVector(r, shelf)
    magSquared = normal .dot. normal

    u = u - TWO * (u .dot. normal) * normal / magSquared
    r = r + NUDGE * u

  end subroutine reflect

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
  elemental subroutine evaluate_slot(self, res, r, shelf)
    class(surfaceSlot), intent(in) :: self
    real(defReal), intent(out)     :: res
    type(vector), intent(in)       :: r
    type(surfaceShelf), intent(in) :: shelf

    call self % slot % evaluate(res, r, shelf)

  end subroutine evaluate_slot

  !!
  !! Reflect from surface inside the slot
  !!
  elemental subroutine reflect_slot(self, r, u, shelf)
    class(surfaceSlot), intent(in)     :: self
    type(vector), intent(inout)        :: r
    type(vector), intent(inout)        :: u
    type(surfaceShelf), intent(in)     :: shelf

    call self % slot % reflect(r, u, shelf)

  end subroutine reflect_slot

  !!
  !! Evaluate distance from the surface in the slot
  !!
  elemental subroutine distance_slot(self, dist, idx, r , u, shelf)
    class(surfaceSlot), intent(in) :: self
    real(defReal), intent(out)     :: dist
    integer(shortInt),intent(out)  :: idx
    type(vector), intent(in)       :: r
    type(vector), intent(in)       :: u
    type(surfaceShelf), intent(in) :: shelf

    call self % slot % distance(dist, idx, r, u, shelf)

  end subroutine distance_slot

  !!
  !! Normal vector into +ve halfspace of the surface inside the slot
  !!
  elemental function normalVector_slot(self, r, shelf) result(normal)
    class(surfaceSlot), intent(in) :: self
    type(vector), intent(in)       :: r
    type(vector)                   :: normal
    type(surfaceShelf),intent(in)  :: shelf

    normal = self % slot % normalVector(r, shelf)

  end function normalVector_slot

  !!
  !! Boundary transform by the surface inside the slot
  !!
  subroutine boundaryTransform_slot(self, r, u, isVacuum, shelf)
    class(surfaceSlot), intent(in)  :: self
    type(vector), intent(inout)     :: r
    type(vector), intent(inout)     :: u
    logical(defBool), intent(out)   :: isVacuum
    type(surfaceShelf), intent(in)  :: shelf

    call self % slot % boundaryTransform(r, u, isVacuum, shelf)

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
  !! Deallocate storage space and returns shelf to uninitialised state
  !!
  subroutine kill_shelf(self)
    class(surfaceShelf), intent(inout) :: self

    if(allocated(self % shelf)) deallocate(self % shelf)
    self % Nmax   = 0
    self % N      = 0
    self % stride = DEF_STRIDE
    self % nextID = 0
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
