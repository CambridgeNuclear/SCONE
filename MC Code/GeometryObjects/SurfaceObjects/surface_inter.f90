module surface_inter
  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, dotProduct
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

  ! Public functions
  public :: printSurfDef  !! Print surface definition string


  !!
  !! Abstract interface for a surface
  !!
  type, abstract, public :: surface
    private
    !! Move the entire public interface into procedures ENCAPSULATION GOD DAMN IT!
!    character(nameLen)          :: name =""
!    integer(shortInt)           :: id = 0
!    integer(shortInt)           :: hash = UNHASHED

    ! Perhaps to be removed -> will be moved to a function (better slot support etc.
    logical(defBool)            :: isCompound   = .FALSE.

  contains
    ! Operators and assignments
    generic                                      :: operator(.same.) => same_surf

    ! Initialisation & Indentification procedures
    procedure(name),deferred                     :: name
    procedure(id),deferred                       :: id
    procedure(hash),deferred                     :: hash
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
!
!  !!
!  !!
!  !!
!  type,public,extends(surface) :: surfaceSlot
!    private
!    class(surface),allocatable :: slot
!  contains
!    ! Initialisation & Identification procedures
!    procedure :: name                  => name_slot
!    procedure :: id                    => id_slot
!    procedure :: hash                  => hash_id
!    procedure :: type                  => type_slot
!    procedure :: getDef                => getDef_slot
!    procedure :: cannotBeBoundary      => cannotBeBoundary_slot
!    procedure :: setBoundaryConditions => setBoundaryConditions_slot
!    procedure :: hashSurfDef           => hashSurfDef_slot
!
!    ! Run time procedures
!    procedure :: evaluate          => evaluate_slot
!    procedure :: reflect           => reflect_slot
!    procedure :: distance          => distance_slot
!    procedure :: normalVector      => normalVector_slot
!    procedure :: boundaryTransform => boundaryTransform_slot
!
!  end type surfaceSlot

  abstract interface
    !!
    !! Returns surface name
    !!
    elemental function name(self)
      import :: nameLen,&
                surface
      class(surface), intent(in) :: self
      character(nameLen)         :: name
    end function name

    !!
    !! Returns surface ID
    !!
    elemental function id(self)
      import :: shortInt,&
                surface
      class(surface), intent(in) :: self
      integer(shortInt)          :: id
    end function id

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
    elemental function evaluate(self, r) result(res)
      import :: surface, &
                defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal)                           :: res
    end function evaluate

    !!
    !! OBSOLETE
    !!
    pure subroutine reflectiveTransform(self, r, u)
      import :: surface, &
                defReal
      class(surface), intent(in)                 :: self
      real(defReal), dimension(3), intent(inout) :: r, u
    end subroutine reflectiveTransform

    !!
    !! Return +ve distance to surface from point r along direction u
    !! Return INFINITY if there is no crossing
    !!
    elemental function distance(self, r, u) result(dist)
      import :: surface, &
                defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r, u
      real(defReal)                           :: dist
    end function distance

    !!
    !! Return vector normal to the surface for a point r on the surface
    !! Vector is pointing into +ve halfspace
    !! No check if r lies on the surface is performed
    !!
    elemental function normalVector(self, r) result(normal)
      import :: surface, &
                defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3)             :: normal
    end function normalVector

    !!
    !! Perform reflection of point r and direction u by the surface
    !! Perform coordinate transformation of a point r with direction u by the surface
    subroutine boundaryTransform(self, r, u, isVacuum)
      use numPrecision
      use genericProcedures
      import :: surface
      class(surface), intent(in)                 :: self
      real(defReal), intent(inout), dimension(3) :: r
      real(defReal), intent(inout), dimension(3) :: u
      logical(defBool), intent(inout)            :: isVacuum
    end subroutine boundaryTransform

  end interface

contains

!!
!! Base surface class procedures
!!
  !!
  !! Determine whether a point occupies the positive or negative halfspace of a surface
  !! Point can also be located on a surface - must include direction to determine halfspace
  !!
  function halfspace(self,r,u) result(position)
    class(surface), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r, &  ! position relative to the surface
                                               u     ! direction of travel (for coincidence cases)
    real(defReal)                           :: res
    logical(defBool)                        :: position

    res = self % evaluate(r)

    ! Point is close to the surface - check direction to determine whether it will be in the
    ! positive or negative halfspace
    if(abs(res) < surface_tol) then
      position = (dotProduct(u, self % normalVector(r)) > ZERO)
      return
    else if (res > ZERO) then
      position = infront
      return
    else
      position = behind
      return
    end if

  end function halfspace

  !!
  !! Reflect a particle incident on a surface and nudge it away from the surface
  !!
  subroutine reflect(self, r, u)
    class(surface), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r, &
                                                  u
    real(defReal), dimension(3)                :: normal
    real(defReal)                              :: magSquared

    normal = self%normalVector(r)
    magSquared = dotProduct(normal,normal)

    u = u - TWO*dotProduct(u,normal)*normal/magSquared
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

    notHashed = (LHS % hash() == UNHASHED) .or. (RHS % hash() == UNHASHED)
    sameHashes = LHS % hash() == RHS % hash()

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
