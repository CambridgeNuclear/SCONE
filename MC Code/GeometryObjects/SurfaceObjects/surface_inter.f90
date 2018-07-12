module surface_inter
  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, dotProduct
  use hashFunctions_func, only : FNV_1

  implicit none
  private

  ! Local Paramethers
  character(*),parameter         :: defTol   = '6'
  character(*),parameter         :: expSize  = '2'
  character(*),parameter         :: width    = '13'  ! defTol + expSize + 5
  integer(shortInt), parameter   :: int_width = 13   ! Needs to be the same as width but as int
  character(*),parameter,public  :: realForm = 'ES'//width//'.'//defTol//'E'//expSize
  integer(shortInt),parameter    :: UNHASHED = 0

  !! Public function to print surface definition string
  public :: printSurfDef


  !!
  !! Abstract interface for a surface
  !!
  type, abstract, public :: surface
    private
    !! Move the entire public interface into procedures ENCAPSULATION GOD DAMN IT!
    character(nameLen)          :: name =""
    integer(shortInt)           :: id = 0
    integer(shortInt)           :: hash = UNHASHED

    ! Perhaps to be removed -> will be moved to a function (better slot support etc.
    logical(defBool)            :: isCompound   = .FALSE.

  contains
    ! Interface to stay
    procedure                                    :: halfspace
    procedure                                    :: reflect
    procedure(evaluate), deferred                :: evaluate
    procedure(distanceToSurface), deferred       :: distanceToSurface
    procedure(normalVector), deferred            :: normalVector
    procedure(boundaryTransform), deferred       :: boundaryTransform

    ! New interface
    generic                                      :: operator(.same.) => same_surf
    procedure                                    :: cannotBeBoundary
    procedure                                    :: setBoundaryConditions
    procedure                                    :: hashSurfDef
    procedure, private                           :: same_surf
    procedure(type),deferred                     :: type
    procedure(getDef),deferred                   :: getDef

    ! Old components
    procedure(name),deferred                     :: name
    procedure(hash),deferred                     :: hash
    procedure(id),deferred                       :: id


    ! To delate
    procedure(whichSurface), deferred            :: whichSurface


  end type surface

!  type, public :: surface_ptr
!    class(surface), pointer :: ptr => null()
!  contains
!    procedure :: halfspace => halfspace_ptr
!    procedure :: reflect => reflect_ptr
!    procedure :: evaluate => evaluate_ptr
!    procedure :: distanceToSurface => distanceToSurface_ptr
!    procedure :: normalVector => normalVector_ptr
!    procedure :: whichSurface => whichSurface_ptr
!    procedure :: isReflective => isReflective_ptr
!    procedure :: isVacuum => isVacuum_ptr
!    procedure :: isPeriodic => isPeriodic_ptr
!    procedure :: periodicTranslation => periodicTranslation_ptr
!    procedure :: setBoundaryConditions => setBoundaryConditions_ptr
!    procedure :: boundaryTransform => boundaryTransform_ptr
!    procedure :: name => name_ptr
!    procedure :: id
!    generic   :: assignment(=) => surface_ptr_assignment, surface_ptr_assignment_target!,surface_ptr_assignment_pointer
!
!    ! New interface
!    procedure :: type_ptr
!    procedure :: getDef_ptr
!    procedure :: cannotBeBoundary_ptr
!
!    procedure,private :: surface_ptr_assignment
!    procedure,private :: surface_ptr_assignment_target
!  end type surface_ptr

  abstract interface

    !!
    !! Return character with name of surface type
    !!
    function type(self)
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
    function evaluate(self, r) result(res)
      import :: surface, &
                defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal)                           :: res
    end function evaluate

    !!
    !! OBSOLETE
    !!
    subroutine reflectiveTransform(self, r, u)
      import :: surface, &
               defReal
      class(surface), intent(in)                 :: self
      real(defReal), dimension(3), intent(inout) :: r, u
    end subroutine reflectiveTransform

    !!
    !! Return +ve distance to surface from point r along direction u
    !! Return INFINITY if there is no crossing
    !!
    function distanceToSurface(self, r, u) result(distance)
      import :: surface, &
                defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r, u
      real(defReal)                           :: distance
    end function

    !!
    !! Return vector normal to the surface for a point r on the surface
    !! Vector is pointing into +ve halfspace
    !! No check if r lies on the surface is performed
    !!
    function normalVector(self, r) result(normal)
      import :: surface, &
                defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3)             :: normal
    end function

    !!
    !! WILL BECOME OBSOLETE
    !!
    function whichSurface(self, r, u) result(surfPointer)
      use numPrecision
      use genericProcedures
      import :: surface
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r, u
      class(surface), pointer                 :: surfPointer
    end function

    !!
    !! Perform reflection of point r and direction u by the surface
    !!
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

!!
!! Surface pointer procedures
!!

  function evaluate_ptr(self, r) result(res)
    class(surface_ptr), intent(in)          :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: res
    res=self%ptr%evaluate(r)
  end function evaluate_ptr

  subroutine reflect_ptr(self, r, u)
   class(surface_ptr), intent(in)             :: self
   real(defReal), dimension(3), intent(inout) :: r, u
   call self%ptr%reflect(r,u)
  end subroutine reflect_ptr

  function distanceToSurface_ptr(self, r, u) result(distance)
    class(surface_ptr), intent(in)          :: self
    real(defReal), dimension(3), intent(in) :: r, u
    real(defReal)                           :: distance
    distance = self%ptr%distanceToSurface(r,u)
  end function distanceToSurface_ptr

  function normalVector_ptr(self, r) result(normal)
    class(surface_ptr), intent(in)          :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: normal
    normal = self%ptr%normalVector(r)
  end function normalVector_ptr

  function halfspace_ptr(self,r,u) result(position)
    class(surface_ptr), intent(in)          :: self
    real(defReal), dimension(3), intent(in) :: r, &  ! position relative to the surface
                                               u     ! direction of travel (for coincidence cases)
    logical(defBool)                        :: position
    position = self%ptr%halfspace(r,u)
  end function halfspace_ptr

  function whichSurface_ptr(self, r, u) result(surfPointer)
    class(surface_ptr), intent(in)          :: self
    real(defReal), dimension(3), intent(in) :: r, u
    class(surface), pointer                 :: surfPointer
    surfPointer => self%ptr%whichSurface(r, u)
  end function

  function isReflective_ptr(self) result(isReflective)
    class(surface_ptr), intent(in) :: self
    logical(defBool)               :: isReflective
    isReflective = self % ptr % isReflective
  end function isReflective_ptr

  function isVacuum_ptr(self) result(isVacuum)
    class(surface_ptr), intent(in) :: self
    logical(defBool)               :: isVacuum
    isVacuum = self % ptr % isVacuum
  end function isVacuum_ptr

  function isPeriodic_ptr(self) result(isPeriodic)
    class(surface_ptr), intent(in) :: self
    logical(defBool)               :: isPeriodic
    isPeriodic = self % ptr % isPeriodic
  end function isPeriodic_ptr

  function id(self) result(ind)
    class(surface_ptr), intent(in) :: self
    integer(shortInt)              :: ind
    ind = self % ptr % id
  end function id

  function type_ptr(self) result(type)
    class(surface_ptr), intent(in) :: self
    character(nameLen)             :: type

    type = self % ptr % type()

  end function type_ptr

  pure subroutine getDef_ptr(self,string)
    class(surface_ptr), intent(in)         :: self
    character(:),allocatable,intent(inout) :: string

    call self % ptr % getDef(string)

  end subroutine getDef_ptr

  function cannotBeBoundary_ptr(self) result(itCant)
    class(surface_ptr), intent(in) :: self
    logical(defBool)               :: itCant

    itCant = self % ptr % cannotBeBoundary()

  end function cannotBeBoundary_ptr

  function periodicTranslation_ptr(self) result(periodicTranslation)
    class(surface_ptr), intent(in) :: self
    real(defReal), dimension(3)    :: periodicTranslation
    periodicTranslation = self % ptr % periodicTranslation
  end function periodicTranslation_ptr

  subroutine setBoundaryConditions_ptr(self,BC)
    class(surface_ptr), intent(in)              :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    call self % ptr % setBoundaryConditions(BC)
  end subroutine setBoundaryConditions_ptr

  subroutine boundaryTransform_ptr(self,r,u,isVacuum)
    class(surface_ptr), intent(in)             :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u
    logical(defBool), intent(inout)            :: isVacuum
    call self % ptr % boundaryTransform(r,u,isVacuum)
  end subroutine boundaryTransform_ptr

  function name_ptr(self) result(name)
    class(surface_ptr), intent(in) :: self
    character(100)                 :: name
    name = self % ptr % name
  end function name_ptr

  subroutine surface_ptr_assignment(LHS,RHS)
    class(surface_ptr), intent(out)  :: LHS
    type(surface_ptr), intent(in)    :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS % ptr
  end subroutine surface_ptr_assignment

  subroutine surface_ptr_assignment_target(LHS,RHS)
    class(surface_ptr), intent(out)        :: LHS
    class(surface), pointer, intent(in)    :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS
  end subroutine surface_ptr_assignment_target

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
