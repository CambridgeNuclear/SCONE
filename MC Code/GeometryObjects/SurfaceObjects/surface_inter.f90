!
! The surface class is a base class for second order surfaces
! These surfaces have the form:
! Axx + Byy + Czz + Dxy + Eyz + Fxz + Gx + Hy + Jz + K
!
! Each surface class should contain the data necessary to define it
! and the ability to locate which halfspace a given point occupies
!
module surface_inter
  use numPrecision
  use universalVariables
  use genericProcedures, only : dotProduct

  implicit none
  private
    
  type, abstract, public :: surface
    logical(defBool)            :: isReflective = .FALSE.
    logical(defBool)            :: isPeriodic   = .FALSE.
    logical(defBool)            :: isVacuum     = .FALSE.
    logical(defBool)            :: isCompound   = .FALSE.
    real(defReal), dimension(3) :: periodicTranslation
    character(100)              :: name =""
    integer(shortInt)           :: id = 0
  contains
    procedure                                    :: halfspace
    procedure                                    :: reflect
    procedure(evaluate), deferred                :: evaluate
    procedure(reflectiveTransform), deferred     :: reflectiveTransform
    procedure(distanceToSurface), deferred       :: distanceToSurface
    procedure(normalVector), deferred            :: normalVector
    procedure(whichSurface), deferred            :: whichSurface
    procedure(setBoundaryConditions), deferred   :: setBoundaryConditions
    procedure(boundaryTransform), deferred       :: boundaryTransform
  end type surface

  type, public :: surface_ptr
    class(surface), pointer :: ptr => null()
  contains
    procedure :: halfspace => halfspace_ptr
    procedure :: reflect => reflect_ptr
    procedure :: evaluate => evaluate_ptr
    procedure :: reflectiveTransform => reflectiveTransform_ptr
    procedure :: distanceToSurface => distanceToSurface_ptr
    procedure :: normalVector => normalVector_ptr
    procedure :: whichSurface => whichSurface_ptr
    procedure :: isReflective => isReflective_ptr
    procedure :: isVacuum => isVacuum_ptr
    procedure :: isPeriodic => isPeriodic_ptr
    procedure :: periodicTranslation => periodicTranslation_ptr
    procedure :: setBoundaryConditions => setBoundaryConditions_ptr
    procedure :: boundaryTransform => boundaryTransform_ptr
    procedure :: name => name_ptr
    procedure :: id
    procedure :: kill
    generic   :: assignment(=) => surface_ptr_assignment, surface_ptr_assignment_target!,surface_ptr_assignment_pointer

    procedure,private :: surface_ptr_assignment
    procedure,private :: surface_ptr_assignment_target
    !procedure,private :: surface_ptr_assignment_pointer
  end type surface_ptr

  abstract interface

    function evaluate(self, r) result(res)
      use numPrecision
      import :: surface
      implicit none
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal)                           :: res
    end function

    subroutine reflectiveTransform(self, r, u)
      use numPrecision
      import :: surface
      implicit none
      class(surface), intent(in)                 :: self
      real(defReal), dimension(3), intent(inout) :: r, u
    end subroutine

    function distanceToSurface(self, r, u) result(distance)
      use numPrecision
      import :: surface
      implicit none
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r, u
      real(defReal)                           :: distance
    end function

    function normalVector(self, r) result(normal)
      use numPrecision
      import :: surface
      implicit none
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3)             :: normal
    end function

    function whichSurface(self, r, u) result(surfPointer)
      use numPrecision
      use genericProcedures
      import :: surface
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r, u
      class(surface), pointer                 :: surfPointer
    end function

    subroutine setBoundaryConditions(self, BC)
      use numPrecision
      use genericProcedures
      import :: surface
      class(surface), intent(inout)               :: self
      integer(shortInt), dimension(6), intent(in) :: BC
    end subroutine setBoundaryConditions

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

  subroutine reflectiveTransform_ptr(self, r, u)
    class(surface_ptr), intent(in)             :: self
    real(defReal), dimension(3), intent(inout) :: r, u
    call self%ptr%reflectiveTransform(r,u)
  end subroutine reflectiveTransform_ptr

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

  function periodicTranslation_ptr(self) result(periodicTranslation)
    class(surface_ptr), intent(in) :: self
    real(defReal), dimension(3)    :: periodicTranslation
    periodicTranslation = self % ptr % periodicTranslation
  end function periodicTranslation_ptr

  subroutine setBoundaryConditions_ptr(self,BC)
    class(surface_ptr), intent(in)              :: self
    integer(shortInt), dimension(6), intent(in) :: BC
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

  subroutine kill(self)
    class(surface_ptr), intent(inout) :: self
    self % ptr => null()
  end subroutine kill

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

end module surface_inter
