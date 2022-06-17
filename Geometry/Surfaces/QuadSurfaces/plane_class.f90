module plane_class

  use numPrecision
  use universalVariables, only : X_AXIS, Y_AXIS, Z_AXIS, INF
  use genericProcedures,  only : fatalError, dotProduct, numToChar
  use dictionary_class,   only : dictionary
  use quadSurface_inter,  only : quadSurface
  use surface_inter,      only : kill_super => kill
  implicit none
  private

  !!
  !! (Axis) Plane
  !!
  !! Planes with a normal which is one of the principal axis (x,y or z)
  !!
  !! F(r) = c1 * x + c2 * y + c3 * z - c4
  !!
  !! Surface tolerance: SURF_TOL
  !!
  !! Sample dictionary input:
  !!  pl { type plane; id 16; coeffs (1.0 2.0 1.0 0.0);}
  !!
  !! Private members:
  !!   norm -> Normal vector (normalised) [c1, c2, c3]
  !!   offset -> Offset, c4 (also normalised)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: plane
    private
    real(defReal), dimension(3) :: norm = ZERO
    real(defReal)               :: offset = ZERO
  contains
    ! Superclass procedures
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill
    
    ! Local procedure
    procedure :: build
  end type plane


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(plane), intent(in)  :: self
    character(:), allocatable  :: str

    str = 'plane'

  end function myType

  !!
  !! Initialise plane from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if id < 0.
  !!   fatalError if is not a plane (coeffcients 1-3 are 0.0)
  !!
  subroutine init(self, dict)
    class(plane), intent(inout)              :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: coeffs
    character(100), parameter :: Here = 'init (plane_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(coeffs,'coeffs')

    ! Check values
    if (id < 1) then
      call fatalError(Here,'Invalid surface id provided. ID must be > 1')

    else if (size(coeffs) /= 4) then
      call fatalError(Here, '4 plane coefficients must be given. There are: '//&
                            numToChar(size(coeffs)))
    else if (all(coeffs(1:3) == ZERO)) then
      call fatalError(Here, 'Invalid plane normal. coefficients 1-3 are all 0.0.')

    end if

    ! Normalise coefficients
    coeffs = coeffs / norm2(coeffs(1:3))

    ! Load data
    call self % setID(id)
    self % norm = coeffs(1:3)
    self % offset = coeffs(4)

  end subroutine init

  !!
  !! Build plane from components
  !!
  !! Args:
  !!   id [in] -> Surface ID
  !!   norm [in]   -> normal vector of plane (normalised if not already)
  !!   offset [in] -> offset of plane
  !!
  !! Errors:
  !!   fatalError if id or radius are -ve
  !!
  subroutine build(self, id, norm, offset)
    class(plane), intent(inout)              :: self
    integer(shortInt), intent(in)            :: id
    real(defReal), dimension(:), intent(in)  :: norm 
    real(defReal), intent(in)                :: offset
    character(100), parameter :: Here = 'build (plane_class.f90)'

    if (id < 1) call fatalError(Here,'Invalid surface id provided. ID must be > 1')

    call self % setID(id)

    self % norm = norm / norm2(norm)
    self % offset = offset

  end subroutine build

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  !! Always returns infinate box (even when aligned with some axis)
  !!
  pure function boundingBox(self) result(aabb)
    class(plane), intent(in)   :: self
    real(defReal), dimension(6) :: aabb

    aabb(1:3) = -INF
    aabb(4:6) = INF


  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(plane), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    c = dotProduct(r, self % norm) - self % offset

  end function evaluate


  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Solve linear equaltion
  !!   d*k + c = 0.0
  !!
  !!   k = F'(r) .dot. u = c1*u1 + c2*u2 +c3*u3
  !!   c = F(r)
  !!
  pure function distance(self, r, u) result(d)
    class(plane), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    real(defReal)                           :: k, c

    k = dotProduct(u, self % norm)
    c = self % evaluate(r)

    if ( k == ZERO .or. abs(c) < self % surfTol()) then ! Parallel or at the surface
      d = INF

    else
      d = -c/k
      if (d <= ZERO .or. d > INF) d = INF

    end if

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  !! Note:
  !!   For parallel direction halfspace is asigned by the sign of `evaluate` result.
  !!
  pure function going(self, r, u) result(halfspace)
    class(plane), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal)                           :: proj

    proj = dotProduct(u, self % norm)
    halfspace = proj > ZERO

    ! Special case of parallel direction
    ! Partilce stays in its current halfspace
    if (proj == ZERO) then
      halfspace = self % evaluate(r) >= ZERO
    end if

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(plane), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % norm = ZERO
    self % offset = ZERO

  end subroutine kill

end module plane_class
