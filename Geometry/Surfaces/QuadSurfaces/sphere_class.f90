module sphere_class

  use numPrecision
  use universalVariables, only : INF, SURF_TOL
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : kill_super => kill
  use quadSurface_inter,  only : quadSurface

  implicit none
  private

  ! Local constants
  character(*), parameter :: TYPE_NAME    = 'sphere'

  !!
  !! Sphere surface
  !!
  !! F(r) = (r1-x0)^2 + (r2-y0)^2 + (r3-z0)^2  - R^2
  !!
  !! Surface tolerance: 2 * R * SURF_TOL
  !!
  !! Sample dictionary input:
  !!   sph { type sphere;
  !!         id 17;
  !!         origin (-1.0 -2.0 0.0);
  !!         radius 5.0;
  !!       }
  !!
  !! Private Members:
  !!   origin -> Location of the midle of the sphere
  !!   r      -> Sphere radius
  !!   r_sq   -> Square of radius r (r^2)
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(quadSurface) :: sphere
    private
    real(defReal), dimension(3) :: origin = ZERO
    real(defReal)               :: r      = ZERO
    real(defReal)               :: r_sq   = ZERO
  contains
    ! Superclass procedures
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: normal
    procedure :: kill
  end type sphere

contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(sphere), intent(in) :: self
    character(:), allocatable  :: str

    str = TYPE_NAME

  end function myType

  !!
  !! Initialise sphere from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if radius or id < 0.
  !!
  subroutine init(self, dict)
    class(sphere), intent(inout)  :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: origin
    character(100), parameter :: Here = 'init (sphere_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(self % r, 'radius')
    call dict % get(origin, 'origin')

    ! Check values
    if (id < 1) then
      call fatalError(Here,'Invalid surface id provided. ID must be > 1')

    else if (size(origin) /= 3) then
      call fatalError(Here,'Origin needs to have size 3. Has: '//numToChar(size(origin)))

    else if ( self % r <= ZERO) then
      call fatalError(Here, 'Radius of sphere must be +ve. Is: '//numToChar(self % r))

    end if

    ! Load data
    self % r_sq = self % r * self % r
    self % origin = origin
    call self % setID(id)

    ! Set surface tolerance
    call self % setTol( TWO * self % r * SURF_TOL)

  end subroutine init

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(sphere), intent(in)   :: self
    real(defReal), dimension(6) :: aabb

    aabb(1:3) = self % origin - [self % r, self % r, self % r]
    aabb(4:6) = self % origin + [self % r, self % r, self % r]

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(c)
    class(sphere), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c
    real(defReal), dimension(3)             :: diff

    diff = r - self % origin

    c = dot_product(diff, diff) - self % r_sq

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Solves quadratic intersection equation
  !!   d^2 + 2kd + c = 0
  !!   c = F(r)
  !!   k = (r1-x0)u1 + (r2-y0)u2 + (r3-z0)u3
  !!
  pure function distance(self, r, u) result(d)
    class(sphere), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    real(defReal)                           :: c, k, delta

    ! Calculate quadratic components
    c = self % evaluate(r)
    k = dot_product(r - self % origin, u)
    delta = k*k - c  ! Technically delta/4

    ! Calculate the distance
    if (delta < ZERO) then ! No intersection
      d = INF

    else if (abs(c) < self % surfTol()) then ! Point at a surface
      if ( k >= ZERO) then
        d = INF
      else
        d = -k + sqrt(delta)
      end if

    else if (c < ZERO) then ! Point inside the surface
      d = -k + sqrt(delta)

    else ! Point outside the surface
      d = -k - sqrt(delta)
      if (d <= ZERO) d = INF

    end if
  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  pure function going(self, r, u) result(halfspace)
    class(sphere), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace

    halfspace = dot_product(r - self % origin, u) >= ZERO

  end function going
  
  !!
  !! Return the normal corresponding to the sphere surface
  !!
  pure function normal(self, r, u) result(n)
    class(sphere), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal), dimension(3)             :: n

    n = r - self % origin
    n = n / norm2(n)

  end function normal


  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(sphere), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % origin = ZERO
    self % r      = ZERO
    self % r_sq   = ZERO

  end subroutine kill

end module sphere_class
