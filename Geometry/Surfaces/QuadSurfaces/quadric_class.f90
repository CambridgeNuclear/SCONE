module quadric_class

  use numPrecision
  use universalVariables, only : SURF_TOL, INF, X_AXIS, Y_AXIS, Z_AXIS
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill

  implicit none
  private

  !!
  !! General quadratic surface - also known as a quadric
  !!
  !! F(x,y,z) = Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fzx + Gx + Hy + Iz + J
  !!
  !! Surface tolerance: 2 * max(coeffs) * SURF_TOL
  !!
  !! Sample dictionary input:
  !!   quad { type quadric; 
  !!         id 3;
  !!         coeffs (1 -2.3 0.4 20 7.7777 0.004 6E9 -4.20 0 11);  
  !!         // correspond to(A B C D E F G H I J)
  !!        }
  !!
  !! Private Members:
  !!   coeffs -> Quadric coefficients
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: quadric
    private
    real(defReal), dimension(10) :: coeffs = ZERO
  contains
    ! Superclass procedures
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill

  end type quadric

contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(quadric), intent(in) :: self
    character(:), allocatable  :: str

    str = 'quadric'

  end function myType

  !!
  !! Initialise quadric from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!   fatalError if id < 0 or incorrect size of coeffs
  !!
  subroutine init(self, dict)
    class(quadric), intent(inout)            :: self
    class(dictionary), intent(in)            :: dict
    integer(shortInt)                        :: id
    real(defReal), dimension(:), allocatable :: coeffs
    character(100), parameter :: Here = 'init (quadric_class.f90)'

    ! Get from dictionary
    call dict % get(id, 'id')
    call dict % get(coeffs, 'coeffs')
    
    ! Check ID validity
    if (id < 1) call fatalError(Here,'Invalid surface id provided. ID must be > 1')

    ! Check origin size
    if (size(coeffs) /= 10) then
      call fatalError(Here,'coeffs needs to have size 10. Has: '//numToChar(size(coeffs)))
    end if

    call self % setID(id)
    self % coeffs = coeffs
   
    ! Set surface tolerance - what should this be?
    ! Could choose alternatively to be equivalent to the
    ! sphere case, i.e., 2 * coeffs(10) = 2*R
    call self % setTol( TWO * maxval(abs(coeffs)) * SURF_TOL)

  end subroutine init

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! Not attempted - may be unbounded, depending on the coefficients.
  !! Would generally require checking the Jacobian of the surface and
  !! finding critical points.
  !! I think this only matters if we wanted to use this as a boundary.
  !!
  !! TODO: This
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(quadric), intent(in)  :: self
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
    class(quadric), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal)                           :: c

    associate(a => self % coeffs, x => r(1), y => r(2), z => r(3))
      c = x * (a(1)*x + a(4)*y + a(6)*z + a(7)) + &
            y * (a(2)*y + a(5)*z + a(8)) + &
            z * (a(3)*z + a(9)) + a(10)
    end associate

  end function evaluate

  !!
  !! Return distance to the surface
  !!
  !! See surface_inter for details
  !!
  !! Converts the general quadric to a 1D quadratic to solve
  !!   ad^2 + 2kd + c = 0
  !!   c = F(r)
  !!   k = A*u1*r1 + B*u2*r2 + C*u3*r3 + 
  !!       (D(u1*r2 + u2*r1) + E(u2*r3 + u3*r2) + F(u3*r1 + u1*r3)
  !!        + G*u1 + H*u2 + I*u3) / 2
  !!   a = A*u1*u1 + B*u2*u2 + C*u3*u3 + D*u1*u2 + E*u2*u3 + F*u1*u3
  !!
  pure function distance(self, r, u) result(d)
    class(quadric), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    real(defReal)                           :: a, c, k, delta

    c = self % evaluate(r)
    associate(cf => self % coeffs, x => r(1), y => r(2), z => r(3), &
                    u1 => u(1), u2 => u(2), u3 => u(3))
      k = cf(1)*u1*x + cf(2)*u2*y + cf(3)*u3*z + &
              HALF * (cf(4) * (u1*y + u2*x) + &
              cf(5) * (u2*z + u3*y) + cf(6) * (u3*x + u1*z) + &
              cf(7) * u1 + cf(8)*u2 + cf(9)*u3)
      a = u1 * (cf(1) * u1 + cf(4) * u2 + cf(6) * u3) + &
              u2 * (cf(2) * u2 + cf(5) * u3) + &
              cf(3) * u3 * u3
    end associate
    delta = k*k - a*c  ! Technically delta/4 - quadratic discriminant

    ! Calculate the distance
    if (delta < ZERO .or. a == ZERO) then ! No intersection
      d = INF

    else if (abs(c) < self % surfTol()) then ! Point at a surface
      if ( k >= ZERO) then
        d = INF
      else
        d = -k + sqrt(delta)
        d = d/a
      end if

    else if (c < ZERO) then ! Point inside the surface
      d = -k + sqrt(delta)
      d = d/a

    else ! Point outside the surface
      d = -k - sqrt(delta)
      d = d/a
      if (d <= ZERO) d = INF

    end if

    ! Cap distance at Infinity
    d = min(d, INF)

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !!
  !! See surface_inter for details
  !!
  pure function going(self, r, u) result(halfspace)
    class(quadric), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace
    real(defReal), dimension(3)             :: grad

    associate(a => self % coeffs)
      grad(1) = TWO * a(1) * r(1) + a(4) * r(2) + a(5) * r(3)
      grad(2) = TWO * a(2) * r(2) + a(4) * r(1) + a(6) * r(3)
      grad(3) = TWO * a(3) * r(3) + a(5) * r(1) + a(6) * r(2)
    end associate
    halfspace = dot_product(grad, u) >= ZERO

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(quadric), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % coeffs = ZERO

  end subroutine kill

end module quadric_class
