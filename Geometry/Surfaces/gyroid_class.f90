module gyroid_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar, swap
  use dictionary_class,   only : dictionary
  use surface_inter,      only : surface, kill_super => kill

  implicit none
  private

  !!
  !! Gyroid surface, a type of triply periodic minimal surface
  !!
  !! Satisfies the equation:
  !! f(x,y,z) = cos(2*pi*a*(x-x0)) * sin(2*pi*b*(y-y0)) + cos(2*pi*b*(y-y0)) * sin(2*pi*d*(z-z0))
  !!           + cos(2*pi*d*(z-z0)) * sin(2*pi*a*(x-x0)) - c = 0
  !!
  !! Sample Dictionary Input:
  !!   myGyroid { type fractal; id 92; origin (-1 -2 4.5); freq (1.0 -3.0 7); c 0.0; }
  !!
  !! Boundary Conditions:
  !!   Does not support boundary conditions.
  !!
  !! Private Members:
  !!   origin -> x, y, and z co-ordinates defining the centre
  !!   freq   -> spatial frequencies in x, y, and z
  !!   c -> 
  !!   BC -> Boundary conditions - not supported
  !!
  !! Interface:
  !!   surface interface
  !!
  type, public, extends(surface) :: gyroid
    private
    real(defReal), dimension(3)     :: origin = ZERO
    real(defReal), dimension(3)     :: freq = ZERO
    real(defReal)                   :: c  = ZERO
    integer(shortInt), dimension(6) :: BC = VACUUM_BC

  contains
    ! Superclass procedures
    procedure :: myType
    procedure :: init
    procedure :: boundingBox
    procedure :: evaluate
    procedure :: distance
    procedure :: going
    procedure :: kill
  end type gyroid


contains

  !!
  !! Return surface type name
  !!
  !! See surface_inter for more details
  !!
  pure function myType(self) result(str)
    class(gyroid), intent(in) :: self
    character(:), allocatable  :: str

    str = 'gyroid'

  end function myType

  !!
  !! Initialise box from a dictionary
  !!
  !! See surface_inter for more details
  !!
  !! Errors:
  !!
  subroutine init(self, dict)
    class(gyroid), intent(inout)  :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt)             :: id, N
    real(defReal), dimension(:), allocatable :: temp
    character(100), parameter :: Here = 'init (gyroid_class.f90)'

    ! Load id
    call dict % get(id,'id')
    if (id <= 0) call fatalError(Here, 'ID must be <=0. Is: '//numToChar(id))
    call self % setID(id)
    
    ! Load origin
    call dict % get(temp,'origin')
    N = size(temp)
    if (N /= 3) call fatalError(Here,'origin must have size 3. Has: '//numToChar(N))
    self % origin = temp

    ! Load spatial frequencies
    call dict % get(temp,'freq')
    N = size(temp)
    if (N /= 3) call fatalError(Here,'freq must have size 3. Has: '//numToChar(N))
    ! Scale by 2*pi for later simplicity
    temp = temp * TWO_PI
    self % freq = temp

    ! Load shift
    call dict % get(self % c,'c')

  end subroutine init

  !!
  !! Return axis-aligned bounding box for the surface
  !!
  !! See surface_inter for details
  !!
  pure function boundingBox(self) result(aabb)
    class(gyroid), intent(in)   :: self
    real(defReal), dimension(6) :: aabb

    aabb(1:3) = -INF
    aabb(4:6) = INF

  end function boundingBox

  !!
  !! Evaluate surface expression c = F(r)
  !!
  !! See surface_inter for details
  !!
  pure function evaluate(self, r) result(f)
    class(gyroid), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3)             :: rb
    real(defReal)                           :: f
    integer(shortInt)                       :: i

    ! Move to origin-frame and scale by frequency
    rb = (r - self % origin) * self % freq

    f = cos(rb(1)) * sin(rb(2)) + &
            cos(rb(2)) * sin(rb(3)) + &
            cos(rb(3)) * sin(rb(1)) - self % c
    
  end function evaluate

  !!
  !! Return distance to the surface
  !! Not defined for the gyroid surface which should only be used with delta tracking.
  !!
  !! See surface_inter for details
  !!
  pure function distance(self, r, u) result(d)
    class(gyroid), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    real(defReal)                           :: d
    real(defReal)                           :: tol, dSmall, dInc
    integer(shortInt)                       :: it
    real(defReal), dimension(3)             :: rb, invFreq, coss, sins, grad
    
    ! Need to do a root-finding of some sort - so we iterate
    dInc = INF
    tol = SURF_TOL

    d = ZERO

    invFreq = TWO_PI / (self % freq)
    dSmall = 0.01_defReal * dot_product(invFreq, u)

    r0 = r
    it = 0

    do while(dInc > tol .and. it < self % maxIt)
    
      it = it + 1
      
      ! Move to origin-frame and scale by frequency
      rb = (r0 - self % origin) * self % freq

      ! Evaluate the function
      do i = 1, 3
        coss(i) = cos(rb(i))
        sins(i) = sin(rb(i))
      end do

      f = coss(1) * sins(2) + coss(2) * sins(3) + coss(3) * sins(1) - self % c
 
      ! Obtain the gradient
      do i = 1, 3   
        grad(i) = self % freq(i) * (coss(i) * coss(mod(i+1,3) + 1) - sins(i) * sins(mod(i,3) + 1)
      end do

      ! How close are we? Should we use ray marching or Newton-Raphson?

      if (doNewton) then

        dInc = -f / dot_product(grad,u)


      else


      end if
      
      d = d + dInc
      r0 = r0 + dInc * u

    end do

  end function distance

  !!
  !! Returns TRUE if particle is going into +ve halfspace
  !! Not defined for the gyroid surface which should only be used with delta tracking.
  !!
  !! See surface_inter for details
  !!
  pure function going(self, r, u) result(halfspace)
    class(gyroid), intent(in)               :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: halfspace

    halfspace = .false.

  end function going

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(gyroid), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Local
    self % origin = ZERO
    self % freq = ZERO
    self % c = ZERO

  end subroutine kill

end module gyroid_class
