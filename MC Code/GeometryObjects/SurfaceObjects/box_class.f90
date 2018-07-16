module box_class
  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError
  use vector_class,      only : vector
  use dictionary_class,  only : dictionary
  use surface_inter,     only : surface, printSurfDef, surfaceShelf
  use xPlane_class,      only : xPlane
  use yPlane_class,      only : yPlane
  use zPlane_class,      only : zPlane

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'box'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.

  !!
  !! Constructor
  !!
  interface box
    module procedure box_fromDict
  end interface

  !!
  !! Cuboid alligned with coordinate axies
  !!  *Planes(1) is above *Planes(2) on the *-axis i.e.
  !!    xPlanes(1) = 1.0 & xPlanes(2) = 0.1
  !!    Indexes of the surfaces in the surfaceShelf are in the same order
  !!
  type, public, extends(surface) :: box
    private
    !type(xPlane), dimension(:), pointer :: xPlanes => null()
    !type(yPlane), dimension(:), pointer :: yPlanes => null()
    !type(zPlane), dimension(:), pointer :: zPlanes => null()

    !! Planes and plane indexes
    type(xPlane), dimension(2)      :: xPlanes
    type(yPlane), dimension(2)      :: yPlanes
    type(zPlane), dimension(2)      :: zPlanes
    integer(shortInt), dimension(6) :: planeIdx  = [1,2,3,4,5,6]


    !! Boundary conditions in order x1, x2, y1, y2, z1, z2
    integer(shortInt), dimension(6) :: BC = noBC

    real(defReal), dimension(3)         :: a    ! the half-width in each direction of the cubic box
    real(defReal), dimension(3)         :: origin

  contains
    ! Initialisation & Indentification procedures
    procedure :: init
    procedure :: type
    procedure :: getDef
    procedure :: cannotBeBoundary
    procedure :: setBoundaryConditions

    ! Runtime procedures
    procedure :: evaluate
    procedure :: distance
    procedure :: normalVector
    procedure :: boundaryTransform

  end type box

contains

  !!
  !! Initialise the box as six plane surfaces
  !! Employ the convention that the first plane
  !! is in front of the origin, the second behind
  !!
  subroutine init(self, origin, a, id)
    class(box), intent(inout)               :: self
    real(defReal), dimension(3), intent(in) :: a, origin
    integer(shortInt), intent(in)           :: id
    integer(shortInt)                       :: i
    real(defReal)                           :: neg
    character(100),parameter :: Here ='init( box_class.f90)'

    self % origin = origin
    self % a = a

    if(any(a < surface_tol)) call fatalError(Here,'Box dimensions must be greater than surface tolerance & +ve ')

    call self % setId(id)

    ! Initialise each plane in each cardinal direction - Top planes
    call self % xPlanes(1) % init(origin(1) + a(1),1)
    call self % yPlanes(1) % init(origin(2) + a(2),1)
    call self % zPlanes(1) % init(origin(3) + a(3),1)

    ! Initialise each plane in each cardinal direction - Bottom planes
    call self % xPlanes(1) % init(origin(1) + a(1),1)
    call self % yPlanes(1) % init(origin(2) + a(2),1)
    call self % zPlanes(1) % init(origin(3) + a(3),1)

    ! Load planes no the shelf and obtain indices
    ! *** IMPLEMENT ***

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns an initialised instance of box from dictionary and name
  !!
  function box_fromDict(dict) result(new)
    class(dictionary), intent(in) :: dict
    type(box)                     :: new
    integer(shortInt)             :: id
    real(defReal),dimension(3)    :: halfwidth, origin
    character(100),parameter :: Here ='box_fromDict ( box_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    halfwidth = dict % getRealArray('halfwidth')
    origin = dict % getRealArray('origin')

    call new % init(origin, halfwidth, id)

  end function box_fromDict

  !!
  !! Evaluate the surface function of the box
  !!
  elemental subroutine evaluate(self,res, r)
    class(box), intent(in)                  :: self
    real(defReal), intent(out)              :: res
    type(vector), intent(in)                :: r
    type(vector)                            :: r_bar

    ! Calculate position in frame at box origin
    r_bar = r - self % origin

    ! Calculate absolute perpendicular distance from the boundary in each cardinal direction
    r_bar = abs(r_bar % v) - self % a

    ! Point is inside the box or on the surface if all r_bar < surface_tol
    if (all(r_bar % v < surface_tol)) then
      ! Return -1 * closest perpendicular distance to the surface
      ! If point is on the surface its magniture is smaller then surface_tol
      res = -minval(abs(r_bar %v))

    else ! Point is outside
      res = maxval(r_bar % v)    ! Return minimum distance

    end if

  end subroutine evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  elemental function type(self)
    class(box), intent(in) :: self
    character(nameLen)     :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(box), intent(in)                 :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % a, self % origin])

  end subroutine getDef

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(box), intent(in) :: self
    logical(defBool)       :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate the distance to the nearest surface of the box
  !!
  !! Use so called slab algorithm
  !! Google: "slab algorithm ray tracing"
  !!
  !! 1) Calculate nearest and furthest intersection point with planes in cardinal direction
  !! 2) Move to next cardinal direction and calculate new nerest and furthest distance
  !! 3) Save maximum nearest distance and minimum furthest distance
  !! 4) Go through all cardinal directions (X,Y,Z)
  !! 5) If at the end furthest < nearest or both furthest and nearest are -ve ray is a miss
  !! 6) Otherwise determine smallest distance to surface and which plane was X-ed
  !!
  !! Returns index of approperiate plane if it was X-ed
  !! If there was not intersection returns box own idx
  !!
  elemental subroutine distance(self, dist, idx, r, u)
    class(box), intent(in)           :: self
    real(defReal), intent(out)       :: dist
    integer(shortInt), intent(out)   :: idx
    type(vector), intent(in)         :: r
    type(vector), intent(in)         :: u
    real(defReal),dimension(3)       :: a_near, a_far
    type(vector)                     :: r_bar
    real(defReal)                    :: t_near, t_far
    real(defReal)                    :: test_near, test_far
    integer(shortInt)                :: near_ax, far_ax
    integer(shortInt), parameter     :: X_axis = 1, Y_axis = 2, Z_axis = 3
    logical(defBool)                 :: farNeg, nearNeg



    ! Transform to frame centered at the slab and choose nearest and furthest vertex
    ! based on octant of direction vector
    r_bar = r - self % origin
    a_far   = sign(self % a, u % v)
    a_near = -a_far

    ! Initialise test distances
    t_near = -huge(t_near)
    t_far  =  huge(t_far)

    ! Intersection with X planes
    if (u % v(1) /= ZERO) then
      test_near = (a_near(1) - r_bar % v(1)) / u % v(1)
      test_far  = (a_far(1)  - r_bar % v(1)) / u % v(1)

      ! Find max(t_near,test_near) & save Axis
      if(t_near < test_near) then
        t_near  = test_near
        near_ax = X_axis
      end if

      ! Find min(t_far,test_far) & save Axis
      if(t_far > test_far) then
        t_far = test_far
        far_ax = X_axis
      end if

    end if

    ! Intersection with Y planes
    if (u % v(2) /= ZERO) then
      test_near = (a_near(2) - r_bar % v(2)) / u % v(2)
      test_far  = (a_far(2) - r_bar % v(2)) / u % v(2)

      ! Find max(t_near,test_near) & save Axis
      if(t_near < test_near) then
        t_near  = test_near
        near_ax = Y_axis
      end if

      ! Find min(t_far,test_far) & save Axis
      if(t_far > test_far) then
        t_far = test_far
        far_ax = Y_axis
      end if

    end if

    ! Intersection with Z planes
    if (u % v(3) /= ZERO) then
      test_near = (a_near(3) - r_bar % v(3)) / u % v(3)
      test_far  =  (a_far(3) - r_bar % v(3)) / u % v(3)

      ! Find max(t_near,test_near) & save Axis
      if(t_near < test_near) then
        t_near  = test_near
        near_ax = Z_axis
      end if

      ! Find min(t_far,test_far) & save Axis
      if(t_far > test_far) then
        t_far = test_far
        far_ax = Z_axis
      end if

    end if

    ! Precalculate if t_far & t_near are -ve
    ! We want to account for the direction of motion when evaluating surface_tol
    ! surface_tol is defined as perpendicular distance, thus we need to divide by cosine
    farNeg  = t_far  < surface_tol / abs(u % v(far_ax))
    nearNeg = t_near < surface_tol / abs(u % v(near_ax))

    ! There is no intersection when either t_far < t_near or both T_far & t_near are -ve
    if( (t_far < t_near) .or. (farNeg .and. nearNeg)) then
      dist = INFINITY
      idx = self % myIdx()

    else
      ! t_near < surf_tol => point is inside surface -> take t_far
      ! When returning make shure to not return dist > INFINITY
      if ( nearNeg) then
        dist = min(t_far,INFINITY)

        ! Calculate plane code and return approperiate Idx
        idx = 2 * far_ax - 1
        if (u % v(far_ax) < ZERO) idx = idx +1
        idx = self % planeIdx(idx)

      else
        dist = min(t_near,INFINITY)

        ! Calculate plane code and return approperiate Idx
        idx = 2 * near_ax - 1
        if (u % v(near_ax) > ZERO) idx = idx +1
        idx = self % planeIdx(idx)

      end if
    end if
  end subroutine distance

  !!
  !! Return normalised vector normal to the surface
  !! Lo and behold! Without single if statement!
  !! Can be reduced to 3 lines in one doesn't care about readability.
  !!
  !! Effectivly returns normal to the closest face of the box (not the surface that contains the face!)
  !!
  !! Transforms position into all +ve octant to make everything easy.
  !! Space is then divided by drawing 45 deg plane through each edge to divide space into three part
  !! Point is transformed so vertex of the octant is the new origin
  !! Normal vector in the direction of maximum component the resulting vector build
  !! Normal vector is transformed back to the original octant
  !!
  elemental function normalVector(self, r) result(normal)
    class(box), intent(in)      :: self
    type(vector), intent(in)    :: r
    type(vector)                :: r_bar, sig_r
    integer(shortInt)           :: maxCom
    type(vector)                :: normal

    ! Move to approperiate octant. And save sign vector to be able to return
    sig_r = sign([ONE, ONE, ONE],r % v)
    r_bar = abs(r % v)

    ! Calculate vector between the octant corner and current point
    r_bar = r_bar - self % a

    ! Identify surface in the octant as the one with the largest component
    maxCom = maxloc(r_bar % v,1)

    ! Build normal vector
    normal = [ZERO, ZERO, ZERO]
    normal % v(maxCom) = ONE

    ! Multiply resulting vector by sign vector to return to original octant
    normal = normal % v * sig_r % v

  end function normalVector

  !!
  !! Apply generic boundary conditions to the box
  !!
  subroutine setBoundaryConditions(self, BC)
    class(box), intent(inout)                   :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditionsBox( box_class.f90)'
!
!    if (size(BC) < 6) call fatalError(Here,'Wrong size of BC string. Must be at least 6')
!
!    ! Positive x boundary
!    if(BC(1) == vacuum) then
!      self % xPlanes(1) % isVacuum = .TRUE.
!    else if(BC(1) == reflective) then
!      self % xPlanes(1) % isReflective = .TRUE.
!    else if(BC(1) == periodic) then
!      if(BC(2) /= periodic) then
!        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
!      else
!        self % xPlanes(1) % isPeriodic = .TRUE.
!        self % xPlanes(1) % periodicTranslation = [-TWO*self % a(1), ZERO, ZERO]
!      end if
!    else
!      call fatalError(Here,'Invalid boundary condition provided')
!    end if
!
!    ! Negative x boundary
!    if(BC(2) == vacuum) then
!      self % xPlanes(2) % isVacuum = .TRUE.
!    else if(BC(2) == reflective) then
!      self % xPlanes(2) % isReflective = .TRUE.
!    else if(BC(2) == periodic) then
!      if(BC(1) /= periodic) then
!        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
!      else
!        self % xPlanes(2) % isPeriodic = .TRUE.
!        self % xPlanes(2) % periodicTranslation = [TWO*self % a(1), ZERO, ZERO]
!      end if
!    else
!      call fatalError(Here,'Invalid boundary condition provided')
!    end if
!
!    ! Positive y boundary
!    if(BC(3) == vacuum) then
!      self % yPlanes(1) % isVacuum = .TRUE.
!    else if(BC(3) == reflective) then
!      self % yPlanes(1) % isReflective = .TRUE.
!    else if(BC(3) == periodic) then
!      if(BC(4) /= periodic) then
!        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
!      else
!        self % yPlanes(1) % isPeriodic = .TRUE.
!        self % yPlanes(1) % periodicTranslation = [ZERO, -TWO*self % a(2), ZERO]
!      end if
!    else
!      call fatalError(Here,'Invalid boundary condition provided')
!    end if
!
!    ! Negative y boundary
!    if(BC(4) == vacuum) then
!      self % yPlanes(2) % isVacuum = .TRUE.
!    else if(BC(4) == reflective) then
!      self % yPlanes(2) % isReflective = .TRUE.
!    else if(BC(4) == periodic) then
!      if(BC(3) /= periodic) then
!        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
!      else
!        self % yPlanes(2) % isPeriodic = .TRUE.
!        self % yPlanes(2) % periodicTranslation = [ZERO, TWO*self % a(2), ZERO]
!      end if
!    else
!      call fatalError(Here,'Invalid boundary condition provided')
!    end if
!
!    ! Positive z boundary
!    if(BC(5) == vacuum) then
!      self % zPlanes(1) % isVacuum = .TRUE.
!    else if(BC(5) == reflective) then
!      self % zPlanes(1) % isReflective = .TRUE.
!    else if(BC(5) == periodic) then
!      if(BC(6) /= periodic) then
!        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
!      else
!        self % zPlanes(1) % isPeriodic = .TRUE.
!        self % zPlanes(1) % periodicTranslation = [ZERO, ZERO, -TWO*self % a(3)]
!      end if
!    else
!      call fatalError(Here,'Invalid boundary condition provided')
!    end if
!
!    ! Negative z boundary
!    if(BC(6) == vacuum) then
!      self % zPlanes(2) % isVacuum = .TRUE.
!    else if(BC(6) == reflective) then
!      self % zPlanes(2) % isReflective = .TRUE.
!    else if(BC(6) == periodic) then
!      if(BC(5) /= periodic) then
!        call fatalError(Here, 'Both positive and negative boundary conditions must be periodic')
!      else
!        self % zPlanes(2) % isPeriodic = .TRUE.
!        self % zPlanes(2) % periodicTranslation = [ZERO, ZERO, TWO*self % a(3)]
!      end if
!    else
!      call fatalError(Here,'Invalid boundary condition provided')
!    end if

  end subroutine setBoundaryConditions

  !! *** NEED to be reworked
  !! Check halfspaces to identify which combination of
  !! boundary conditions to apply to a point
  !!
  subroutine boundaryTransform(self, r, u)
    class(box), intent(in)                     :: self
    type(vector), intent(inout)                :: r
    type(vector), intent(inout)                :: u
    logical(defBool)                           :: front, back, left, right, above, below
    character(100),parameter :: Here ='boundaryTransform( box_class.f90)'

!    ! Point can be within one of 27 regions
!    ! Locate which region and then apply BCs as appropriate
!    front = self % xPlanes(1) % halfspace(r, u)
!    back = .NOT. self % xPlanes(2) % halfspace(r, u)
!    left = self % yPlanes(1) % halfspace(r, u)
!    right = .NOT. self % yPlanes(2) % halfspace(r, u)
!    above = self % zPlanes(1) % halfspace(r, u)
!    below = .NOT. self % zPlanes(2) % halfspace(r, u)
!
!    ! Point is in the upper front left corner
!    if (front .AND. left .AND. above) then
!      call self % xPlanes(1) % boundaryTransform(r, u)
!      call self % yPlanes(1) % boundaryTransform(r, u)
!      call self % zPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is in the upper back left corner
!    else if (back .AND. left .AND. above) then
!      call self % xPlanes(2) % boundaryTransform(r, u)
!      call self % yPlanes(1) % boundaryTransform(r, u)
!      call self % zPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is above and to the left
!    else if (left .AND. above) then
!      call self % yPlanes(1) % boundaryTransform(r, u)
!      call self % zPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is in the upper front right corner
!    else if (front .AND. right .AND. above) then
!      call self % xPlanes(1) % boundaryTransform(r, u)
!      call self % yPlanes(2) % boundaryTransform(r, u)
!      call self % zPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is in the upper back right corner
!    else if (back .AND. right .AND. above) then
!      call self % xPlanes(2) % boundaryTransform(r, u)
!      call self % yPlanes(2) % boundaryTransform(r, u)
!      call self % zPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is above and to the right
!    else if (right .AND. above) then
!      call self % yPlanes(2) % boundaryTransform(r, u)
!      call self % zPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is above
!    else if (above) then
!      call self % zPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is above and in front
!    else if (front .AND. above) then
!      call self % xPlanes(1) % boundaryTransform(r, u)
!      call self % zPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is above and behind
!    else if (back .AND. above) then
!      call self % xPlanes(2) % boundaryTransform(r, u)
!      call self % zPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is in the lower front left
!    else if (front .AND. left .AND. below) then
!      call self % xPlanes(1) % boundaryTransform(r, u)
!      call self % yPlanes(1) % boundaryTransform(r, u)
!      call self % zPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is in the lower back left
!    else if (back .AND. left .AND. below) then
!      call self % xPlanes(2) % boundaryTransform(r, u)
!      call self % yPlanes(1) % boundaryTransform(r, u)
!      call self % zPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is in the lower left
!    else if (left .AND. below) then
!      call self % yPlanes(1) % boundaryTransform(r, u)
!      call self % zPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is in the lower front right
!    else if (front .AND. right .AND. below) then
!      call self % xPlanes(1) % boundaryTransform(r, u)
!      call self % yPlanes(2) % boundaryTransform(r, u)
!      call self % zPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is in the lower back right
!    else if (back .AND. right .AND. below) then
!      call self % xPlanes(2) % boundaryTransform(r, u)
!      call self % yPlanes(2) % boundaryTransform(r, u)
!      call self % zPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is in the lower right
!    else if (right .AND. below) then
!      call self % yPlanes(2) % boundaryTransform(r, u)
!      call self % zPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is below
!    else if (below) then
!      call self % zPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is in the front left
!    else if (front .AND. left) then
!      call self % xPlanes(1) % boundaryTransform(r, u)
!      call self % yPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is in the back left
!    else if (back .AND. left) then
!      call self % xPlanes(2) % boundaryTransform(r, u)
!      call self % yPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is on the left
!    else if (left) then
!      call self % yPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is on the front right
!    else if (front .AND. right) then
!      call self % xPlanes(1) % boundaryTransform(r, u)
!      call self % yPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is on the back right
!    else if (back .AND. right) then
!      call self % xPlanes(2) % boundaryTransform(r, u)
!      call self % yPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is on the front
!    else if (front) then
!      call self % xPlanes(1) % boundaryTransform(r, u)
!
!    ! Point is on the back
!    else if (back) then
!      call self % xPlanes(2) % boundaryTransform(r, u)
!
!    ! Point is in the box
!    else
!      call fatalError(Here,'Cannot apply boundary transformation: point is already within the boundary')
!    end if

  end subroutine boundaryTransform

end module box_class
