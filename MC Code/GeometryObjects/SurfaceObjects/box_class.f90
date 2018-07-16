module box_class
  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, numToChar
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
  integer(shortInt),parameter  :: ABOVE = 1, IN = 0, BELOW = -1

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
    integer(shortInt), dimension(1:6) :: BC = noBc

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

    ! Private procedures
    procedure, private :: applyBC

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
    logical(defBool), dimension(3)              :: top_periodic
    character(100),parameter :: Here ='setBoundaryConditionsBox( box_class.f90)'

    ! Load BC codes
    self % BC = BC(1:6)

    ! Verify BC periodic BC
    !*** Sorry for the oneliner but TBH its not that horrible.
    !*** I'm just to tired to name a halper variable in a good way... -MAK
    if(.not.all( (self % BC([1,3,5] ) == periodic) .eqv. (self % BC([2,4,6]) == periodic))) then
      call fatalError(Here,'Periodic BC need to be applied to oposite surfaces')

    end if

  end subroutine setBoundaryConditions

  !!
  !! Apply boundary conditions specified by BC array given at initialisation
  !!
  subroutine boundaryTransform(self, r, u)
    class(box), intent(in)               :: self
    type(vector), intent(inout)          :: r
    type(vector), intent(inout)          :: u
    real(defReal), dimension(3)          :: pos
    real(defReal), dimension(3)          :: a_bar
    integer(shortInt), dimension(3)      :: Ri
    integer(shortInt)                    :: N_trans, i, dir,flag
    character(100),parameter :: Here ='boundaryTransform( box_class.f90)'

    ! Calculate dimension reduced by surface tolerance
    ! This is needed to apply BC for particles at the surface
    a_bar = self % a - surface_tol

    ! Calculate maximum distance in terms of coordinate transforms N_trans
    Ri = ceiling(abs(r % v - self % origin)/ a_bar)/2
    N_trans = max(Ri(1), Ri(2), Ri(3))

    ! Loop over number of transformation required to fold back to inside geometry
    do i=1,N_trans

      ! Calculate position of the point wrt to orgin
      pos = r % v - self % origin

      ! Loop over cardinal directions and apply BC
      do dir =1,3
        if( pos(dir) > a_bar(dir) ) then       ! ABOVE BOX
          flag = ABOVE

        else if( pos(dir) < -a_bar(dir) ) then ! BELOW BOX
          flag = BELOW

        else
          flag = IN

        end if
        ! Apply BC
        call self % applyBC(r,u,flag,dir)

      end do
    end do
  end subroutine boundaryTransform

  !!
  !! Macro to apply BC given direction and flag which specifies position wrt planes
  !!
  subroutine applyBC(self,r,u,flag,dir)
    class(box), intent(in)        :: self
    type(vector), intent(inout)   :: r
    type(vector), intent(inout)   :: u
    integer(shortInt), intent(in) :: flag
    integer(shortInt), intent(in) :: dir
    integer(shortInt)             :: BC, idx
    real(defReal)                 :: off, disp
    character(100),parameter      :: Here ='applyBC( box_class.f90)'

    ! Flag can be 1, 0, -1 corresponding to Above, Inside, Below
    ! dir can be 1, 2, 3 and specifies cardinal direction

    ! Note that flag 0 is eqivalent to application of vacuum BC

    ! Recover BC
    if (flag == IN) then
      BC = vacuum

    else
      idx = dir * 2
      if( flag == ABOVE) idx = idx - 1
      BC = self % BC(idx)

    end if

    ! Apply BC
    select case(BC)
      case(vacuum)
        ! Do nothing

      case(reflective)
        ! Find point offset
        off = self % a(dir)
        if (flag == BELOW) off = -off
        off = self % origin(dir) + off

        ! Calculate displacement
        disp = r % v(dir) - off

        ! Perform reflection
        r % v(dir) = r % v(dir) - TWO * disp
        u % v(dir) = -u % v(dir)

      case(periodic)
        ! Calculate offset
        off = self % a(dir)
        if (flag == BELOW) off = -off

        ! Perform transition
        r % v(dir) = r % v(dir) - TWO * off

      case default
        call fatalError(Here,'Unrecognised BC: '// numToChar(BC) )

    end select
  end subroutine applyBC

end module box_class
