module xSquareCylinder_class

  use numPrecision
  use universalVariables
  use genericProcedures, only : fatalError, numToChar
  use vector_class,      only : vector
  use dictionary_class,  only : dictionary

  use surface_inter,     only : surface, printSurfDef, surfaceShelf
  use yPlane_class,      only : yPlane
  use zPlane_class,      only : zPlane

  implicit none
  private

  !!
  !! Constants describing surface properties
  !!
  character(nameLen),parameter :: TYPE_NAME    = 'xSquareCylinder'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.
  integer(shortInt),parameter  :: ABOVE = 1, IN = 0, BELOW = -1

  !!
  !! Constructor
  !!
  interface xSquareCylinder
    module procedure xSquareCylinder_fromDict
  end interface

  !!
  !! Square cylinder aligned with x-axis
  !! ***
  !! All Square cylinders should be combined into single code with diffrent axis mapping
  !!
  type, public, extends(surface) :: xSquareCylinder
    private

    !! Planes and plane indexes
    integer(shortInt), dimension(6) :: planeIdx  = [1,2,3,4,5,6]

    !! Boundary conditions in order x1, x2, y1, y2, z1, z2
    integer(shortInt), dimension(6) :: BC = noBc

    real(defReal), dimension(3)         :: a                ! the half-width in each direction of the cylinder (0 in x)
    real(defReal), dimension(3)         :: origin = [ZERO, ZERO, ZERO]

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

  end type xSquareCylinder


contains

  !!
  !! Initialise the square cylinder from components
  !!
  subroutine init(self, origin, a, id)
    class(xSquareCylinder), intent(inout)   :: self
    real(defReal), dimension(3), intent(in) :: a
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in)           :: id

    self % a = a
    if((a(2) < surface_tol) .OR. (a(3) < surface_tol)) &
    call fatalError('init, xSquareCylinder','Widths must be greater than surface tolerance')
    self % origin = origin

    call self % setId(id)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns and initialised instance of xSquareCylinder from dictionary and name
  !!
  function xSquareCylinder_fromDict(dict) result(new)
    class(dictionary), intent(in)  :: dict
    type(xSquareCylinder)          :: new
    integer(shortInt)              :: id
    real(defReal), dimension(3)    :: origin, halfwidth
    character(100),parameter :: Here ='xSquareCylinder_fromDict ( xSquareCylinder_class.f90)'

    id = dict % getInt('id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    halfwidth = dict % getRealArray('halfwidth')
    origin = dict % getRealArray('origin')

    call new % init(origin, halfwidth, id)

  end function xSquareCylinder_fromDict

  !!
  !! Evaluate the surface function of the square cylinder
  !!
  elemental subroutine evaluate(self,res, r)
    class(xSquareCylinder), intent(in)      :: self
    real(defReal), intent(out)              :: res
    type(vector), intent(in)                :: r
    real(defReal),dimension(2)              :: r_bar

    ! Calculate position in frame at box origin (ignore x axis)
    r_bar = r % v(2:3) - self % origin(2:3)

    ! Calculate absolute perpendicular distance from the boundary in each cardinal direction(except x)
    r_bar = abs(r_bar) - self % a(2:3)

    ! Point is inside the box or on the surface if all r_bar < surface_tol
    if (all(r_bar < surface_tol)) then
      ! Return -1 * closest perpendicular distance to the surface
      ! If point is on the surface its magniture is smaller then surface_tol
      res = -minval(abs(r_bar))

    else ! Point is outside
      res = maxval(r_bar)    ! Return perpendicular distance

    end if

  end subroutine evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  elemental function type(self)
    class(xSquareCylinder), intent(in) :: self
    character(nameLen)                 :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(xSquareCylinder), intent(in)     :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % a, self % origin])

  end subroutine getDef

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(xSquareCylinder), intent(in) :: self
    logical(defBool)                   :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate the distance to the nearest surface of the box
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the box plane surfaces
  !!
  elemental subroutine distance(self, dist, idx, r, u)
    class(xSquareCylinder), intent(in) :: self
    real(defReal), intent(out)          :: dist
    integer(shortInt), intent(out)      :: idx
    type(vector), intent(in)            :: r
    type(vector), intent(in)            :: u
    real(defReal),dimension(3)          :: a_near, a_far
    type(vector)                        :: r_bar
    real(defReal)                       :: t_near, t_far
    real(defReal)                       :: test_near, test_far
    integer(shortInt)                   :: near_ax, far_ax
    integer(shortInt), parameter        :: X_axis = 1, Y_axis = 2, Z_axis = 3
    logical(defBool)                    :: farNeg, nearNeg

    ! Transform to frame centered at the slab and choose nearest and furthest vertex
    ! based on octant of direction vector
    r_bar = r - self % origin
    a_far   = sign(self % a, u % v)
    a_near = -a_far

    ! Initialise test distances
    t_near = -huge(t_near)
    t_far  =  huge(t_far)

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
  !!
  !!
  elemental function normalVector(self, r) result(normal)
    class(xSquareCylinder), intent(in)  :: self
    type(vector), intent(in)            :: r
    real(defReal),dimension(2)          :: sig_r, r_bar
    integer(shortInt)                   :: maxCom
    type(vector)                        :: normal

    ! Move to approperiate octant. And save sign vector to be able to return
    sig_r = sign( [ONE, ONE],r % v(2:3))
    r_bar = abs(r % v(2:3))

    ! Calculate vector between the octant corner and current point
    r_bar = r_bar - self % a(2:3)

    ! Identify surface in the octant as the one with the largest component
    maxCom = maxloc(r_bar ,1)

    ! Translate maxCom to cmponent of 3-D vector
    ! Use select case for clarity. (Will have the same form for problematic ySquareCylinder)
    ! Hope that compiler makes it branchless (it should)
    select case(maxCom)
      case(1)
        maxCom = 2

      case(2)
        maxCom = 3

    end select

    ! Build normal vector
    normal = [ZERO ,ZERO, ZERO]
    normal % v(maxCom) = ONE

    ! Multiply resulting vector by sign vector to return to original octant
    normal % v(2:3) = normal % v(2:3) * sig_r

  end function normalVector


  !!
  !! Set boundary conditions for an xSquareCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditions(self, BC)
    class(xSquareCylinder), intent(inout)       :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions ( xSquareCylinder_class.f90)'

    if(size(BC) < 6) call fatalError(Here,'Wrong size of BC string. Must be at least 6')

    ! Load BC codes
    self % BC = BC(1:6)

    ! Verify BC periodic BC
    ! *** I am not apologising this time! - MAK
    if(.not.all( (self % BC([3,5] ) == periodic) .eqv. (self % BC([4,6]) == periodic))) then
      call fatalError(Here,'Periodic BC need to be applied to both oposite surfaces')

    end if


  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u)
    class(xSquareCylinder), intent(in)   :: self
    type(vector), intent(inout)          :: r
    type(vector), intent(inout)          :: u
    real(defReal), dimension(3)          :: pos
    real(defReal), dimension(3)          :: a_bar
    integer(shortInt), dimension(3)      :: Ri
    integer(shortInt)                    :: N_trans, i, dir,flag
    character(100),parameter :: Here ='boundaryTransform ( xSquareCylinder_class.f90)'

    ! Calculate dimension reduced by surface tolerance
    ! This is needed to apply BC for particles at the surface
    a_bar = self % a - surface_tol

    ! Calculate maximum distance in terms of coordinate transforms N_trans

    Ri(2:3) = ceiling(abs(r % v(2:3)- self % origin(2:3))/ a_bar(2:3))/2
    N_trans = maxval(Ri(2:3))

    ! Loop over number of transformation required to fold back to inside geometry
    do i=1,N_trans

      ! Calculate position of the point wrt to orgin
      pos = r % v - self % origin

      ! Loop over cardinal directions and apply BC
      do dir =2,3
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
    class(xSquareCylinder), intent(in) :: self
    type(vector), intent(inout)        :: r
    type(vector), intent(inout)        :: u
    integer(shortInt), intent(in)      :: flag
    integer(shortInt), intent(in)      :: dir
    integer(shortInt)                  :: BC, idx
    real(defReal)                      :: off, disp
    character(100),parameter           :: Here ='applyBC( xSquareCylinder_class.f90)'

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

end module xSquareCylinder_class
