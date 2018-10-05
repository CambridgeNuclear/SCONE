module squareCylinder_class

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
  character(nameLen),parameter :: TYPE_NAME    = 'squareCylinder'
  logical(defBool),parameter   :: ACCEPTS_BC   = .true.
  integer(shortInt),parameter  :: ABOVE = 1, IN = 0, BELOW = -1

  !!
  !! Constructor interfaces
  !!
  public :: xSquareCylinder
  interface xSquareCylinder
    module procedure xSquareCylinder_fromDict
  end interface

  public :: ySquareCylinder
  interface ySquareCylinder
    module procedure ySquareCylinder_fromDict
  end interface

  public :: zSquareCylinder
  interface zSquareCylinder
    module procedure zSquareCylinder_fromDict
  end interface

  !!
  !! Square cylinder aligned with x-axis, y-axis or z-axis
  !!
  !! Square cylinder is infinate in one axis (like real cylinder!)
  !! Square cylinder is like real cylinder but its base is a ractangle
  !!
  !! Despite the fact that we need to store only two halfwidths a we store them in 3-D array
  !! in order to be able to switch between diffrently aligned squareCylinders by changing axis
  !! mapping in component "dirs"
  !!
  !! Dummy direction components of origin and a are set to ZERO to avoid duplicate definitions
  !! of same surface.
  !!
  type, public, extends(surface) :: squareCylinder
    private

    !! Plane indexes and array that contains direction mapping for diffrent cylinders
    integer(shortInt), dimension(6) :: planeIdx  = [1,2,3,4,5,6]
    integer(shortInt), dimension(2) :: dirs      = -88

    !! Boundary conditions in order x1, x2, y1, y2, z1, z2
    integer(shortInt), dimension(6) :: BC = noBc

    !! Halfwidths in each direction of the cylinder (only two approperiate one are used)
    !! Also cylinder origin vector
    real(defReal), dimension(3)         :: a
    real(defReal), dimension(3)         :: origin = [ZERO, ZERO, ZERO]

  contains
    ! Initialisation & Indentification procedures
    procedure :: init
    procedure :: type
    procedure :: getDef
    procedure :: boundingBox
    procedure :: cannotBeBoundary
    procedure :: setBoundaryConditions

    ! Runtime procedures
    procedure :: evaluate
    procedure :: distance
    procedure :: normalVector
    procedure :: boundaryTransform

    ! Private procedures
    procedure, private :: applyBC

  end type squareCylinder


contains

  !!
  !! Initialise the square cylinder from components
  !!
  subroutine init(self, origin, a, id, type)
    class(squareCylinder), intent(inout)    :: self
    real(defReal), dimension(3), intent(in) :: a
    real(defReal), dimension(3), intent(in) :: origin
    integer(shortInt), intent(in)           :: id
    character(*)                            :: type
    character(100), parameter  :: Here ='init (squareCylinder_class.f90)'

    select case(type)
      case('xSquareCylinder')
        ! Load half lengths & origin
        self % a = a
        self % origin = origin

        if((a(2) < surface_tol) .or. (a(3) < surface_tol)) then
          call fatalError(Here,'Y & Z widths must be greater than surface tolerance')
        end if

        ! Select directions and zero dummy components of a and origin
        self % dirs = [Y_AXIS,Z_AXIS]
        self % a(X_AXIS)      = ZERO
        self % origin(X_AXIS) = ZERO

      case('ySquareCylinder')
        ! Load half lengths & origin
        self % a = a
        self % origin = origin

        if((a(1) < surface_tol) .or. (a(3) < surface_tol)) then
          call fatalError(Here,'X & Z widths must be greater than surface tolerance')
        end if

        ! Select directions and zero dummy components of a and origin
        self % dirs = [X_AXIS,Z_AXIS]
        self % a(Y_AXIS)      = ZERO
        self % origin(Y_AXIS) = ZERO

      case('zSquareCylinder')
        ! Load half lengths & origin
        self % a = a
        self % origin = origin

        if((a(1) < surface_tol) .or. (a(2) < surface_tol)) then
          call fatalError(Here,'X & Z widths must be greater than surface tolerance')
        end if

        ! Select directions and zero dummy components of a and origin
        self % dirs = [X_AXIS,Y_AXIS]
        self % a(Z_AXIS)      = ZERO
        self % origin(Z_AXIS) = ZERO

      case default
        call fatalError(Here,'Unrecognised type of square cylinder: ' // type)

    end select

    call self % setId(id)

    ! Hash and store surface definition
    call self % hashSurfDef()

  end subroutine init

  !!
  !! Returns and initialised instance of squareCylinder along X from dictionary
  !!
  function xSquareCylinder_fromDict(dict) result(new)
    class(dictionary), intent(in)              :: dict
    type(squareCylinder)                       :: new
    integer(shortInt)                          :: id
    real(defReal), dimension(:),allocatable    :: origin, halfwidth
    character(100),parameter :: Here ='xSquareCylinder_fromDict ( squareCylinder_class.f90)'

    call dict % get(id, 'id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    call dict % get(halfwidth, 'halfwidth')
    call dict % get(origin, 'origin')

    if(size(halfwidth) /= 3) then
      call fatalError(Here,'Halfwidth must have size 3')
    else if(size(origin) /=3) then
      call fatalError(Here,'Origin must have size 3')
    end if

    call new % init(origin, halfwidth, id,'xSquareCylinder')

  end function xSquareCylinder_fromDict

  !!
  !! Returns and initialised instance of squareCylinder along Y from dictionary
  !!
  function ySquareCylinder_fromDict(dict) result(new)
    class(dictionary), intent(in)              :: dict
    type(squareCylinder)                       :: new
    integer(shortInt)                          :: id
    real(defReal), dimension(:),allocatable    :: origin, halfwidth
    character(100),parameter :: Here ='ySquareCylinder_fromDict ( squareCylinder_class.f90)'

    call dict % get(id, 'id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    call dict % get(halfwidth, 'halfwidth')
    call dict % get(origin, 'origin')

    if(size(halfwidth) /= 3) then
      call fatalError(Here,'Halfwidth must have size 3')
    else if(size(origin) /=3) then
      call fatalError(Here,'Origin must have size 3')
    end if

    call new % init(origin, halfwidth, id,'ySquareCylinder')

  end function ySquareCylinder_fromDict

  !!
  !! Returns and initialised instance of squareCylinder along Y from dictionary
  !!
  function zSquareCylinder_fromDict(dict) result(new)
    class(dictionary), intent(in)              :: dict
    type(squareCylinder)                       :: new
    integer(shortInt)                          :: id
    real(defReal), dimension(:),allocatable    :: origin, halfwidth
    character(100),parameter :: Here ='zSquareCylinder_fromDict ( squareCylinder_class.f90)'

    call dict % get(id, 'id')
    if(id < 1) call fatalError(Here,'Invalid surface id provided')

    call dict % get(halfwidth, 'halfwidth')
    call dict % get(origin, 'origin')

    if(size(halfwidth) /= 3) then
      call fatalError(Here,'Halfwidth must have size 3')
    else if(size(origin) /=3) then
      call fatalError(Here,'Origin must have size 3')
    end if

    call new % init(origin, halfwidth, id,'zSquareCylinder')

  end function zSquareCylinder_fromDict

  !!
  !! Evaluate the surface function of the square cylinder
  !!
  elemental subroutine evaluate(self,res, r)
    class(squareCylinder), intent(in)      :: self
    real(defReal), intent(out)              :: res
    type(vector), intent(in)                :: r
    real(defReal),dimension(2)              :: r_bar

    associate( di => self % dirs)
      ! Calculate position in frame at box origin (ignore one axis)
      r_bar = r % v(di) - self % origin(di)

      ! Calculate absolute perpendicular distance from the boundary
      ! in each cardinal direction except one axis
      r_bar = abs(r_bar) - self % a(di)

      ! Point is inside the box or on the surface if all r_bar < surface_tol
      if (all(r_bar < surface_tol)) then
        ! Return -1 * closest perpendicular distance to the surface
        ! If point is on the surface its magniture is smaller then surface_tol
        res = -minval(abs(r_bar))

      else ! Point is outside
        res = maxval(r_bar)    ! Return perpendicular distance

      end if
    end associate
  end subroutine evaluate

  !!
  !! Return parameter character containing TYPE NAME
  !!
  elemental function type(self)
    class(squareCylinder), intent(in) :: self
    character(nameLen)                :: type

    type = TYPE_NAME

  end function type

  !!
  !! Return string with definition of this surface
  !!
  pure subroutine getDef(self,string)
    class(squareCylinder), intent(in)     :: self
    character(:),allocatable,intent(inout) :: string

    string = printSurfDef(TYPE_NAME, [self % a, self % origin])

  end subroutine getDef

  !!
  !! Returns an axis alligned bouding box of surface -ve halfspace
  !!
  pure subroutine boundingBox(self,origin, halfwidth)
    class(squareCylinder), intent(in)       :: self
    real(defReal), dimension(3),intent(out) :: origin
    real(defReal), dimension(3),intent(out) :: halfwidth

    origin    = self % origin
    halfwidth = self % a

    ! Change ZERO in dummy direction to INFINITY
    where (halfwidth == ZERO)
      halfwidth = INFINITY
    end where

  end subroutine boundingBox

  !!
  !! Override base type function to returns .false.
  !!
  function cannotBeBoundary(self) result(itCant)
    class(squareCylinder), intent(in) :: self
    logical(defBool)                   :: itCant

    itCant = .not.ACCEPTS_BC

  end function cannotBeBoundary

  !!
  !! Calculate the distance to the nearest surface of the box
  !! Requires checking that only real surfaces are intercepted,
  !! i.e., not the extensions of the box plane surfaces
  !!
  elemental subroutine distance(self, dist, idx, r, u)
    class(squareCylinder), intent(in) :: self
    real(defReal), intent(out)          :: dist
    integer(shortInt), intent(out)      :: idx
    type(vector), intent(in)            :: r
    type(vector), intent(in)            :: u
    real(defReal),dimension(3)          :: a_near, a_far
    type(vector)                        :: r_bar
    real(defReal)                       :: t_near, t_far
    real(defReal)                       :: test_near, test_far
    integer(shortInt)                   :: near_ax, far_ax, di
    logical(defBool)                    :: farNeg, nearNeg

    ! Transform to frame centered at the slab and choose nearest and furthest vertex
    ! based on octant of direction vector
    r_bar  = r - self % origin
    a_far  = sign(self % a, u % v)
    a_near = -a_far

    ! Initialise test distances
    t_near = -huge(t_near)
    t_far  =  huge(t_far)

    ! Intersection with 1st direction
    di = self % dirs(1)
    if (u % v(di) /= ZERO) then
      test_near = (a_near(di) - r_bar % v(di)) / u % v(di)
      test_far  = (a_far(di) - r_bar % v(di)) / u % v(di)

      ! Find max(t_near,test_near) & save Axis
      if(t_near < test_near) then
        t_near  = test_near
        near_ax = di
      end if

      ! Find min(t_far,test_far) & save Axis
      if(t_far > test_far) then
        t_far = test_far
        far_ax = di
      end if

    end if

    ! Intersection with 2nd direction
    di = self % dirs(2)
    if (u % v(di) /= ZERO) then
      test_near = (a_near(di) - r_bar % v(di)) / u % v(di)
      test_far  =  (a_far(di) - r_bar % v(di)) / u % v(di)

      ! Find max(t_near,test_near) & save Axis
      if(t_near < test_near) then
        t_near  = test_near
        near_ax = di
      end if

      ! Find min(t_far,test_far) & save Axis
      if(t_far > test_far) then
        t_far = test_far
        far_ax = di
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
    class(squareCylinder), intent(in)  :: self
    type(vector), intent(in)            :: r
    real(defReal),dimension(2)          :: sig_r, r_bar
    integer(shortInt)                   :: maxCom
    type(vector)                        :: normal

    associate( di => self % dirs)
      ! Move to approperiate quadrant. And save sign vector to be able to return
      sig_r = sign( [ONE, ONE],r % v(di))
      r_bar = abs(r % v(di))

      ! Calculate vector between the quadrant corner and current point
      r_bar = r_bar - self % a(di)

      ! Identify surface in the octant as the one with the largest component
      maxCom = maxloc(r_bar ,1)

      ! Translate maxCom to cmponent of 3-D vector
      maxCom = di(maxCom)

      ! Build normal vector
      normal = [ZERO ,ZERO, ZERO]
      normal % v(maxCom) = ONE

      ! Multiply resulting vector by sign vector to return to original octant
      normal % v(di) = normal % v(di) * sig_r

    end associate

  end function normalVector


  !!
  !! Set boundary conditions for an squareCylinder: may only be vacuum
  !!
  subroutine setBoundaryConditions(self, BC)
    class(squareCylinder), intent(inout)       :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here ='setBoundaryConditions ( squareCylinder_class.f90)'

    if(size(BC) < 6) call fatalError(Here,'Wrong size of BC string. Must be at least 6')

    ! Load BC codes
    self % BC = BC(1:6)

    ! Verify BC periodic BC
    ! *** I am not apologising this time! - MAK
    ! *** Technically it should allow mixed BC for the dummy axis.
    ! *** Add it to TODO list!
    if(.not.all( (self % BC([1,3,5] ) == periodic) .eqv. (self % BC([2,4,6]) == periodic))) then
      call fatalError(Here,'Periodic BC need to be applied to both oposite surfaces')

    end if

  end subroutine setBoundaryConditions

  !!
  !! Apply boundary transformation
  !!
  subroutine boundaryTransform(self, r, u)
    class(squareCylinder), intent(in)   :: self
    type(vector), intent(inout)          :: r
    type(vector), intent(inout)          :: u
    real(defReal), dimension(3)          :: pos
    real(defReal), dimension(3)          :: a_bar
    integer(shortInt), dimension(3)      :: Ri
    integer(shortInt)                    :: N_trans, i, j, dir, flag
    character(100),parameter :: Here ='boundaryTransform ( squareCylinder_class.f90)'

    associate (di => self % dirs)
      ! Calculate dimension reduced by surface tolerance
      ! This is needed to apply BC for particles at the surface
      a_bar = self % a - surface_tol

      ! Calculate maximum distance in terms of coordinate transforms N_trans
      ! Number of times we need to apply BCs to bring particle back into geometry with
      ! reflections and translations.
      Ri(di) = ceiling(abs(r % v(di)- self % origin(di))/ a_bar(di))/2
      N_trans = max(Ri(di(1)), Ri(di(2)))

      ! Loop over number of transformation required to fold back to inside geometry
      do i=1,N_trans

        ! Calculate position of the point wrt to orgin
        pos = r % v - self % origin

        ! Loop over cardinal directions and apply BC
        do j =1,2
          dir = di(j)
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
    end associate
  end subroutine boundaryTransform


  !!
  !! Macro to apply BC given direction and flag which specifies position wrt planes
  !!
  subroutine applyBC(self,r,u,flag,dir)
    class(squareCylinder), intent(in) :: self
    type(vector), intent(inout)        :: r
    type(vector), intent(inout)        :: u
    integer(shortInt), intent(in)      :: flag
    integer(shortInt), intent(in)      :: dir
    integer(shortInt)                  :: BC, idx
    real(defReal)                      :: off, disp
    character(100),parameter           :: Here ='applyBC( squareCylinder_class.f90)'

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

end module squareCylinder_class
