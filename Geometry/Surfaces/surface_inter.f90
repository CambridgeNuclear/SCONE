module surface_inter

  use numPrecision
  use universalVariables, only : SURF_TOL, VACUUM_BC
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary

  implicit none
  private

  !!
  !! Extandable superclass procedures
  !!
  public :: kill

  !!
  !! Abstract interface for all surfaces
  !!
  !! All surfaces may be represented by some equation F(r) = 0
  !!
  !! By default all surfaces support only VACUUM boundary conditions given
  !! by first entry in BC string. Override relevant procedures in subclasses to change
  !! this behaviour!
  !!
  !! Magnitide of surface tolerance is a property of the surface. By default it is
  !! equal to SURF_TOL parameter.
  !!
  !! Private Members:
  !!   surfId -> Surface ID for this surface
  !!
  !! Interface:
  !!   setId       -> Set surface ID
  !!   id          -> Return surface ID
  !!   setTol      -> Set surface tolerance
  !!   surfTol     -> Get value of surface tolerance
  !!   setBC       -> Load boundary conditions in surface-specific order
  !!   myType      -> Returns a string with surface type name
  !!   init        -> Initialise surface from a dictionary
  !!   boundingBox -> Return definition of axis-aligned bounding box over the surface
  !!   kill        -> Return to unitinitialised state
  !!   halfspace   -> Return halfspace ocupied by a particle
  !!   evaluate    -> Return remainder of the surface equation c = F(r)
  !!   distance    -> Return distance to the surface
  !!   going       -> Determine to which halfspace particle is currently going
  !!   explicitBC  -> Apply explicit BCs
  !!   transformBC -> Apply transform BCs
  !!
  type, public, abstract :: surface
    private
    integer(shortInt)  :: surfId = -1
    real(defReal)      :: surf_tol = SURF_TOL

  contains
    ! Initialisation procedures
    procedure                        :: setId
    procedure                        :: id
    procedure                        :: setTol
    procedure, non_overridable       :: surfTol
    procedure                        :: setBC
    procedure(myType), deferred      :: myType
    procedure(init), deferred        :: init
    procedure(boundingBox), deferred :: boundingBox
    procedure                        :: kill

    ! Runtime procedures
    procedure                     :: halfspace
    procedure(evaluate), deferred :: evaluate
    procedure(distance), deferred :: distance
    procedure(going), deferred    :: going
    procedure                     :: explicitBC
    procedure                     :: transformBC
  end type surface

  abstract interface

    !!
    !! Return surface type name
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Allocatable string with surface type name without leading or trailing blanks
    !!
    pure function myType(self) result(str)
      import :: surface
      class(surface), intent(in) :: self
      character(:), allocatable  :: str
    end function myType

    !!
    !! Initialise surface from a dictionary
    !!
    !! Args:
    !!   dict [in] -> Dictionary with surface definition
    !!
    subroutine init(self, dict)
      import :: surface, dictionary
      class(surface), intent(inout) :: self
      class(dictionary), intent(in) :: dict
    end subroutine init

    !!
    !! Return axis-aligned bounding box
    !!
    !! Note:
    !!   If bounding box is infinate in any axis its smalles value is -INF and highest INF
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Smallest axis-aligned bounding box that contains the entire surface.
    !!   Real array of size 6
    !!   [x_min, y_min, z_min, x_max, y_max, z_max ]
    !!
    pure function boundingBox(self) result(aabb)
      import :: surface, defReal
      class(surface), intent(in)  :: self
      real(defReal), dimension(6) :: aabb
    end function boundingBox

    !!
    !! Evaluate surface expression c = F(r)
    !!
    !! Args:
    !!  r [in] -> Position of the point
    !!
    !! Result:
    !!   Remainder of the surface expression c
    !!
    pure function evaluate(self, r) result(c)
      import :: surface, defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal)                           :: c
    end function evaluate

    !!
    !! Return distance to the surface
    !!
    !! Respects surface transparency.
    !! If position r is such that |F(r)| = |c| < SURF_TOL,
    !! ignore closest crossing in terms of absulute distance |d|
    !! (so for d1 = -eps and d2 = 2*eps d1 would be ignored)
    !!
    !! Args:
    !!  r [in] -> Position of the particle
    !!  u [in] -> DIrection of the particle. Assume norm2(u) = 1.0
    !!
    !! Result:
    !!   +ve distance to the next crossing. If there is no crossing
    !!   in +ve direction, returns INF
    !!
    pure function distance(self, r, u) result(d)
      import :: surface, defReal
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3), intent(in) :: u
      real(defReal)                           :: d
    end function distance

    !!
    !! Returns TRUE if particle is going into +ve halfspace
    !!
    !! Is used to determine to which halfspace particle is going when on a surface.
    !!
    !! Args:
    !!   r [in] -> Position of the particle. As the surface (F(r) ~= 0)
    !!   u [in] -> Direction of the partcle. Assume norm2(u) = 1.0
    !!
    !! Result:
    !!   If particle is moving into +ve halfspace (F(r+eps*u0 > 0)) return true
    !!   If particle is moving into -ve halfspace retunr false
    !!
    pure function going(self, r, u) result(halfspace)
      import :: surface, defReal, defBool
      class(surface), intent(in)              :: self
      real(defReal), dimension(3), intent(in) :: r
      real(defReal), dimension(3), intent(in) :: u
      logical(defBool)                        :: halfspace
    end function going
  end interface

contains

  !!
  !! Set surface ID
  !!
  !! Args:
  !!   id [in] -> Surface ID > 0
  !!
  !! Errors:
  !!   fatalError if ID is <= 0
  !!
  subroutine setId(self, id)
    class(surface), intent(inout) :: self
    integer(shortInt), intent(in) :: id
    character(100), parameter :: Here ='setId (surface_inter.f90)'

    if (id <= 0) then
      call fatalError(Here, 'Trying to set not +ve surface ID: '//numToChar(id)//&
                            ' on surface '//self % myType())
    end if

    self % surfId = id

  end subroutine setId

  !!
  !! Return surface ID
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   +ve Integer surface ID of the surface
  !!
  !! Errors:
  !!   fatalError if surfId has not been set
  !!
  function id(self)
    class(surface), intent(in) :: self
    integer(shortInt)          :: id
    character(100), parameter :: Here = 'id (surface_inter.f90)'

    id = self % surfId
    if (id <= 0) then
      call fatalError(Here, 'Id has not been set for this surface of type '//self % myType())
    end if

  end function id

  !!
  !! Set surface tolerance
  !!
  !! Surface tolerance must be a property of the surface
  !! so the SURF_TOL parameter may represent approximetly a width of
  !! ambigous range
  !!
  !! By default is equal to SURF_TOL parameter
  !!
  !! Args:
  !!   tol [in] -> +ve value of the surface tolaerance
  !!
  !! Errors:
  !!   fatalError is tolaerance is not +ve
  !!
  subroutine setTol(self, tol)
    class(surface), intent(inout) :: self
    real(defReal), intent(in)     :: tol
    character(100), parameter :: Here = 'setTol (surface_inter.f90)'

    if (tol <= ZERO) then
      call fatalError(Here, 'Tolerance for surface: '//self % myType()//' must be +ve is: '//&
                             numToChar(tol))
    end if

    self % surf_tol = tol

  end subroutine setTol

  !!
  !! Return value of the surface tolerance for the surface
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Value of the surface tolerance for the surface
  !!
  pure function surfTol(self) result(tol)
    class(surface), intent(in) :: self
    real(defReal)              :: tol

    tol = self % surf_tol

  end function surfTol

  !!
  !! Set boundary conditions
  !!
  !! All surfaces support single vacuum BC for entire surface by default.
  !! To change this override this procedure in a subclass.
  !!
  !! Allowable BCs:
  !!   BC(1) -> entire surface (only vacuum)
  !!
  !! Args:
  !!   BC [in] -> Integer array with BC flags for diffrent faces. Order is determined
  !!     by a surface subclass. Length is arbitrary
  !!
  !! Errors:
  !!   fatalError BC is not supported
  !!
  subroutine setBC(self, BC)
    class(surface), intent(inout)               :: self
    integer(shortInt), dimension(:), intent(in) :: BC
    character(100),parameter :: Here = 'setBC (surface_inter.f90)'

    if (size(BC) == 0) then
      call fatalError(Here, 'At least one entry in the BC string is required!')
    end if

    if (BC(1) /= VACUUM_BC) then
      call fatalError(Here, self % myType()//' supports only VACUUM BCs. Was given: '//&
                            numToChar(BC(1)))
    end if

    ! Nothing to be done

  end subroutine setBC

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(surface), intent(inout) :: self

    self % surfId = -1
    self % surf_tol = SURF_TOL

  end subroutine kill

  !!
  !! Return true if particle is in +ve halfspace
  !!
  !! For c = F(r) halfspace is:
  !!   +ve  if c > 0
  !!   -ve  if c <= 0
  !!
  !! Use surface tolerance if |c| = |F(r)| < SURF_TOL then halfspace is
  !! determined by direction of the particle (+ve if is is moving into +ve halfspace;
  !! -ve otherwise; determined with `going` procedure)
  !!
  !! Args:
  !!   r [in] -> Particle location
  !!   u [in] -> Particle direction (assume norm2(u) = 1.0 )
  !!
  !! Result:
  !!   True if position is in +ve halfspace. False if it is in -ve
  !!
  pure function halfspace(self, r, u) result(hs)
    class(surface), intent(in)              :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: hs
    real(defReal)                           :: c

    c = self % evaluate(r)
    hs = c > ZERO

    ! Apply surface tolarance
    if (abs(c) < self % surfTol()) then
      hs = self % going(r, u)
    end if

  end function halfspace

  !!
  !! Apply explicit BCs
  !!
  !! Vacuum by default. Override in a subclass to change it!
  !!
  !! Args:
  !!  r [inout] -> Position pre and post BC. Assume that (F(r) ~= 0)
  !!  u [inout] -> Direction pre and post BC. Assume that norm2(u) = 1.0
  !!
  subroutine explicitBC(self, r, u)
    class(surface), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u

    ! Do nothing. Vacuum BC

  end subroutine explicitBC

  !!
  !! Apply co-ordinate transform BC
  !!
  !! Vacuum by default. Override in a subclass to change it!
  !!
  !! Args:
  !!  r [inout] -> Position pre and post BC.
  !!  u [inout] -> Direction pre and post BC. Assume that norm2(u) = 1.0
  !!
  subroutine transformBC(self, r, u)
    class(surface), intent(in)                 :: self
    real(defReal), dimension(3), intent(inout) :: r
    real(defReal), dimension(3), intent(inout) :: u

    ! Do nothing. Vacuum BC

  end subroutine transformBC


end module surface_inter
