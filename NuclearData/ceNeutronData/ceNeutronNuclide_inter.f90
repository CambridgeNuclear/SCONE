module ceNeutronNuclide_inter

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use RNG_class,         only : RNG

  ! Nuclear Data Handles
  use nuclideHandle_inter,     only : nuclideHandle
  use neutronXsPackages_class, only : neutronMicroXSs

  ! CE Neutron Interfaces
  use ceNeutronDatabase_inter, only : ceNeutronDatabase

  ! Cache
  use ceNeutronCache_mod,      only : nuclideCache

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public ceNeutronNuclide_CptrCast

  !!
  !! Procedures to extend in subclasses
  !!
  public :: kill

  !!
  !! An abstract class that represents all CE Neutron Nuclide data
  !!
  !! Exist mainly in order to decouple caching logic from the database implementation
  !! so there is no need to repeat it in every database type. Thus it will be easier to
  !! mantain and optimise.
  !!
  !! Private Members
  !!   nucIdx  -> nucIdx for this nuclide
  !!   data    -> pointer to a database to request update of XSs
  !!   fissile -> flag that specifies if the nuclide is fissile
  !!   DBRC    -> Doppler Broadening Rejection Correction flag
  !!
  !! Interface:
  !!   nuclideHandle Interface
  !!   getTotalXS      -> Returns total XS for the nuclide
  !!   getMainXSs      -> Returns a XS package with main neutron XSs
  !!   set             -> set nucIdx, data pointer and isFissile flag
  !!   getNucIdx       -> get nuclide index of the nuclide
  !!   isFissile       -> Return .true. if nuclide can fission
  !!   invertInelastic -> Selects type of inelastic neutron scattering
  !!   xsOf            -> Returns microscopic XS given MT number
  !!   elScatteringXS  -> Returns elastic scattering XS for the nuclide
  !!   hasDBRC         -> Returns the value of the hasDBRC flag
  !!
  type, public, abstract, extends(nuclideHandle) :: ceNeutronNuclide
    private
    integer(shortInt)                :: nucIdx  =  0
    class(ceNeutronDatabase),pointer :: data    => null()
    logical(defBool)                 :: fissile = .false.
    real(defReal)                    :: mass    =  ZERO
    real(defReal)                    :: kT      =  ZERO

    ! DBRC nuclide flag
    logical(defBool)                 :: DBRC = .false.

  contains

    procedure, non_overridable :: getTotalXS
    procedure, non_overridable :: getMicroXSs
    procedure, non_overridable :: set
    procedure, non_overridable :: getNucIdx
    procedure, non_overridable :: isFissile
    procedure, non_overridable :: hasDBRC
    procedure                  :: getMass
    procedure                  :: getkT
    procedure                  :: kill

    ! Procedures for specific implementations
    procedure(invertInelastic),deferred :: invertInelastic
    procedure(xsOf), deferred           :: xsOf
    procedure(elScatteringXS), deferred :: elScatteringXS
    procedure(needsSabEl), deferred     :: needsSabEl

  end type ceNeutronNuclide

  abstract interface

    !!
    !! Invert PDF of inelastic stattering
    !!
    !! Sample outgoing reaction channel provided the reaction is lumped in
    !! inelasticScattering on the XS package.
    !!
    !! Args:
    !!   E [in]       -> Energy for the inversion [MeV]
    !!   rand [inout] -> Random Number Generator
    !!
    !! Result:
    !!   MT number of the sampled reaction
    !!
    !! Errors:
    !!   FatalError is energy E is out-of-bounds for avalible data.
    !!
    function invertInelastic(self, E, rand) result(MT)
      import :: ceNeutronNuclide, RNG, shortInt, defReal
      class(ceNeutronNuclide), intent(in) :: self
      real(defReal), intent(in)           :: E
      class(RNG), intent(inout)           :: rand
      integer(shortInt)                   :: MT
    end function invertInelastic

    !!
    !! Return Cross-Section of reaction MT at energy E
    !!
    !! Args:
    !!   MT [in] -> MT number of the required reaction
    !!   E [in]  -> Energy for the cross-section [MeV]
    !!
    !! Result:
    !!   Microscopic reaction cross-section [barn]
    !!
    !! Error:
    !!   fatalError if Energy is out-of-bounds for available data
    !!   fatalError if reaction under MT number is not present
    !!
    function xsOf(self, MT, E) result(xs)
      import :: ceNeutronNuclide, shortInt, defReal
      class(ceNeutronNuclide), intent(in) :: self
      integer(shortInt), intent(in)       :: MT
      real(defReal), intent(in)           :: E
      real(defReal)                       :: xs

    end function xsOf

    !!
    !! Return elastic scattering XS for the nuclide
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!
    !! Result:
    !!   Elastic scattering nuclide microscopic cross-section [barn]
    !!
    !! Errors:
    !!   fatalError if E is out-of-bounds of the present data
    !!   Invalid idx beyond array bounds -> undefined behaviour
    !!   Invalid f (outside [0;1]) -> incorrect value of XS
    !!
    function elScatteringXS(self, E) result(xs)
      import :: ceNeutronNuclide, defReal
      class(ceNeutronNuclide), intent(in) :: self
      real(defReal), intent(in)           :: E
      real(defReal)                       :: xs

    end function elScatteringXS

    !!
    !! Function that checks whether this nuclide at the provided energy should
    !! have S(a,b) elastic scattering data or not
    !!
    !! Args:
    !!   E [in] -> incident neutron energy
    !!
    !! Result:
    !!    True or false
    !!
    elemental function needsSabEl(self, E) result(doesIt)
      import :: ceNeutronNuclide, defReal, defBool
      class(ceNeutronNuclide), intent(in)  :: self
      real(defReal), intent(in)            :: E
      logical(defBool)                     :: doesIt

    end function needsSabEl

  end interface

contains

  !!
  !! Return Total XS for the nuclide
  !!
  !! Args:
  !!   E [in]       -> required energy [MeV]
  !!   rand [inout] -> random number generator
  !!
  !! Result:
  !!   Total nuclide microscopic cross-section [barn]
  !!
  !! Errors:
  !!   fatalError if E is out-of-bounds of the present data
  !!
  function getTotalXS(self, E, rand) result(xs)
    class(ceNeutronNuclide), intent(in) :: self
    real(defReal), intent(in)           :: E
    class(RNG), intent(inout)           :: rand
    real(defReal)                       :: xs

    ! Check Cache and update if needed
    if (nuclideCache(self % nucIdx) % E_tot /= E) then
      call self % data % updateTotalNucXS(E, self % nucIdx, rand)
    end if

    xs = nuclideCache(self % nucIdx) % xss % total

  end function getTotalXS

  !!
  !! Return Microscopic XSs for the nuclide
  !!
  !! Args:
  !!   xss [out]    -> Cross section package to store the data
  !!   E [in]       -> Requested energy [MeV]
  !!   rand [inout] -> Random Number Generator
  !!
  !! Errors:
  !!   fatalError if E is out-of-bounds for the stored data
  !!
  subroutine getMicroXSs(self, xss, E, rand)
    class(ceNeutronNuclide), intent(in) :: self
    type(neutronMicroXSs), intent(out)  :: xss
    real(defReal), intent(in)           :: E
    class(RNG), intent(inout)           :: rand

    ! Check Cache and update if needed
    if (nuclideCache(self % nucIdx) % E_tail /= E) then
      call self % data % updateMicroXSs(E, self % nucIdx, rand)
    end if

    xss = nuclideCache(self % nucIdx) % xss

  end subroutine getMicroXSs

  !!
  !! Set nucIdx and pointer to a database
  !!
  !! All arguments are optional. Use with keyword association e.g.
  !!   call nuc % set( isFissile=.false., nucIdx = 7)
  !!
  !! Use this pprocedure ONLY during build. NEVER during transport.
  !! IT IS NOT THREAD SAFE!
  !!
  !! Args:
  !!   nucIdx [in]    -> nuclide index
  !!   database [in]  -> pointer to a database that updates XSs on the ceNeutronCache
  !!   fissile [in] -> flag indicating whether fission data is present
  !!   mass [in]      -> Mass of nuclide in neutron masses
  !!   kT [in]        -> Temperature [MeV] of the data in the nuclide
  !!
  !! Error:
  !!   fatalError if kT <= 0.0
  !!   fatalError if mass <= 0.0
  !!
  subroutine set(self, nucIdx, database, fissile, mass, kT, dbrc)
    class(ceNeutronNuclide), intent(inout)                  :: self
    integer(shortInt), intent(in),optional                  :: nucIdx
    class(ceNeutronDatabase), pointer, optional, intent(in) :: database
    logical(defBool), intent(in), optional                  :: fissile
    real(defReal), intent(in), optional                     :: mass
    real(defReal), intent(in), optional                     :: kT
    logical(defBool), intent(in), optional                  :: dbrc
    character(100), parameter :: Here = 'set (ceNuetronNuclide_inter.f90)'

    if (present(nucIdx))    self % nucIdx  = nucIdx
    if (present(database))  self % data    => database
    if (present(fissile))   self % fissile = fissile
    if (present(dbrc))      self % DBRC    = dbrc

    if (present(mass)) then
      if (mass <= ZERO) call fatalError(Here,"Mass of nuclide cannot be -ve: "//numToChar(mass))
      self % mass = mass
    end if

    if (present(kT)) then
      if (kT < ZERO) call fatalError(Here, "Temperature of nuclide cannot be -ve: "//numToChar(kT))
      self % kT = kT
    end if


  end subroutine set

  !!
  !! Get nuclide index
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Nuclide index of the nuclide
  !!
  elemental function getNucIdx(self) result(nucIdx)
    class(ceNeutronNuclide), intent(in) :: self
    integer(shortInt)                   :: nucIdx

    nucIdx = self % nucIdx

  end function getNucIdx

  !!
  !! Return .true. if nuclide is fissile
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   .TRUE. if fissile, .FALSE. otherwise
  !!
  !! Errors:
  !!   None
  !!
  pure function isFissile(self) result(isIt)
    class(ceNeutronNuclide), intent(in) :: self
    logical(defBool)                    :: isIt

    isIt = self % fissile

  end function isFissile

  !!
  !! Return .true. if the nuclide needs to use DBRC
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   .TRUE. if DBRC is on, .FALSE. otherwise
  !!
  !! Errors:
  !!   None
  !!
  pure function hasDBRC(self) result(hasIt)
    class(ceNeutronNuclide), intent(in) :: self
    logical(defBool)                    :: hasIt

    hasIt = self % DBRC

  end function hasDBRC

    !!
    !! Return a mass of the nuclide
    !!
    !! See nuclideHandle documentation
    !!
    pure function getMass(self) result(M)
      class(ceNeutronNuclide), intent(in) :: self
      real(defReal)                       :: M

      M = self % mass

    end function getMass

    !!
    !! Return nuclide temperature
    !!
    !! See nuclideHandle documentation
    !!
    pure function getkT(self) result(kT)
      class(ceNeutronNuclide), intent(in) :: self
      real(defReal)                       :: kT

      kT = self % kT

    end function getkT

  !!
  !! Cast nuclideHandle pointer to ceNeutronNuclide pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclideHandle
  !!
  !! Result:
  !!   Null is source is not of ceNeutronNuclide
  !!   Pointer to source if source is ceNuclearDatabase class
  !!
  pure function ceNeutronNuclide_CptrCast(source) result(ptr)
    class(nuclideHandle), pointer, intent(in) :: source
    class(ceNeutronNuclide), pointer          :: ptr

    select type(source)
      class is(ceNeutronNuclide)
        ptr => source

      class default
        ptr => null()
    end select

  end function ceNeutronNuclide_CptrCast

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(ceNeutronNuclide), intent(inout) :: self

    self % nucIdx = 0
    self % data => null()
    self % fissile = .false.

  end subroutine kill

end module ceNeutronNuclide_inter
