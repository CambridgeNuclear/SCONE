module ceNeutronNuclide_inter

  use numPrecision
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG

  ! Nuclear Data Handles
  use nuclideHandle_inter,   only : nuclideHandle

  use neutronXsPackages_class, only : neutronMicroXSs

  ! Cache
  use ceNeutronCache_mod,    only : nuclideCache

  implicit none
  private


  !!
  !! An abstract class that represents all CE Neutron Nuclide data
  !!
  !! Exist mainly in order to decouple caching logic from the database implementation
  !! so there is no need to repeat it in every database type. Thus it will be easier to
  !! mantain and optimise.
  !!
  !! Interface:
  !!   nuclideHandle Interface
  !!   getTotalXS      -> Returns total XS for the nuclide
  !!   getMainXSs      -> Returns a XS package with main neutron XSs
  !!   setNucIdx       -> set nuclide index for the nuclide
  !!   getNucIdx       -> get nuclide index of the nuclide
  !!   isFissile       -> Return .true. if nuclide can fission
  !!   invertInelastic -> Selects type of inelastic neutron scattering
  !!   xsOf            -> Returns microscopic XS given MT number
  !!   updateMainXSs   -> update main XSs set for the nuclide on CE Neutron Cache
  !!   updateTotalXS   -> update total XS for the nuclide on CE Neutron Cache
  !!
  type, public, abstract, extends(nuclideHandle) :: ceNeutronNuclide
    private
    integer(shortInt) :: nucIdx = 0
  contains

    procedure, non_overridable :: getTotalXS
    procedure, non_overridable :: getMicroXSs
    procedure, non_overridable :: setNucIdx
    procedure, non_overridable :: getNucIdx

    ! Procedures for specific implementations
    procedure(isFissile),deferred       :: isFissile
    procedure(invertInelastic),deferred :: invertInelastic
    procedure(xsOf), deferred           :: xsOf
    procedure(updateMicroXSs),deferred  :: updateMicroXSs
    procedure(updateTotalXS),deferred   :: updateTotalXs

  end type ceNeutronNuclide

  abstract interface
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
      import :: ceNeutronNuclide, defBool
      class(ceNeutronNuclide), intent(in) :: self
      logical(defBool)                    :: isIt
    end function isFissile

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
      import :: ceNeutronNuclide, RNG, shortInt
      class(ceNeutronNuclide), intent(in) :: self
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
    !! Ensure that microscopic cross-sections are up-to-date for energy E in ceNeutronCache
    !! ANY CHANGE in ceNeutronChache is POSSIBLE
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   rand [inout] -> Random Number Generator [MeV]
    !!
    !! Errors:
    !!   fatalError if Energy is out-of-bounds for available data
    !!
    subroutine updateMicroXSs(self, E, rand)
      import :: ceNeutronNuclide, defReal, RNG
      class(ceNeutronNuclide), intent(in) :: self
      real(defReal), intent(in)           :: E
      class(RNG), intent(inout)           :: rand
    end subroutine updateMicroXSs

    !!
    !! Ensure that total cross-section are up-to-date for energy E in ceNeutronCache
    !! ANY CHANGE in ceNeutronChache is POSSIBLE
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   rand [inout] -> Random Number Generator [MeV]
    !!
    !! Errors:
    !!   fatalError if Energy is out-of-bounds for available data
    !!
    subroutine updateTotalXS(self, E, rand)
      import :: ceNeutronNuclide, defReal, RNG
      class(ceNeutronNuclide), intent(in) :: self
      real(defReal), intent(in)           :: E
      class(RNG), intent(inout)           :: rand
    end subroutine updateTotalXS


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
    if (nuclideCache(self % getNucIdx()) % E_tot /= E) call self % updateTotalXS(E, rand)

    xs = nuclideCache(self % getNucIdx()) % xss % total

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
    if(nuclideCache(self % getNucIdx()) % E_tail /= E) call self % updateMicroXSs(E, rand)

    xss = nuclideCache(self % getNucIdx()) % xss

  end subroutine getMicroXSs

  !!
  !! Set nucIdx of the nuclide
  !!
  !! Args:
  !!   nucIdx [in] -> nuclide index
  !!
  elemental subroutine setNucIdx(self, nucIdx)
    class(ceNeutronNuclide), intent(inout) :: self
    integer(shortInt), intent(in)          :: nucIdx

    self % nucIdx = nucIdx

  end subroutine setNucIdx

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


end module ceNeutronNuclide_inter
