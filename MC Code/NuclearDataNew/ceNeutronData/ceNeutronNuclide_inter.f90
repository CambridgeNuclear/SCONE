module ceNeutronNuclide_inter

  use numPrecision
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG

  ! Nuclear Data Handles
  use nuclideHandle_inter,   only : nuclideHandle

  use neutronXsPackages_class, only : neutronMicroXSs

  ! CE Neutron Interfaces
  use ceNeutronDatabase_inter, only : ceNeutronDatabase

  ! Cache
  use ceNeutronCache_mod,    only : nuclideCache

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public ceNeutronNuclide_CptrCast


  !!
  !! An abstract class that represents all CE Neutron Nuclide data
  !!
  !! Exist mainly in order to decouple caching logic from the database implementation
  !! so there is no need to repeat it in every database type. Thus it will be easier to
  !! mantain and optimise.
  !!
  !! Private Members
  !!   nucIdx    -> nucIdx for this nuclide
  !!   data      -> pointer to a database to request update of XSs
  !!   fissile   -> flag that specifies if the nuclide is fissile
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
  !!
  type, public, abstract, extends(nuclideHandle) :: ceNeutronNuclide
    private
    integer(shortInt)                :: nucIdx = 0
    class(ceNeutronDatabase),pointer :: data => null()
    logical(defBool)                 :: fissile =.false.
  contains

    procedure, non_overridable :: getTotalXS
    procedure, non_overridable :: getMicroXSs
    procedure, non_overridable :: set
    procedure, non_overridable :: getNucIdx
    procedure, non_overridable :: isFissile

    ! Procedures for specific implementations
    procedure(invertInelastic),deferred :: invertInelastic
    procedure(xsOf), deferred           :: xsOf

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
    if (nuclideCache(self % getNucIdx()) % E_tot /= E) then
      call self % data % updateTotalXS(E, self % nucIdx, rand)
    end if

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
    if(nuclideCache(self % getNucIdx()) % E_tail /= E) then
      call self % data % updateMicroXSs(E, self % nucIdx, rand)
    end if

    xss = nuclideCache(self % getNucIdx()) % xss

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
  !!   isFissile [in] -> flag indicating whether fission data is present
  !!
  subroutine set(self, nucIdx, database, fissile)
    class(ceNeutronNuclide), intent(inout)                :: self
    integer(shortInt), intent(in),optional                :: nucIdx
    class(ceNeutronDatabase),pointer, optional,intent(in) :: database
    logical(defBool),intent(in), optional                 :: fissile

    if(present(nucIdx))    self % nucIdx  = nucIdx
    if(present(database))  self % data    => database
    if(present(fissile))   self % fissile = fissile

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


end module ceNeutronNuclide_inter
