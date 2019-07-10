module ceNeutronDatabase_inter

  use numPrecision
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG
  use particle_class,    only : particle, P_NEUTRON, printType
  use charMap_class,     only : charMap

  ! Nuclear Data Handles
  use nuclideHandle_inter,   only : nuclideHandle
  use materialHandle_inter,  only : materialHandle
  use reactionHandle_inter,  only : reactionHandle
  use nuclearDatabase_inter, only : nuclearDatabase

  ! Cache
  use ceNeutronCache_mod,    only : materialCache, majorantCache

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public ceNeutronDatabase_CptrCast

  !!
  !! An abstract base class for all nulcear databases that support CE Neutron
  !!
  !! Its primary goal is to contain CE Neutron caching logic so there is
  !! no need to reproduce it in each database implementation.
  !!
  !! It is also used by material and nuclide handles for CE Neutron data to order an
  !! update of XSs on the cache
  !!
  !! Interface:
  !!   nuclearDatabase Interface
  !!   updateTotalMatXS -> update Total Material XS on CE Neutron Cache
  !!   updateMajorantXS -> update Majorant XS on CE Neutron Cache
  !!   updateMacroXSs   -> update Macroscopic XSs for a selected material
  !!   updateTotalXS    -> update Total XS for a selected nuclide
  !!   updateMicroXSs   -> update Microscopic XSs for a selected nuclide
  !!
  type, public, abstract, extends(nuclearDatabase) :: ceNeutronDatabase

  contains
    ! nuclearDatabase Interface Implementation
    procedure, non_overridable :: getTransMatXS
    procedure, non_overridable :: getTotalMatXS
    procedure, non_overridable :: getMajorantXS

    ! Procedures implemented by a specific CE Neutron Database
    procedure(updateTotalMatXS),deferred :: updateTotalMatXS
    procedure(updateMajorantXS),deferred :: updateMajorantXS
    procedure(updateMacroXSs),deferred   :: updateMacroXSs
    procedure(updateTotalXS),deferred    :: updateTotalNucXS
    procedure(updateMicroXSs),deferred   :: updateMicroXSs
  end type ceNeutronDatabase

  abstract interface
    !!
    !! Make sure that totalXS of material with matIdx is at energy E
    !! in ceNeutronChache
    !!
    !! ANY CHANGE in ceNeutronChache is POSSIBLE
    !!   E.G. All material XSs may be updated to energy E
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateTotalMatXS(self, E, matIdx, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: matIdx
      class(RNG), intent(inout)            :: rand
    end subroutine updateTotalMatXS

    !!
    !! Make sure that the majorant of ALL Active materials is at energy E
    !! in ceNeutronChache
    !!
    !! ANY CHANGE in ceNeutronChache is POSSIBLE
    !!   E.G. All material XSs may be updated to energy E
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   rand [inout] -> random number generator
    !!
    subroutine updateMajorantXS(self, E, rand)
      import :: ceNeutronDatabase, defReal, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      class(RNG), intent(inout)            :: rand
    end subroutine updateMajorantXS

    !!
    !! Make sure that the macroscopic XSs for the material with matIdx are set
    !! to energy E in ceNeutronCache
    !!
    !! ANY CHANGE in ceNeutronChache is POSSIBLE
    !!   E.G. Extra materials may be set to energy E as well
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   matIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateMacroXSs(self, E, matIdx, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: matIdx
      class(RNG), intent(inout)            :: rand
    end subroutine updateMacroXSs

    !!
    !! Make sure that totalXS of nuclide with nucIdx is at energy E
    !! in ceNeutronChache
    !!
    !! ANY CHANGE in ceNeutronChache is POSSIBLE
    !!   E.G. All nuclid XSs may be updated to energy E
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   nucIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateTotalXS(self, E, nucIdx, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: nucIdx
      class(RNG), intent(inout)            :: rand
    end subroutine updateTotalXS

    !!
    !! Make sure that the microscopic XSs for the nuclide with nucIdx are set
    !! to energy E in ceNeutronCache
    !!
    !! ANY CHANGE in ceNeutronChache is POSSIBLE
    !!   E.G. Extra nuclides may be set to energy E as well
    !!
    !! Assume that call to this procedure implies that data is NOT up-to-date
    !!
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   nucIdx [in]  -> material index that needs to be updated
    !!   rand [inout] -> random number generator
    !!
    subroutine updateMicroXSs(self, E, nucIdx, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      integer(shortInt), intent(in)        :: nucIdx
      class(RNG), intent(inout)            :: rand
    end subroutine updateMicroXSs
  end interface

contains

  !!
  !! Return transport XS for material matIdx
  !!
  !! See nuclearDatabase_inter for details!
  !!
  !! Error:
  !!   fatalError if particle is not CE Neutron
  !!
  function getTransMatXS(self, p, matIdx) result(xs)
    class(ceNeutronDatabase), intent(inout) :: self
    class(particle), intent(in)             :: p
    integer(shortInt), intent(in)           :: matIdx
    real(defReal)                           :: xs

    xs = self % getTotalMatXS(p, matIdx)

  end function getTransMatXS

  !!
  !! Return Total XS for matIdx
  !!
  !! See nuclearDatabase_inter for details!
  !!
  !! Error:
  !!   fatalError if particle is not CE Neutron
  !!
  function getTotalMatXS(self, p, matIdx) result(xs)
    class(ceNeutronDatabase), intent(inout) :: self
    class(particle), intent(in)             :: p
    integer(shortInt), intent(in)           :: matIdx
    real(defReal)                           :: xs
    character(100),parameter :: Here = 'getTotalMatXS (ceNeutronDatabase_inter.f90)'

    ! Check dynamic type of the particle
    if(p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Dynamic type of the partcle is not CE Neutron but:'//p % typeToChar())
    end if

    ! Check Cache and update if needed
    if(materialCache(matIdx) % E_tot /= p % E) call self % updateTotalMatXS(p % E, matIdx, p % pRNG)

    ! Return Cross-Section
    xs = materialCache(matIdx) % xss % total

  end function getTotalMatXS

  !!
  !! Return Majorant XS
  !!
  !! See nuclearDatabase_inter for details
  !!
  !! Error:
  !!   fatalError if particle is not CE Neutron
  !!
  function getMajorantXS(self, p) result(xs)
    class(ceNeutronDatabase), intent(inout) :: self
    class(particle), intent(in)             :: p
    real(defReal)                           :: xs
    character(100),parameter :: Here = 'getMajorantXS (ceNeutronDatabase_inter.f90)'

    ! Check dynamic type of the particle
    if(p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Dynamic type of the partcle is not CE Neutron but:'//p % typeToChar())
    end if

    ! Check Cache and update if needed
    if(majorantCache(1) % E /= p % E) call self % updateMajorantXS(p % E, p % pRNG)

    ! Return Cross-Section
    xs = majorantCache(1) % xs

  end function getMajorantXS

  !!
  !! Cast nuclearDatabase pointer to ceNeutronDatabase pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null is source is not of ceNuclearDatabase class
  !!   Target points to source if source is ceNuclearDatabase class
  !!
  pure function ceNeutronDatabase_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(ceNeutronDatabase), pointer           :: ptr

    select type(source)
      class is(ceNeutronDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function ceNeutronDatabase_CptrCast

    
end module ceNeutronDatabase_inter
