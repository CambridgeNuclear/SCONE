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
  !! An abstract base class for all nulcear databases that support CE Neutron
  !!
  !! Its primary goal is to contain CE Neutron caching logic so there is
  !! no need to reproduce it in each database implementation.
  !!
  !! Interface:
  !!   nuclearDatabase Interface
  !!   updateTotalMatXS -> update Total Material XS on CE Neutron Cache
  !!   updateMajorantXS -> update Majorant XS on CE Neutron Cache
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
  end type ceNeutronDatabase

  abstract interface
    !!
    !! Make sure that totalXS of material with matIdx is at energy E
    !! in ceNeutronChache
    !!
    !! ANY CHANGE in ceNeutronChache is POSSIBLE
    !!   E.G. All material XSs may be updated to energy E
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
    !! Args:
    !!   E [in]       -> required energy [MeV]
    !!   rand [inout] -> random number generator
    !!
    subroutine updateMajorantXS(self, E, rand)
      import :: ceNeutronDatabase, defReal, shortInt, RNG
      class(ceNeutronDatabase), intent(in) :: self
      real(defReal), intent(in)            :: E
      class(RNG), intent(inout)            :: rand
    end subroutine updateMajorantXS
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

    
end module ceNeutronDatabase_inter
