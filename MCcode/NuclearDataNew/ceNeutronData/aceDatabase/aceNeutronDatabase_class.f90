module aceNeutronDatabase_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary 
  use RNG_class,         only : RNG
  use charMap_class,     only : charMap 

  ! Nuclear Data Interfaces
  use nuclearData_inter,            only : nuclearData
  use materialHandle_inter,         only : materialHandle
  use nuclideHandle_inter,          only : nuclideHandle 
  use reactionHandle_inter,         only : reactionHandle
  use ceNeutronDatabase_inter,      only : ceNeutronDatabase
  use neutronXSPackages_class,      only : neutronMicroXSs

  ! CE NEUTRON CACHE
  use ceNeutronCache_mod,           only : nuclideCache

  implicit none
  private

  !!
  !! A CE Neutron Database based on ACE file format 
  !!
  !! For now the simplest possible implementation. 
  !!
  !! Public Members: 
  !!
  !!
  !! Interface:
  !!   nuclearData Interface
  !!   ceNeutronDatabase Interface
  !!   
  type, public, extends(ceNeutronDatabase) :: aceNeutronDatabase

  contains
    ! nuclearData Procedures 
    procedure :: kill
    procedure :: matNamesMap 
    procedure :: getMaterial 
    procedure :: getNuclide 
    procedure :: getReaction

    ! ceNeutronDatabase Procedures 
    procedure :: energyBounds 
    procedure :: updateTotalMatXS
    procedure :: updateMajorantXS
    procedure :: updateMacroXSs
    procedure :: updateTotalNucXS
    procedure :: updateMicroXSs

    ! This type procedures 
    procedure :: init      
  end type aceNeutronDatabase



contains 

  !!
  !! Return to uninitialised state 
  !!
  elemental subroutine kill(self) 
    class(aceNeutronDatabase), intent(inout) :: self 
  end subroutine kill

  !!
  !! Return pointer to material names map
  !!
  !! See nuclearData_inter for  more details 
  !!
  function matNamesMap(self) result(map)
    class(aceNeutronDatabase), intent(in) :: self
    type(charMap), pointer                :: map
  
    map => null()

  end function matNamesMap

  !!
  !! Return pointer to material in a database
  !!
  !! See nuclearData_inter for  more details 
  !!
  function getMaterial(self, matIdx) result(mat)
    class(aceNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)         :: matIdx
    class(materialHandle), pointer        :: mat

    mat => null()

  end function getMaterial

  !!
  !! Return pointer to nuclide in a database
  !!
  !! See nuclearData_inter for  more details 
  !!
  function getNuclide(self, nucIdx) result(nuc)
    class(aceNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)         :: nucIdx
    class(nuclideHandle), pointer         :: nuc
  
    nuc => null()

  end function getNuclide

  !!
  !! Return a pointer to a reaction
  !!
  !! See nuclearData_inter for  more details 
  !!
  function getReaction(self, MT, idx) result(reac)
    class(aceNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)         :: MT
    integer(shortInt), intent(in)         :: idx
    class(reactionHandle),pointer         :: reac

    reac => null()

  end function getReaction


  !!
  !! Return energy bounds for data in the database
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine energyBounds(self, E_min, E_max)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(out)            :: E_min
    real(defReal), intent(out)            :: E_max
  
    E_min = ONE
    E_max = ONE

  end subroutine energyBounds

  !!
  !! Make sure that totalXS of material with matIdx is at energy E
  !! in ceNeutronChache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateTotalMatXS(self, E, matIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    class(RNG), intent(inout)             :: rand
  end subroutine updateTotalMatXS

  !!
  !! Make sure that the majorant of ALL Active materials is at energy E
  !! in ceNeutronChache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateMajorantXS(self, E, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    class(RNG), intent(inout)             :: rand
  end subroutine updateMajorantXS

  !!
  !! Make sure that the macroscopic XSs for the material with matIdx are set
  !! to energy E in ceNeutronCache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateMacroXSs(self, E, matIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    class(RNG), intent(inout)             :: rand
  end subroutine updateMacroXSs

  !!
  !! Make sure that totalXS of nuclide with nucIdx is at energy E
  !! in ceNeutronChache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateTotalNucXS(self, E, nucIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: nucIdx
    class(RNG), intent(inout)             :: rand
  end subroutine updateTotalNucXS

  !!
  !! Make sure that the microscopic XSs for the nuclide with nucIdx are set
  !! to energy E in ceNeutronCache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateMicroXSs(self, E, nucIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: nucIdx
    class(RNG), intent(inout)             :: rand
  end subroutine updateMicroXSs

  !!
  !! Initialise Database from dictionary and pointer to self 
  !!
  !!
  !! Args: 
  !!   dict [in] -> Dictionary with the settings
  !!   ptr  [in] -> Pointer to self (instance of the aceNeutronDatabase beeing build) 
  !!     of type nuclearData 
  !!
  !! Errors 
  !!   FatalError is ptr is not assosiated with self 
  !!
  subroutine init(self, dict, ptr) 
    class(aceNeutronDatabase), target, intent(inout) :: self 
    class(dictionary), intent(in)                    :: dict 
    class(nuclearData), pointer, intent(in)          :: ptr 
    character(100), parameter :: Here = 'init (aceNeutronDatabase_class.f90)'

    ! Verify pointer 

 !   if(.not.associated(ptr, self)) then
   !   call fatalError(Here,"Pointer needs to be associated with the self")
 !   end if

    ! Create list of all nuclides 



    ! Build nuclide definitions 

    ! Build Material definitions 



  end subroutine init 








end module aceNeutronDatabase_class
