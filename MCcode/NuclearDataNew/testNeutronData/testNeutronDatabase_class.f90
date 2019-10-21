module testNeutronDatabase_class

  use numPrecision
  use particle_class,        only : particle
  use dictionary_class,      only : dictionary
  use charMap_class,         only : charMap

  ! Nuclear Data Interfaces
  use nuclearDatabase_inter, only : nuclearDatabase
  use materialHandle_inter,  only : materialHandle
  use nuclideHandle_inter,   only : nuclideHandle
  use reactionHandle_inter,  only : reactionHandle

  ! Other Test Neutron Objects
  use testNeutronMaterial_class, only : testNeutronMaterial


  implicit none
  private

  !!
  !! Nuclear Database for Testing
  !! Mimics Neutron Data
  !!
  !! Should be created with "Build" method.
  !! By default sets all XSs to xsVal
  !! Specific XSs can be changed with optional arguments of "Build"
  !!
  !! Public Members:
  !!   xsVal -> Default Value of all XSs
  !!   mat   -> Pointer to testNeutronMaterial
  !!
  !! SAMPLE INPUT DICTIONARY:
  !!
  !! testData {
  !!   xsVal 1.3;
  !! }
  !!
  type, public, extends(nuclearDatabase) :: testNeutronDatabase
    real(defReal)                      :: xsVal = ZERO
    type(testNeutronMaterial), pointer :: mat => null()
    !
  contains
    ! Superclass Interface
    procedure :: init
    procedure :: activate
    procedure :: getTransMatXS
    procedure :: getTotalMatXS
    procedure :: getMajorantXS
    procedure :: matNamesMap
    procedure :: getMaterial
    procedure :: getNuclide
    procedure :: getReaction
    procedure :: kill

    ! Local Procedures
    procedure :: build
  end type testNeutronDatabase

contains

  !!
  !! Build testNeutronDatabase from the individual XSs Values
  !!
  !! Args:
  !!   xsVal [in]       -> Default Value of all XSs
  !!   eScatterXS [in]  -> Optional. Value of Elastic Scatter XS
  !!   ieScatterXS [in] -> Oprional. Value of Inelastic Scatter XS
  !!   captureXS [in]   -> Optional. Value of Capture XS
  !!   fissionXS [in]   -> Optional. Value of Fission XS
  !!   nuFissionXS [in] -> Optional Value of nuFission
  !!
  !! Errors:
  !!   None
  !!
  subroutine build(self, xsVal, eScatterXS, ieScatterXS ,captureXS, fissionXS, nuFissionXS)
    class(testNeutroNDatabase), intent(inout) :: self
    real(defReal), intent(in)                 :: xsVal
    real(defReal), intent(in),optional        :: eScatterXS
    real(defReal), intent(in),optional        :: ieScatterXS
    real(defReal), intent(in),optional        :: captureXS
    real(defReal), intent(in),optional        :: fissionXS
    real(defReal), intent(in),optional        :: nuFissionXS

    self % xsVal = xsVal

    if(associated(self % mat)) deallocate(self % mat)
    allocate(self % mat)

    ! Load Total XS
    self % mat % xss % total = xsVal

    ! Elastic Scattering
    if(present(eScatterXS)) then
      self % mat % xss % elasticScatter = eScatterXS
    else
      self % mat % xss % elasticScatter = xsVal
    end if

    ! Inelastic Scattering
    if(present(ieScatterXS)) then
      self % mat % xss % inelasticScatter = ieScatterXS
    else
      self % mat % xss % inelasticScatter = xsVal
    end if

    ! Capture
    if(present(captureXS)) then
      self % mat % xss % capture = captureXS
    else
      self % mat % xss % capture = xsVal
    end if

    ! Fission
    if(present(fissionXS)) then
      self % mat % xss % fission = fissionXS
    else
      self % mat % xss % fission = xsVal
    end if

    ! nu*Fission
    if(present(nuFissionXS)) then
      self % mat % xss % nuFission = nuFissionXS
    else
      self % mat % xss % nuFission = xsVal
    end if

  end subroutine build


  !!
  !! Initialise Database from dictionary and pointer to self
  !!
  !! See nuclearDatabase_inter for details
  !!
  subroutine init(self, dict, ptr, silent)
    class(testNeutronDatabase), target, intent(inout) :: self
    class(dictionary), intent(in)                     :: dict
    class(nuclearDatabase), pointer, intent(in)       :: ptr
    logical(defBool), optional, intent(in)            :: silent

    ! Allocate Material
    if(associated(self % mat)) deallocate(self % mat)
    allocate(self % mat)

    ! Read xsVal
    call dict % get(self % xsVal, 'xsVal')

  end subroutine

  !!
  !! Activate this nuclearDatabase
  !!
  !! See nuclearDatabase_inter for details
  !!
  subroutine activate(self, activeMat)
    class(testNeutronDatabase), intent(inout)   :: self
    integer(shortInt), dimension(:), intent(in) :: activeMat

    ! Do nothing

  end subroutine activate

  !!
  !! Return value of Material Transport XS for a particle
  !!
  !! See nuclearDatabase_inter for details
  !!
  function getTransMatXS(self, p, matIdx) result(xs)
    class(testNeutronDatabase), intent(inout) :: self
    class(particle), intent(in)               :: p
    integer(shortInt), intent(in)             :: matIdx
    real(defReal)                             :: xs

    xs = self % xsVal

  end function getTransMatXS

  !!
  !! Return value of Material Total XS for a particle
  !!
  !! See nuclearDatabase_inter for details
  !!
  function getTotalMatXS(self, p, matIdx) result(xs)
    class(testNeutronDatabase), intent(inout) :: self
    class(particle), intent(in)               :: p
    integer(shortInt), intent(in)             :: matIdx
    real(defReal)                             :: xs

    xs = self % xsVal

  end function getTotalMatXS

  !!
  !! Return value of Majorant XS for a particle
  !!
  !! See nuclearDatabase_inter for details
  !!
  function getMajorantXS(self, p) result(xs)
    class(testNeutronDatabase), intent(inout) :: self
    class(particle), intent(in)               :: p
    real(defReal)                             :: xs

    xs = self % xsVal

  end function getMajorantXS

  !!
  !! Return pointer to material names map
  !!
  !! See nuclearDatabase_inter for details
  !!
  function matNamesMap(self) result(map)
    class(testNeutronDatabase), intent(in) :: self
    type(charMap), pointer                 :: map

    map => null()

  end function matNamesMap

  !!
  !! Return pointer to material in a database
  !!
  !! See nuclearDatabase_inter for details
  !!
  function getMaterial(self, matIdx) result(mat)
    class(testNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)          :: matIdx
    class(materialHandle), pointer         :: mat

    mat => self % mat

  end function getMaterial

  !!
  !! Return pointer to nuclide in a database
  !!
  !! See nuclearDatabase_inter for details
  !!
  function getNuclide(self, nucIdx) result(nuc)
    class(testNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)          :: nucIdx
    class(nuclideHandle), pointer          :: nuc

    nuc => null()

  end function getNuclide

  !!
  !! Return a pointer to a reaction
  !!
  !! See nuclearDatabase_inter for details
  !!
  function getReaction(self, MT, idx) result(reac)
    class(testNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)          :: MT
    integer(shortInt), intent(in)          :: idx
    class(reactionHandle),pointer          :: reac

    reac => null()

  end function getReaction

  !!
  !! Return to uninitialised state
  !!
  !! See nuclearDatabase_inter for details
  !!
  elemental subroutine kill(self)
    class(testNeutronDatabase), intent(inout) :: self

    if(associated(self % mat)) deallocate(self % mat)
    self % xsVal = ZERO

  end subroutine kill


end module testNeutronDatabase_class
