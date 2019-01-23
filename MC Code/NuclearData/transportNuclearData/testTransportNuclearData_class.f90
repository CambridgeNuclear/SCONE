module testTransportNuclearData_class

  use numPrecision
  use particle_class,             only : particle
  use dictionary_class,           only : dictionary
  use xsMacroSet_class,           only : xsMacroSet_ptr, xsMacroSet
  use transportNuclearData_inter, only : transportNuclearData

  implicit none
  private

  !!
  !! Extreamly simple implementation of transportNuclearData
  !! Always returns a constant, given value of a XS
  !! Can return a pointer to a xsMacroSet filled with values set in "build" subroutine
  !! DOES NOT check for any consistency of XSs.
  !! Ignores matIdx
  !! Should be used for TESTING only
  !!
  type, public,extends(transportNuclearData) :: testTransportNuclearData
    private
    real(defReal)            :: xsVal = 0.0
    type(xsMacroSet),pointer :: macroXS => null()
  contains
    ! Class specific procedures
    procedure :: build

    ! Nuclear Data interface
    procedure :: getName
    procedure :: getIdx
    procedure :: kill
    procedure :: init

    ! Transport Nuclear Data Interface
    procedure :: getTransXS_p
    procedure :: getMajorantXS_p
    procedure :: getTotalMatXS_p
    procedure :: getMatMacroXS_p

    procedure :: isFissileMat
    procedure :: initFissionSite
    procedure :: setActiveMaterials
  end type testTransportNuclearData

contains

  !!
  !! Build instance of testTransportNuclearData
  !!
  subroutine build(self, xsVal, scatterXS, captureXS, fissionXS, nuFissionXS)
    class(testTransportNuclearData), intent(inout) :: self
    real(defReal), intent(in)                      :: xsVal
    real(defReal), intent(in),optional             :: scatterXS
    real(defReal), intent(in),optional             :: captureXS
    real(defReal), intent(in),optional             :: fissionXS
    real(defReal), intent(in),optional             :: nuFissionXS

    self % xsVal = xsVal

    allocate(self % macroXS)
    self % macroXS % totalXS = xsVal

    ! Scattering
    if(present(scatterXS)) then
      self % macroXS % scatterXS = scatterXS
    else
      self % macroXS % scatterXS = xsVal
    end if

    ! Capture
    if(present(captureXS)) then
      self % macroXS % captureXS = captureXS
    else
      self % macroXS % captureXS = xsVal
    end if

    ! Fission
    if(present(fissionXS)) then
      self % macroXS % fissionXS = fissionXS
    else
      self % macroXS % fissionXS = xsVal
    end if

    ! nu*Fission
    if(present(nuFissionXS)) then
      self % macroXS % nuFissionXS = nuFissionXS
    else
      self % macroXS % nuFissionXS = xsVal
    end if


  end subroutine build

  !!
  !! Initialise testTransportNuclearData
  !! from material dictionary
  !! and an array of material names, that specifies material indices
  !!
  subroutine init(self, dict, matNames)
    class(testTransportNuclearData), intent(inout) :: self
    class(dictionary), intent(in)                  :: dict
    character(nameLen), dimension(:), intent(in)   :: matNames
    real(defReal)                                  :: xsVal

    call dict % get(xsVal,'xsValue')
    call self % build(xsVal)

  end subroutine init

  !!
  !! Returns xsVal always
  !!
  function getTransXS_p(self,p,matIdx) result (xs)
    class(testTransportNuclearData), intent(inout) :: self
    class(particle), intent(in)                    :: p
    integer(shortInt), intent(in)                  :: matIdx
    real(defReal)                                  :: xs

    xs = self % xsVal

  end function getTransXS_p

  !!
  !! Returns xsVal always
  !!
  function getMajorantXS_p(self,p) result(xs)
    class(testTransportNuclearData), intent(inout) :: self
    class(particle), intent(in)                :: p
    real(defReal)                              :: xs

    xs = self % xsVal

  end function getMajorantXS_p

  !!
  !! Get total XS of a given material
  !!
  function getTotalMatXS_p(self,p,matIdx) result (xs)
    class(testTransportNuclearData), intent(inout) :: self
    class(particle), intent(in)                :: p
    integer(shortInt), intent(in)              :: matIdx
    real(defReal)                              :: xs

    xs = self % xsVal

  end function getTotalMatXS_p

  !!
  !! Point to local value package
  !!
  subroutine getMatMacroXS_p(self,macroXS,p,matIdx)
    class(testTransportNuclearData), intent(inout)  :: self
    type(xsMacroSet_ptr),intent(inout)              :: macroXS
    class(particle), intent(in)                     :: p
    integer(shortInt),intent(in)                    :: matIdx

    macroXS = self % macroXS

  end subroutine getMatMacroXS_p

  !!
  !! Return .false.
  !!
  function isFissileMat(self,matIdx) result(isIt)
    class(testTransportNuclearData), intent(in) :: self
    integer(shortInt), intent(in)           :: matIdx
    logical(defBool)                        :: isIt

    isIt = .false.

  end function isFissileMat

  !!
  !! Does nothing
  !!
  subroutine initFissionSite(self,p,r)
    class(testTransportNuclearData), intent(in) :: self
    class(particle), intent(inout)          :: p
    real(defReal),dimension(3), intent(in)  :: r
  end subroutine initFissionSite

  !!
  !! Does nothing
  !!
  subroutine setActiveMaterials(self,matIdxList)
    class(testTransportNuclearData), intent(inout) :: self
    integer(shortInt),dimension(:), intent(in) :: matIdxList
  end subroutine

  !!
  !! Return empty character
  !!
  function getName(self, matIdx) result(matName)
    class(testTransportNuclearData), intent(in) :: self
    integer(shortInt), intent(in)               :: matIdx
    character(nameLen)                          :: matName

    matName = ''

  end function getName

  !!
  !! Always return 0
  !!
  function getIdx(self,matName) result(matIdx)
    class(testTransportNuclearData), intent(in) :: self
    character(*),intent(in)                     :: matName
    integer(shortInt)                           :: matIdx

    matIdx = 0

  end function getIdx

  !!
  !! Do nothing
  !!
  elemental subroutine kill(self)
    class(testTransportNuclearData), intent(inout) :: self

    deallocate(self % macroXS)

  end subroutine kill

end module testTransportNuclearData_class
