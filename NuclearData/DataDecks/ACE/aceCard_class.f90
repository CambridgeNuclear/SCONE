module aceCard_class

  use numPrecision
  use endfConstants
  use dataDeck_inter,    only : dataDeck
  use genericProcedures, only : fatalError, openToRead, isInteger, linFind,&
                                targetNotFound, searchError, numToChar

  implicit none
  private

  !!
  !! Function that reads scalar or rank 1 real that contains an integer (e.g. 2.0)
  !! Converts it to integer; Returns error if the real does not contain integer (e.g. 2.1)
  !!
  interface real2Int
    module procedure real2Int_scalar
    module procedure real2Int_array
  end interface

  integer(shortInt), parameter :: unINIT = -17 ! Value of uninitialised integer in this module

  !!
  !! Private type to store all data related to an MT reaction in the ACE card
  !!
  type,private :: MTreaction
    integer(shortInt) :: MT        = unInit   ! ENDF MT number of reaction
    real(defReal)     :: Q         = -999.0   ! Q-value for reaction
    integer(shortInt) :: TY        = unInit   ! Neutron release for reaction +ve integer. 19 means fission
    logical(defBool)  :: isCapture = .true.   ! For internal use. Indicates capture reactions
    logical(defBool)  :: CMframe   = .false.  ! Is data for reaction given in CM frame
    integer(shortInt) :: IE        = unINIT   ! First energy grid index
    integer(shortInt) :: N_xs      = unINIT   ! Number of consecutive XS data points
    integer(shortInt) :: XSp       = unINIT   ! Location of cross sections on XSS table
    integer(shortInt) :: LOCB      = unINIT   ! Raw ACE LOCB parameter from MCNP manual
    integer(shortInt) :: ANDp      = unINIT   ! Location of angular distribution data in XSS table
    integer(shortInt) :: LOCC      = unINIT   ! Offset of energy distribution data in XSS table
    logical(defBool)  :: isotropic = .false.  ! Indicates isotropic scattering at all energies
    logical(defBool)  :: correl    = .false.  ! Correlated angle-energy scattering
  contains
    procedure :: print => print_MTreaction
  end type MTreaction

  !!
  !! ACE DATA CARD PARSER FOR A SINGLE NUCLIDE (e.g. 92238.03c)
  !!
  !! Basic interaction is by reading enteries from advancing read head:
  !! Procedures exist that set read head to the beggining of relevant data
  !! After head is set entery of interest (e.g. elastic scattering mu PDF data) should be read as
  !! a stream of data. For safety head should be reset (resetHead procedure) after reading
  !! of a data block.
  !! For reading it is necessary to specify whether entery at this position should be an integer
  !! or a real. All data in XSS table is stored as reals but based on ACE specification some
  !! positions are known to contains integers. If head is misaligned and under an integer entry
  !! there is a real value (e.g. 2.1) procedures for "readInt" will return an error.
  !!
  !! For specyfication of ACE data format please refer to MCNP Manual Appendix F.
  !! Manual for MCNP 4 from 18. December 2000 was used to create this parser.
  !! Please note that in Appendix F, that was used obsolete flags for interpolation scheme in
  !! angular data are given (Table F-12). In modern ACE histogram interpolation uses flag 1
  !! (instead of 0), and linear interpolation uses flag 2 (instead of 0). In other words ACE flags
  !! are consistant with standard ENDF interpolation flags.
  !!
  !! Please note that header of the ACE DATA CARD is stored in public components of type(aceCard)
  !!
  !! Procedures for MT accept elastic scattering which has MT = N_N_elastic (from endfConstants)
  !!
  type, public, extends(dataDeck) :: aceCard
    private

    ! Public Components
    character(nameLen),public :: ZAID = ''   ! 10 character name ZZZAAA.nnC
    real(defReal),public      :: AW   = -ONE ! Atomic weight ratio. Atomic weight divided by the neutron mass
    real(defReal),public      :: TZ   = -ONE ! Temperature at which data were processed [MeV]
    character(10),public      :: HD   = ' '  ! 10 character date when data wre processed
    character(70),public      :: HK   = ' '  ! 70 character comment
    character(10),public      :: HM   = ' '  ! 10 character MAT indentifier

    ! Private Components
    integer(shortInt)                          :: head = 0    ! Current read location on XSS
    type(MTreaction),dimension(:),allocatable  :: MTdata

    ! Fission related data
    logical(defBool)  :: isFiss       = .false. ! Flag is true if JXS(2) /= 0
    logical(defBool)  :: hasFISBlock  = .false. ! Flag is true if JXS(21) /= 0
    integer(shortInt) :: fissIE       = unINIT  ! First energy grid index for fission XSs
    integer(shortInt) :: fissNE       = unINIT  ! Number of fission XS points
    integer(shortInt) :: fissXSp      = unINIT  ! Location of first fission xs point on XSS table
    integer(shortInt) :: promptNUp    = unINIT  ! Location of prompt NU data
    integer(shortInt) :: delayNUp     = unINIT  ! Location of delay NU data
    integer(shortInt) :: totalNUp     = unINIT  ! Location of total NU data

    ! RAW ACE TABLES *** PUBLIC in DEBUG * WILL BE PRIVATE
    integer(shortInt),dimension(16)        :: NXS
    integer(shortInt),dimension(32)        :: JXS
    real(defReal),dimension(:),allocatable :: XSS

  contains
    ! Superclass procedures
    procedure :: myType

    ! XSs directly from blocks
    procedure, non_overridable :: gridSize
    procedure, non_overridable :: ESZ_XS
    procedure :: ESZblock            ! Returns ESZ block XS (see detailed description) !*LEGACY

    ! Procedures to enquire about MT reaction and set read head to angle or energy data
    !
    procedure :: getMTs              ! Get array of MT numbers present
    procedure :: getScatterMTs
    procedure :: getCaptureMTs
    procedure :: getFissionMTs
    procedure :: numMT               ! Get number of MT reactions present
    procedure :: numMTscatter        ! Number of scattering MT reactions
    procedure :: firstIdxMT          ! Get first MT xs point index
    procedure :: numXsPointsMT       ! Get number of MT XS points
    procedure :: xsMT                ! Get MT XS array
    procedure :: neutronReleaseMT    ! Return TY for MT
    procedure :: isCaptureMT         ! returns .true. if MT is capture reaction
    procedure :: QforMT              ! Return Q for MT
    procedure :: isCMframe           ! Return Ref frame for MT
    procedure :: LOCBforMT           ! Return LOCB for MT
    procedure :: setToAngleMT        ! Set head to beggining of angular for MT reaction
    procedure :: setToEnergyMT       ! Set head to energy for MT
    procedure :: LOCCforMT           ! Get Offset of energy law for MT

    ! Procedures releted to Elastic Scattering
    procedure :: LOCBforEscatter     ! Return LOCB for Elastic Scattering
    procedure :: setToAngleEscatter  ! Set head to beginning of angle data for Elastic Scatter

    ! Procedures releated to Fission data
    procedure :: firstIdxFiss        ! Get first fission Idx
    procedure :: numXsPointsFiss     ! Get number of fission XS points
    procedure :: xsFiss              ! Get fission XS array
    procedure :: isFissile           ! Returns true if nuclide is fissile
    procedure :: hasFIS              ! Returns true if card has FIS block
    procedure :: precursorGroups     ! Return number of precursor groups

    procedure :: hasNuPrompt         ! Is prompt NU present
    procedure :: hasNuTotal          ! Is total NU present
    procedure :: hasNuDelayed        ! Is delayed NU present

    procedure :: setToNuPrompt       ! Set head to prompt NU
    procedure :: setToNuTotal        ! Set head to total NU
    procedure :: setToNuDelayed      ! Set head to delayed Nu
    procedure :: setToPrecursors     ! Set head to the precursors data
    procedure :: setToPrecursorEnergy ! Set Head to energy distribution for a given precursor group
    procedure :: LOCCforPrecursor    ! Return location of delayed fission spectrum relative to JXS(27)

    ! Procedures related to reading under head and head status
    procedure :: readInt             ! Read int under head and advance
    procedure :: readReal            ! Read real under head and advance
    procedure :: readIntArray        ! Read N ints beginning under head and advance by N
    procedure :: readRealArray       ! Read N ints beginning under head and advance by N

    procedure :: readIntNotAdvance   ! Read int under head and NOT ADVANCE
    procedure :: readRealNotAdvance  ! Read real under head and NOT ADVANCE

    procedure :: advanceHead         ! Move head by an integer +ve forward -ve backward
    procedure :: resetHead           ! Move head back to position 0

    procedure :: setToAnglePdf       ! Sets head to single energy mu pdf. Adress relative to JXS(9)
    procedure :: setToEnergyLaw      ! Sets head to single energy law. Adress relative to JXS(11)

    procedure :: getRootAddress      ! Get Root adress do diffrent blocks
    procedure :: setRelativeTo       ! Sets to position given by root and offset (root + offset -1)

    ! Procedures related to probability tables
    procedure :: hasProbTab          ! Returns true if the nuclide has URR probability tables
    procedure :: setToProbTab        ! Set head to probability table data

    ! Initialisation and display procedures
    procedure :: readFromFile
    procedure :: print => print_aceCard

    ! Private procedures
    procedure,private :: setMTdata
    procedure,private :: setFissionData
    procedure,private :: getMTidx

  end type aceCard

contains
  !!
  !! Return String with type name
  !!
  !! See dataDeck_inter for details
  !!
  pure function myType(self) result(type)
    class(aceCard), intent(in) :: self
    character(:),allocatable :: type

    type = 'aceCard_class'

  end function myType

  !!
  !! Determines if the scalar real is integer. If it is converts it to shortInt.
  !! If it isn't returns error
  !! NOTE: The most important function in this module. Protects against errors due to offset
  !!       of a pointer at XSS table
  !!
  function real2Int_scalar(r,Where) result (i)
    real(defReal), intent(in)  :: r
    character(*), intent(in)   :: Where
    integer(shortInt)          :: i

    if (.not.isInteger(r)) call fatalError(Where,' Tried to read real into integer variable')
    i = int(r,shortInt)

  end function real2Int_scalar

  !!
  !! Determines if the scalar real is integer. If it is converts it to shortInt.
  !! If it isn't returns error
  !! NOTE: The most important function in this module. Protects against errors due to offset
  !!       of a pointer at XSS table
  !!
  function real2Int_array(r,Where) result(i)
    real(defReal), dimension(:), intent(in) :: r
    character(*), intent(in)                :: Where
    integer(shortInt),dimension(size(r))    :: i

    if (all(isInteger(r))) then
      i = int(r,shortInt)
    else
      call fatalError(Where,'Tried to read reals into integer array')
    end if

  end function real2Int_array

  !!
  !! Returns size of the energy grid
  !!
  pure function gridSize(self)
    class(aceCard), intent(in) :: self
    integer(shortInt)          :: gridSize

    gridSize = self % NXS(3)

  end function gridSize

  !!
  !! Returns data from ESZ block
  !!
  !! Hopefully will not blow up the stack
  !!
  !! Following request can be delivered (CASE SENSITIVE):
  !!  'energyGrid'   -> energy grid
  !!  'totalXS'      -> total cross-section
  !!  'absorptionXS' -> total absorbtion cross-section (without fission)
  !!  'elasticXS'    -> elastic scattering cross-section
  !!  'heatingNumber'-> average heating number
  !!
  !! Args:
  !!   request [in] -> character with the XS request
  !!
  !! Error:
  !!   fatalError if request is unrecognised
  !!
  function ESZ_XS(self, request) result(xs)
    class(aceCard), intent(in)             :: self
    character(*), intent(in)               :: request
    real(defReal),dimension(self % NXS(3)) :: xs
    integer(shortInt)                      :: ptr, N
    character(100), parameter :: Here = 'ESZ_XS (aceCard_class.f90)'

    ! Obtain size of the XS grid
    N = self % NXS(3)

    select case(request)
      case('energyGrid')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1)

      case('totalXS')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + N

      case('absorptionXS')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + 2*N

      case('elasticXS')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + 3*N

      case('heatingNumber')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + 4*N

      case default
        call fatalError(Here,'Unrecognised request string: |' // request //'|')
        ptr = 0
      end select

      ! Return data
      xs = self % XSS(ptr : ptr + N-1)

  end function ESZ_XS

  !!
  !! Returns data from ESZ block
  !! request string specifies which array is returned:
  !!  'energyGrid'   -> energy grid
  !!  'totalXS'      -> total cross-section
  !!  'absorptionXS' -> total absorbtion cross-section (without fission)
  !!  'elasticXS'    -> elastic scattering cross-section
  !!  'heatingNumber'-> average heating number
  !!
  subroutine ESZblock(self,data,request)
    class(aceCard),intent(in)                          :: self
    real(defReal),dimension(:),allocatable,intent(out) :: data
    character(*), intent(in)                           :: request
    integer(shortInt)                                  :: ptr, N
    character(100),parameter  :: Here ='ESZblock (aceCard_class.f90)'

    ! Obtain size of the XS grid
    N = self % NXS(3)

    select case(request)
      case('energyGrid')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1)

      case('totalXS')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + N

      case('absorptionXS')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + 2*N

      case('elasticXS')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + 3*N

      case('heatingNumber')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + 4*N

      case default
        call fatalError(Here,'Unrecognised request string: |' // request //'|')
        ptr = 0
      end select

      ! Return data
      data = self % XSS(ptr : ptr + N-1)

  end subroutine ESZblock

  !!
  !! Returns array of present MT reactions
  !!
  function getMTs(self) result(MT)
    class(aceCard), intent(in)                        :: self
    integer(shortInt),dimension(size(self % MTdata))  :: MT

    MT = self % MTdata(:) % MT

  end function

  !!
  !! Returns array of scattering MTs. (Secondary particles but not fission)
  !!
  function getScatterMTs(self) result(MT)
    class(aceCard), intent(in)                      :: self
    logical(defBool),dimension(size(self % MTdata)) :: mask
    integer(shortInt),dimension(:),allocatable          :: MT

    mask = (self % MTdata % TY /= 19) .and. (.not.self % MTdata % isCapture)
    MT = pack(self % MTdata % MT, mask)

  end function getScatterMTs

  !!
  !! Returns array of capture MTs
  !!
  function getCaptureMTs(self) result(MT)
    class(aceCard), intent(in)                      :: self
    logical(defBool),dimension(size(self % MTdata)) :: mask
    integer(shortInt),dimension(:),allocatable          :: MT

    mask = self % MTdata % isCapture
    MT = pack(self % MTdata % MT, mask)

  end function getCaptureMTs

  !!
  !! Returns array of fission MTs
  !!
  function getFissionMTs(self) result(MT)
    class(aceCard), intent(in)                      :: self
    logical(defBool),dimension(size(self % MTdata)) :: mask
    integer(shortInt),dimension(:),allocatable          :: MT

    mask = self % MTdata % TY == 19
    MT = pack(self % MTdata % MT, mask)

  end function getFissionMTs

  !!
  !! Returns number of MT reactions in a nuclide
  !!
  function numMT(self) result(N)
    class(aceCard), intent(in) :: self
    integer(shortInt)          :: N

    N = size(self % MTdata)

  end function numMT

  !!
  !! Returns number of scattering MT reactions.
  !!
  function numMTscatter(self) result(N)
    class(aceCard), intent(in) :: self
    integer(shortInt)          :: N

    N = self % NXS(5)

  end function numMTscatter

  !!
  !! For a given MT number returns first index of its XS at the energy grid
  !!
  function firstIdxMT(self, MT) result(IE)
    class(aceCard),intent(in)    :: self
    integer(shortInt),intent(in) :: MT
    integer(shortInt)            :: IE
    integer(shortInt)            :: idx

    ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      IE = 1
      return
    end if

    ! Get index of MT number on MTdata
    idx = self % getMTidx(MT)

    ! Return first index on an energy grid
    IE = self % MTdata(idx) % IE

  end function firstIdxMT

  !!
  !! For a given MT number returns total number of XS data points
  !!
  function numXSPointsMT(self, MT) result(N)
    class(aceCard),intent(in)    :: self
    integer(shortInt),intent(in) :: MT
    integer(shortInt)            :: N
    integer(shortInt)            :: idx

    ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      N = self % NXS(3)
      return
    end if

    ! Get index of MT number on MTdata
    idx = self % getMTidx(MT)

    ! Return number of energy points
    N = self % MTdata(idx) % N_xs

  end function numXSPointsMT

  !!
  !! For a given MT number returns array of its XS data point
  !!
  function xsMT(self, MT) result (xs)
    class(aceCard), intent(in)             :: self
    integer(shortInt), intent(in)          :: MT
    real(defReal),dimension(:),allocatable :: xs
    integer(shortInt)                      :: idx, N_xs, ptr

    ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      call self % ESZblock(xs,'elasticXS')
      return
    end if

    ! Get index of MT number on MTdata
    idx = self % getMTidx(MT)

    ! Get number of energy points
    N_xs = self % MTdata(idx) % N_xs

    ! Get Xs data pointer
    ptr = self % MTdata(idx) % XSp

    ! Get XS data
    xs = self % XSS(ptr : ptr + N_xs-1)

  end function xsMT

  !!
  !! Returns integer neutron release
  !! This is abs(TY) where TY is defined in MCNP manual
  !! TY = 19 indicates fission and NU data to be used
  !! TY > 100 indicates anothe NU like table for release (energy dependant yield)
  !!
  function neutronReleaseMT(self,MT) result(N)
    class(aceCard), intent(in)   :: self
    integer(shortInt),intent(in) :: MT
    integer(shortInt)            :: N
    integer(shortInt)            :: idx

   ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      N = ONE
      return
    end if

    idx = self % getMTidx(MT)

    N = self % MTdata(idx) % TY

  end function neutronReleaseMT

  !!
  !! Returns .true. if the reaction under MT is capture
  !!
  function isCaptureMT(self,MT) result(isIt)
    class(aceCard),intent(in)    :: self
    integer(shortInt),intent(in) :: MT
    logical(defBool)             :: isIt
    integer(shortInt)            :: idx

    ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      isIt = .false.
      return
    end if

    ! Get index of MT number on MTdata
    idx = self % getMTidx(MT)

    ! Check if MT represents capture
    isIt = self % MTdata(idx) % isCapture

  end function isCaptureMT


  !!
  !! Resturns Q-value for reaction MT
  !!
  function QforMT(self,MT) result(Q)
    class(aceCard), intent(in)    :: self
    integer(shortInt), intent(in) :: MT
    real(defReal)                 :: Q
    integer(shortInt)             :: idx

   ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      Q = ZERO
      return
    end if

    idx = self % getMTidx(MT)

    Q = self % MTdata(idx) % Q

  end function QforMT

  !!
  !! Returns .true. is secondary neutron data for a MT reaction
  !! is in Center-of-Mass (CM) frame
  !!
  function isCMframe(self,MT) result (isIt)
    class(aceCard), intent(in)    :: self
    integer(shortInt), intent(in) :: MT
    logical(defBool)              :: isIt
    integer(shortInt)             :: idx

   ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      isIt = .true.
      return
    end if

    idx = self % getMTidx(MT)

    isIt = self % MTdata(idx) % CMframe

  end function isCMframe

  !!
  !! Returns raw LOCB parameter for MT reaction as defined in MCNP Manual
  !! Returns error is reaction is capture
  !!
  function LOCBforMT(self,MT) result(LOCB)
    class(aceCard), intent(in)    :: self
    integer(shortInt), intent(in) :: MT
    integer(shortInt)             :: LOCB
    integer(shortInt)             :: idx
    character(100), parameter :: Here='LOCBforMT (aceCard_class.f90)'

   ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      LOCB = self % LOCBforEscatter()
      return
    end if

    idx = self % getMTidx(MT)

    if(self % MTdata(idx) % isCapture) call fatalError(Here,'MT reaction is capture. No LOCB data.')

    LOCB = self % MTdata(idx) % LOCB

  end function LOCBforMT

  !!
  !! Sets read head to beginning of angular data for reaction MT
  !! Returns error is reaction is capture
  !!
  subroutine setToAngleMT(self,MT)
    class(aceCard), intent(inout) :: self
    integer(shortInt), intent(in) :: MT
    integer(shortInt)             :: idx
    character(100), parameter :: Here='setToAngleMT (aceCard_class.f90)'

   ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      call self % setToAngleEscatter()
      return
    end if

    idx = self % getMTidx(MT)

    if(self % MTdata(idx) % isCapture) call fatalError(Here,'MT reaction is capture. Angle data &
                                                             & does not exist')

    self % head = self % MTdata(idx) % ANDp

  end subroutine setToAngleMT

  !!
  !! Sets read head to beginning of energy data for reaction MT
  !! Returns error is reaction is capture
  !!
  subroutine setToEnergyMT(self,MT)
    class(aceCard), intent(inout) :: self
    integer(shortInt), intent(in) :: MT
    integer(shortInt)             :: idx
    character(100), parameter :: Here='setToEnergyMT (aceCard_class.f90)'

   ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      call fatalError(Here,'Elastic scattering has no energy law')
      return
    end if

    idx = self % getMTidx(MT)

    if(self % MTdata(idx) % isCapture) call fatalError(Here,'MT reaction is capture. Energy data &
                                                             & does not exist')

    self % head = self % MTdata(idx) % LOCC + self % JXS(11) - 1

  end subroutine setToEnergyMT

  !!
  !! Returns LOCC paramter for energy law related to MT reaction
  !!
  !! LOCC is position relative to JXS(11)
  !!
  !! Args:
  !!   MT [in] -> Requested MT number
  !!
  !! Errors:
  !!   fatalError if requested MT is not present
  !!
  function LOCCforMT(self, MT) result(locc)
    class(aceCard),intent(in)      :: self
    integer(shortInt), intent(in)  :: MT
    integer(shortInt)              :: locc
    integer(shortInt)              :: idx
    character(100), parameter :: Here = 'LOCCforMT (aceCard_class.f90)'

       ! Special case for elastic scattering
    if( MT == N_N_elastic) then
      call fatalError(Here,'Elastic scattering has no energy law and LOCC')
      locc = 0 ! Avoid Compiler Warning
      return
    end if

    idx = self % getMTidx(MT)
    locc = self % MTdata(idx) % LOCC

  end function LOCCforMT


  !!
  !! Returns LOCB for elastic scattering
  !!
  function LOCBforEscatter(self) result(LOCB)
    class(aceCard), intent(in) :: self
    integer(shortInt)          :: LOCB
    integer(shortInt)          :: ptr
    character(100), parameter :: Here = 'LOCBforEscatter (aceCard_class.f90)'

    ptr = self % JXS(8)

    LOCB = real2Int( self % XSS(ptr),Here)

  end function LOCBforEscatter

  !!
  !! Sets read head to beginning of angular data for ELASTIC SCATTERING
  !!
  subroutine setToAngleEscatter(self)
    class(aceCard), intent(inout) :: self

    self % head = self % JXS(9)

  end subroutine setToAngleEscatter

  !!
  !! Returns first index of fission XS data on the energy grid
  !! Returns error is nuclide is not fissile
  !!
  function firstIdxFiss(self) result(IE)
    class(aceCard), intent(in) :: self
    integer(shortInt)          :: IE
    character(100), parameter  :: Here='firstIdxFiss (aceCard_class.f90)'

    if(.not. self % isFiss) call fatalError(Here,'Nuclide: ' // self % ZAID //' is not fissile')
    IE = self % fissIE

  end function firstIdxFiss

  !!
  !! Returns number of fission xs data points
  !! Returns error is nuclide is not fissile
  !!
  function numXsPointsFiss(self) result(N)
    class(aceCard), intent(in) :: self
    integer(shortInt)          :: N
    character(100), parameter  :: Here='numXsPointsFiss (aceCard_class.f90)'

    if(.not.self % isFiss) call fatalError(Here,'Nuclide: ' // self % ZAID //' is not fissile')
    N = self % fissNE

  end function numXsPointsFiss

  !!
  !! Returns number of fission xs data points
  !! Returns error is nuclide is not fissile
  !!
  function xsFiss(self) result(xs)
    class(aceCard), intent(in) :: self
    real(defReal),dimension(:),allocatable :: xs
    integer(shortInt)                      :: ptr, N
    character(100), parameter  :: Here='xsFiss(aceCard_class.f90)'

    if(self % isFiss) then
      ! Get number of XS points
      N = self % fissNE

      ! Set pointer to XS data
      ptr = self % fissXSp

      ! Read xs data
      xs = self % XSS(ptr : ptr + N-1)

    else
      call fatalError(Here,'Nuclide: ' // self % ZAID //' is not fissile')

    end if

  end function xsFiss

  !!
  !! Returns .true. if nuclide is fissile
  !!
  function isFissile(self) result(isIt)
    class(aceCard), intent(in) :: self
    logical(defBool)           :: isIt

    isIt = self % isFiss

  end function isFissile

  !!
  !! Return .true. if ACE card has FIS block
  !!
  function hasFIS(self) result(isIt)
    class(aceCard), intent(in) :: self
    logical(defBool)           :: isIt

    isIt = self % hasFISBlock

  end function hasFIS

  !!
  !! Return number of precursor groups
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of precursor groups (NXS(8))
  !!
  function precursorGroups(self) result(N)
    class(aceCard), intent(in) :: self
    integer(shortInt)          :: N

    N = self % NXS(8)

  end function precursorGroups


  !!
  !! Returns .true. if nuclide has NU prompt data
  !! If only one prompt/total table is given it is regardes as total
  !! Returns error is nuclide is not fissile
  !!
  function hasNuPrompt(self) result(doesIt)
    class(aceCard), intent(in) :: self
    logical(defBool)           :: doesIt
    character(100), parameter  :: Here='hasNuPrompt(aceCard_class.f90)'

    if(.not. self % isFiss) call fatalError(Here,'Nuclide: ' // self % ZAID //' is not fissile')
    doesIt = (self % promptNUp /= unINIT)

  end function hasNuPrompt

  !!
  !! Returns .true. if nuclide has NU total data
  !! If only one prompt/total table is given it is regardes as total
  !! Returns error is nuclide is not fissile
  !!
  function hasNuTotal(self) result(doesIt)
    class(aceCard), intent(in) :: self
    logical(defBool)           :: doesIt
    character(100), parameter  :: Here='hasNuPrompt(aceCard_class.f90)'

    if(.not. self % isFiss) call fatalError(Here,'Nuclide: ' // self % ZAID //' is not fissile')
    doesIt = (self % totalNUp /= unINIT)

  end function hasNuTotal

  !!
  !! Returns .true. if nuclide has NU delayed data
  !! Returns error is nuclide is not fissile
  !!
  function hasNuDelayed(self) result(doesIt)
    class(aceCard), intent(in) :: self
    logical(defBool)           :: doesIt
    character(100), parameter  :: Here='hasNuPrompt(aceCard_class.f90)'

    if(.not. self % isFiss) call fatalError(Here,'Nuclide: ' // self % ZAID //' is not fissile')
    doesIt = (self % delayNUp /= unINIT)

  end function hasNuDelayed

  !!
  !! Sets read head to beginning of NU prompt Data
  !! If only one prompt/total table is given it is regardes as total
  !! If there is no prompt data returns an error
  !! Returns error is nuclide is not fissile
  !!
  subroutine setToNuPrompt(self)
    class(aceCard), intent(inout) :: self
    character(100),parameter :: Here ='setToNuPromt (aceCard_class.f90)'

    if(self % hasNuPrompt()) then
      self % head = self % promptNUp

    else
      call fatalError(Here,'Should never happen')

    end if

  end subroutine setToNuPrompt

  !!
  !! Sets read head to beginning of NU total Data
  !! If only one prompt/total table is given it is regardes as total
  !! If there is no total data returns an error
  !! Returns error is nuclide is not fissile
  !!
  subroutine setToNuTotal(self)
    class(aceCard), intent(inout) :: self
    character(100),parameter :: Here ='setToNuTotal (aceCard_class.f90)'

    if(self % hasNuTotal()) then
      self % head = self % totalNUp

    else
      call fatalError(Here,'Should never happen')

    end if

  end subroutine setToNuTotal

  !!
  !! Sets read head to beginning of NU delayed Data
  !! If there is no total data returns an error
  !! Returns error is nuclide is not fissile
  !!
  subroutine setToNuDelayed(self)
    class(aceCard), intent(inout) :: self
    character(100),parameter :: Here ='setToNuDelayed (aceCard_class.f90)'

    if(self % hasNuDelayed()) then
      self % head = self % delayNUp

    else
      call fatalError(Here,'Should never happen')

    end if

  end subroutine setToNuDelayed

  !!
  !! Sets head to beginning of precursor data for delayed fission neutrons
  !!
  !!
  !! Errors:
  !!   fatalError if there is no data for precursors
  !!
  subroutine setToPrecursors(self)
    class(aceCard), intent(inout) :: self
    character(100), parameter :: Here = 'setToPrecursors (aceCard_class.f90)'

    if(self % JXS(25) /= 0) then
      self % head = self % JXS(25)

    else
      call fatalError(Here, 'Missing Fission Data. Cannot Locate Precursor PDF. JXS(25) == 0')

    end if

  end subroutine setToPrecursors

  !!
  !! Set Head of ACE Card to beginning of energy data for emission from
  !! precursor group N
  !!
  !! Args:
  !!   N [in] -> Precursor Group index
  !!
  !! Errors:
  !!   fatalError if N is invalid [-ve, 0, or > NXS(8)]
  !!
  subroutine setToPrecursorEnergy(self, N)
    class(aceCard), intent(inout) :: self
    integer(shortInt), intent(in) :: N
    character(100),parameter :: Here = 'setToPrecursorEnergy (aceCard_class.f90)'

    ! Validate N
    if( N < 1) then
      call fatalError(Here,trim(self % ZAID) //' Precursor group index must be non-negative. &
                           & Was: '//numToChar(N))
    else if( N > self % NXS(8)) then
      call fatalError(Here,trim(self % ZAID) // 'N is to large. Was given' // numToChar(N) // &
                           ' with only ' // numToChar(self % NXS(8)) // ' groups present')
    end if

    ! Set address
    self % head = int(self % XSS(self % JXS(26) + N - 1) + self % JXS(27))

  end subroutine setToPrecursorEnergy

  !!
  !! Return LOCC for a precursor group
  !!
  !! Args:
  !!   N [in] -> Precursor group index
  !!
  !! Result:
  !!   LOCC - loccation of enrgy LAW relative to JXS(27)
  !!
  !! Erros:
  !!   fatalError if N is invalid [-ve, 0 or > NXS(8)]
  !!
  function LOCCforPrecursor(self, N) result(locc)
    class(aceCard), intent(inout) :: self
    integer(shortInt), intent(in) :: N
    integer(shortInt)             :: locc
    character(100), parameter :: Here = 'LOCCforPrecursor (aceCard_class.f90)'

    ! Check invalid N
    if ( N < 1 .or. N > self % NXS(8)) then
      call fatalError(Here, trim(self % ZAID)//' Invalid precursor group: '//numToChar(N))
      locc = 0
    end if

    locc = real2Int(self % XSS(self % JXS(26) + N -1), Here)

  end function LOCCforPrecursor

  !!
  !! Get address to a root for different blocks
  !!
  !! Allowed Blocks (Case Sensitive):
  !!   'angleLaws'
  !!   'energyLawsMT'
  !!   'energyLawsPrecursors'
  !!
  !! Args:
  !!   block [in] -> Character that specifies requested block
  !!
  !! Result:
  !!   Address of the block on the XSS array
  !!
  !! Errors:
  !!   fatalError if the block is unrecognised
  !!
  function getRootAddress(self, block) result(addr)
    class(aceCard), intent(in) :: self
    character(*), intent(in)   :: block
    integer(shortInt)          :: addr
    character(100), parameter :: Here = 'getRootAddress (aceCard_class.f90)'

    select case(block)
      case('angleLaws')
        addr = self % JXS(9)

      case('energyLawsMT')
        addr = self % JXS(11)

      case('energyLawsPrecursors')
        addr = self % JXS(27)

      case default
        call fatalError(Here, 'Unrecognised block: '//block)
        addr = 0
    end select

  end function getRootAddress

  !!
  !! Sets head to position given by Root and and offset
  !!
  !! Head = root + offset - 1
  !!
  !! Args:
  !!   root   [in] -> Integer thet specifies the root
  !!   offset [in] -> Integer that specifies an offset
  !!
  !! Errors:
  !!   fatalError if the head is outside the bounds of XSS
  !!   fatalError if root is outside the bounds of XSS
  !!
  subroutine setRelativeTo(self, root, offset)
    class(aceCard), intent(inout) :: self
    integer(shortInt), intent(in) :: root
    integer(shortInt), intent(in) :: offset
    integer(shortInt)             :: pos
    character(100), parameter :: Here = 'setRelativeTo (aceCard_class.f90)'

    ! Validate root
    if (root < 1 .or. root > self % NXS(1)) then
      call fatalError(Here,trim(self % ZAID)//' root is outside bounds: '//numToChar(root))
    end if

    pos = root + offset - 1

    ! Validate position
    if(pos < 1 .or. pos > self % NXS(1)) then
      call fatalError(Here,trim(self % ZAID)//' position reached outside the bounds')
    end if

    self % head = pos

  end subroutine setRelativeTo

  !!
  !! Read single integer and advance read head
  !! Return error if under head is not integer
  !!
  function readInt(self) result(i)
    class(aceCard), intent(inout) :: self
    integer(shortInt)             :: i
    character(100),parameter :: Here ='readInt (aceCard_class.f90)'

    ! Check head status
    if(self % head <= 0) call fatalError(Here,'Reading with unset (not +ve) head')

    ! Read value
    i = real2Int(self % XSS(self % head),Here)

    ! Advance head
    self % head = self % head + 1

  end function readInt

  !!
  !! Read single real and advance read head
  !!
  function readReal(self) result(r)
    class(aceCard), intent(inout) :: self
    real(defReal)                 :: r
    character(100),parameter :: Here ='readReal (aceCard_class.f90)'

    ! Check head status
    if(self % head <= 0) call fatalError(Here,'Reading with unset (not +ve) head')

    ! Read value
    r = self % XSS(self % head)

    ! Advance head
    self % head = self % head + 1

  end function readReal

  !!
  !! Read array of N integers and advance read head by N
  !! Return error if any content under head is not an integer
  !!
  function readIntArray(self,N) result(i)
    class(aceCard), intent(inout)  :: self
    integer(shortInt), intent(in)  :: N
    integer(shortInt),dimension(N) :: i
    integer(shortInt)              :: ptr
    character(100),parameter :: Here ='readIntArray (aceCard_class.f90)'

    ! Check head status
    if(self % head <= 0) call fatalError(Here,'Reading with unset (not +ve) head')

    ! Read value
    ptr = self % head
    i = real2Int( self % XSS(ptr : ptr + N-1),Here)

    ! Advance head
    self % head = self % head + N

  end function readIntArray

  !!
  !! Read array of N integers and advance read head by N
  !!
  function readRealArray(self,N) result(r)
    class(aceCard), intent(inout)  :: self
    integer(shortInt), intent(in)  :: N
    real(defReal),dimension(N)     :: r
    integer(shortInt)              :: ptr
    character(100),parameter :: Here ='readRealArray (aceCard_class.f90)'

    ! Check head status
    if(self % head <= 0) call fatalError(Here,'Reading with unset (not +ve) head')

    ! Read value
    ptr = self % head
    r = self % XSS(ptr : ptr + N-1)

    ! Advance head
    self % head = self % head + N

  end function readRealArray

  !!
  !! Reads integer under head and does NOT advance
  !!
  function readIntNotAdvance(self) result(i)
    class(aceCard),intent(in)  :: self
    integer(shortInt)          :: i
    character(100),parameter   :: Here ='readIntNotAdvance (aceCard_class.f90)'

    ! Check head status
    if(self % head <= 0) call fatalError(Here,'Reading with unset (not +ve) head')

    ! Read value
    i = real2Int(self % XSS(self % head),Here)

  end function readIntNotAdvance

  !!
  !! Reads real under head and does NOT advance
  !!
  function readRealNotAdvance(self) result(r)
    class(aceCard), intent(in) :: self
    real(defReal)              :: r
    character(100),parameter   :: Here ='readRealNotAdvance (aceCard_class.f90)'

    ! Check head status
    if(self % head <= 0) call fatalError(Here,'Reading with unset (not +ve) head')

    ! Read value
    r = self % XSS(self % head)

  end function readRealNotAdvance

  !!
  !! Advances head by N
  !! Moves forward for +ve N
  !! Moves backwards for -ve N
  !!
  subroutine advanceHead(self,N)
    class(aceCard), intent(inout) :: self
    integer(shortInt), intent(in) :: N

    self % head = self % head + N

  end subroutine advanceHead

  !!
  !! Resets head position back to 0
  !!
  subroutine resetHead(self)
    class(aceCard), intent(inout) :: self

    self % head = 0

  end subroutine resetHead

  !!
  !! Sets head to position of a single energy mu pdf
  !! Position is supplied by addr( LC(J) in Table F-12 of Appendix F )
  !!
  subroutine setToAnglePdf(self,addr)
    class(aceCard), intent(inout) :: self
    integer(shortInt), intent(in) :: addr

    self % head = self % JXS(9) + addr - 1

  end subroutine setToAnglePdf

  !!
  !! Sets head to position of a single energy law
  !! Position is supplied by addrr ( LNW in Table F-14 of Appendix F )
  !!                               ( L,LC in Table F-14 l. of Appendix F )
  !!
  subroutine setToEnergyLaw(self,addr)
    class(aceCard), intent(inout) :: self
    integer(shortInt), intent(in) :: addr

    self % head = self % JXS(11) + addr -1

  end subroutine setToEnergyLaw

  !!
  !! Returns .true. if nuclide has probability table data
  !!
  function hasProbTab(self) result(doesIt)
    class(aceCard), intent(in) :: self
    logical(defBool)           :: doesIt

    doesIt = (self % JXS(23) /= 0)

  end function hasProbTab

  !!
  !! Sets head to position of probability table data
  !!
  subroutine setToProbTab(self)
    class(aceCard), intent(inout) :: self

    self % head = self % JXS(23)

  end subroutine setToProbTab

  !!
  !! Load data for every MT reaction type
  !! NOTE: in ACE format reactions with no secondary neutrons are at the end of MT numbers list.
  !!       As the result it is safe to read data with MTdata(1:NMTs).
  !!
  subroutine setMTdata(self)
    class(aceCard),intent(inout) :: self
    integer(shortInt)            :: NMT  ! number of MT reactions
    integer(shortInt)            :: NMTs ! number of MT reactions with secondary neutrons
    integer(shortInt)            :: ptr  ! pointer on XSS table
    integer(shortInt)            :: i
    integer(shortInt)            :: LOCB ! local LOCB
    character(100), parameter :: Here ='setMTdata (aceCard_class.f90)'

    ! Clean Previous MT-related data
    ! Necessary to ensure that reloading is correct
    if(allocated(self % MTdata)) deallocate(self % MTdata)

    ! Get number of MT reaction
    NMT  = self % NXS(4)
    NMTs = self % NXS(5)

    ! Allocate space
    allocate(self % MTdata(NMT))

    ! Load all MT numbers
    ptr = self % JXS(3)
    self % MTdata(:) % MT = real2Int( self % XSS(ptr : ptr+NMT-1),Here )

    ! Load all Q-values
    ptr = self % JXS(4)
    self % MTdata(:) % Q = self % XSS(ptr : ptr+NMT-1)

    ! Load all raw TY values
    ptr = self % JXS(5)
    self % MTdata(:) % TY = real2Int( self % XSS(ptr : ptr+NMT-1),Here )

    ! Set coordinate frame for MT
    ! If raw TY is -ve, set CMframe to .true. and convert TY to +ve
    do i=1,NMT
      if( self % MTdata(i) % TY < 0 ) then
        self % MTdata(i) % TY        = abs(self % MTdata(i) % TY)
        self % MTdata(i) % CMframe   = .true.
        self % MTdata(i) % isCapture = .false.

      else if(self % MTdata(i) % TY == 0) then
        self % MTdata(i) % isCapture = .true.

      else
        self % MTdata(i) % CMframe   = .false.
        self % MTdata(i) % isCapture = .false.

      end if
    end do

    ! Load raw pointers to XS data (LOCA in MCNP Manual)
    ptr = self % JXS(6)
    self % MTdata(:) % XSp = real2Int( self % XSS(ptr : ptr+NMT-1),Here )

    ! Change relation if raw XS pointers to XSS table (from SIG block)
    self % MTdata(:) % XSp = self % MTdata(:) % XSp + self % JXS(7)

    ! Load number of xs points
    do i=1,NMT
      self % MTdata(i) % N_xs = real2Int( self % XSS( self % MTdata(i) % XSp),Here)
    end do

    ! Load first energy grid index for MT reactions
    do i=1,NMT
      self % MTdata(i) % IE = real2Int( self % XSS( self % MTdata(i) % XSp-1),Here)
    end do

    ! Move XS pointer to point to first entery with XS data
    self % MTdata(:) % XSp = self % MTdata(:) % XSp + 1

    ! Read LOCB data in LAND block
    ptr = self % JXS(8) + 1
    self % MTdata(1:NMTs) % LOCB = real2Int( self % XSS(ptr : ptr+NMTs-1),Here )

    ! Assign pointers to angular distribution data
    do i=1,NMTs
      ! Read LOCB
      LOCB = self % MTdata(i) % LOCB

      select case(LOCB)
        case(0)  ! Isotropic scattering, No data is given
          self % MTdata(i) % isotropic = .true.

        case(-1) ! Correleted scattering
          self % MTdata(i) % correl = .true.

        case(1:huge(i)) ! Positive integer. Can assign pointer to angular data
          self % MTdata(i) % ANDp = self % JXS(9) + LOCB - 1

        case(unINIT) ! default case for absorbtion reactions
          ! Do nothing

        case default ! Should never happen
          call fatalError(Here,'For some reason LOCB is -ve and diffrent from unINIT. WTF?')

      end select
    end do

    ! Read raw LOCC parameters (locators of energy distibutions on DLW block - MCNP manual)
    ptr = self % JXS(10)
    self % MTdata(1:NMTs) % LOCC = real2Int( self % XSS(ptr : ptr+NMTs-1),Here )

  end subroutine setMTdata

  !!
  !! Loads data related to fission
  !!
  subroutine setFissionData(self)
    class(aceCard), intent(inout) :: self
    integer(shortInt)             :: ptr
    integer(shortInt)             :: KNU          ! Pointer to prompt/total Nu data
    character(100), parameter :: Here ='setFissionData (aceCard_class.f90)'

    ! Clean Previous fission-related data
    ! Necessary to ensure that reloading is correct
    self % isFiss       = .false.
    self % hasFISBlock  = .false.
    self % fissIE       = unINIT
    self % fissNE       = unINIT
    self % fissXSp      = unINIT
    self % promptNUp    = unINIT
    self % delayNUp     = unINIT
    self % totalNUp     = unINIT

    ! Check is the nuclide is fissile
    self % isFiss = (self % JXS(2) /= 0)

    ! Check if it has FIS block with XS
    self % hasFISBlock = (self % JXS(21) /= 0)

    ! Exit if nuclide is not fissile
    if (.not. self % isFiss) return

    ! Read data releted to FIS block
    if (self % hasFisBlock) then
      ptr = self % JXS(21)
      self % fissIE  = real2Int( self % XSS(ptr),Here)
      self % fissNE  = real2Int( self % XSS(ptr+1),Here)
      self % fissXSp = self % JXS(21) + 2

    end if

    ! Read KNU pointer to prompt/total NU data
    ptr = self % JXS(2)
    KNU = real2Int(self % XSS(ptr),Here)

    if(KNU > 0 ) then ! Only single promp/total NU data is given
      self % totalNUp = self % JXS(2)

    else if(KNU < 0) then ! Both total and prompt NU data is given
      self % promptNUp = self % JXS(2) + 1
      self % totalNUp  = self % JXS(2) + abs(KNU) + 1

    else
      call fatalError(Here,'KNU is equal to 0. Function should have alrady returned.')

    end if

    ! Read data related to deleyed neutron emissions
    if(self % JXS(24) > 0 ) then ! Delayd NU data is present
      self % delayNUp = self % JXS(24)

    end if

  end subroutine setFissionData

  !!
  !! Read ACE card from file in provided filePath that begins at provided lineNum
  !!
  subroutine readFromFile(self,filePath,lineNum)
    class(aceCard), intent(inout)  :: self
    character(*), intent(in)       :: filePath
    integer(shortInt), intent(in)  :: lineNum
    integer(shortInt)              :: aceFile = 8
    integer(shortInt)              :: i, xssLen
    character(13)                  :: skip
    character(100),parameter :: Here ='readFromFile (aceCard_class.f90)'

    ! Open file to read data
    ! If a path has whitespace at LHS, file will fail to open. We need to trim the whitespace.
    call openToRead(aceFile, adjustl(filePath))

    ! Skip lines
    select case(lineNum)
      case(:0)
        call fatalError(Here,'Invalid line number. 0 or -ve' )

      case(1)
        !Do nothing

      case(2)
        read(aceFile,*)

      case(3:)
        write(skip,'(A1,I10,A2)') '(',lineNum-2,'/)'
        read(aceFile,skip)

      end select

    ! Old skip lines for reference with a do loop
      !if (lineNum > 1) then
      !  do i = 1, lineNum-1
      !    read(aceFile,*)
      !  end do
      !endif
    ! Using write(unit,('numLines/')) to skip lines is much faster.
    ! Needs to be lineNum-2 becouse it skips lineNum-2 lines and then reads another line
    ! with read statement. lineNum-1 is the target line becouse in Fortran indices start with 1

    ! Read Header information
    read(aceFile,'(A10, F12.6, E12.0, 1X, A10)') self % ZAID, self % AW, self % TZ, self % HD
    read(aceFile,'(A70, A10)') self % HK, self % HM

    ! Skip enteries for IZ(I) and AW(I) tabels -> they are legacy empty entery
    do i=1,4
      read(aceFile,*)
    end do

    ! Read NXS, JXS and XSS data tables
    read(aceFile, '(8I9)') self % NXS
    read(aceFile, '(8I9)') self % JXS

    ! Make an attempt to reuse the memory for XSS

    ! Get current length of XSS array
    if(allocated(self % XSS)) then
      xssLen = size(self % XSS)
    else
      xssLen = 0
    end if

    if(self % NXS(1) > xssLen) then
      if(allocated(self % XSS)) deallocate(self % XSS)
      allocate(self % XSS(self % NXS(1)))
    end if

    ! Read XSS array
    read(aceFile,*) self % XSS(1:self % NXS(1))

    ! Close input file
    close(aceFile)

    call self % setMTdata()
    call self % setFissionData()

  end subroutine

  !!
  !! Print contents to the screen.
  !! For DEBUG. Is not pretty.
  !!
  subroutine print_aceCard(self)
    class(aceCard),intent(in) :: self
    integer(shortInt)         :: i

    print *, 'MT REACTION DATA:'
    do i=1,size(self % MTdata)
      call self % MTdata(i) % print()
    end do

    print *, 'NUCLIDE DATA: '
    print *, 'ACE HEADER DATA: '
    print *, self % ZAID, self % AW, self % TZ, self % HD, self % HK, self % HM
    print *, 'FISSION related DATA: '
    print *, self % isFiss, self % fissIE, self % fissNE, self % fissXSp, self % promptNUp, &
             self % delayNup, self % totalNup

  end subroutine print_aceCard

  !!
  !! Print contents of MT reaction data to screen
  !! For DEBUG. Is not pretty.
  !!
  subroutine print_MTreaction(self)
    class(MTreaction),intent(in) :: self

    print *, self % MT, self % Q, self % TY, self % isCapture ,self % CMframe, self % IE, &
             self % N_xs, self % XSp, self % LOCB, self % ANDp, self % LOCC, self % isotropic, self % correl

  end subroutine print_MTreaction

  !!
  !! Returns index of a given MT reaction in MTdata table
  !!
  function getMTidx(self,MT) result(idx)
    class(aceCard),intent(in)     :: self
    integer(shortInt), intent(in) :: MT
    integer(shortInt)             :: idx
    character(100),parameter :: Here = 'getMTidx ( aceCard_class.f90)'

    idx = linFind( self % MTdata(:) % MT, MT)
    if(idx == targetNotFound) call fatalError(Here,'Given MT is not present in ACE card')
    call searchError(idx,Here)

  end function getMTidx

end module aceCard_class
