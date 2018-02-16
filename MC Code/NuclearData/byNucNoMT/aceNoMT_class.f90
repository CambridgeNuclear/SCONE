module aceNoMT_class

  use numPrecision
  use endfConstants
  use genericProcedures ,      only : openToRead, binarySearch, interpolate, searchError
  use RNG_class,               only : RNG
  use emissionFromACE_func,    only : emissionFromACE
  use emissionENDF_class,      only : emissionENDF, emissionENDF_ptr
  use xsEnergyPointNoMT_class, only : xsEnergyPointNoMT
  use materialDataNoMT_class,  only : materialDataNoMT


  implicit none
  private

  integer(shortInt), parameter :: captureIdx  = 1, &
                                  escatterIdx = 2, &
                                  fissionIdx  = 3
  !!
  !! Class that stores nuclide reaction cross-sections read from ACE data card.
  !!
  type, public :: aceNoMT
    real(defReal),public      :: atomWeight    !! Atomic weight ratio A/m_neutron
    real(defReal),public      :: temp          !! Temperature of nucleide [MeV]
    character(zzIdLen),public :: zzId          !! Nuclide ZZid ie. ZZZAAA.nnC
    character(10)             :: date          !! Date the data were processed
    character(70)             :: comment       !! Quick! Guess what is it!
    character(10)             :: MATid         !! MAT indentifier (see p. F-9 Appendic F MCNP 4 manual)
!    integer(shortInt),public  :: nReact        !! number of reaction channels


    logical(defBool)   :: isFissile = .false.

!    integer(shortInt),dimension(:), pointer         :: xsMT  => null()  !! MT numbers of cross-sections in xs
    type(emissionENDF_ptr),dimension(3)             :: emissionData
    real(defReal),dimension(:),allocatable          :: energyGrid         !! Energy grid for xs data
    type(xsEnergyPointNoMT),dimension(:),pointer    :: xsData => null()


  contains
    procedure :: init
    procedure :: energyIdxFor

    procedure,private :: readAceLibrary
    procedure,private :: readXS

  end type aceNoMT

contains

  !!
  !! Searches for the given energy and returns index of the bottom of the interval
  !!
  function energyIdxFor(self,E) result(idx)
    class(aceNoMT), intent(in)  :: self
    real(defReal), intent(in)   :: E
    integer(shortInt)           :: idx
    character(100), parameter   :: Here = 'searchFor (aceNoMT_class.f90)'

    idx = binarySearch(self % energyGrid,E)
    call searchError(idx,Here)

  end function energyIdxFor

  !!
  !! Load reaction cross-section data from ACE file.
  !!
  subroutine init(self,filePath,line)
    class(aceNoMT), intent(inout)           :: self
    character(*), intent(in)                :: filePath  !! Path to file with ACE data card
    integer(shortInt), intent(in)           :: line      !! Line at which the ACE data card begins in the file
    integer(shortInt), dimension(16)        :: NXS
    integer(shortInt), dimension(32)        :: JXS
    real(defReal),dimension(:), allocatable :: XSS

    call self % readAceLibrary(filePath, line, NXS, JXS, XSS)
    call self % readXS(NXS,JXS,XSS)

  end subroutine init

  !!
  !! Read ACE library tabels from the provided filepath and first line number
  !!
  subroutine readAceLibrary(self,filePath,line,NXS,JXS,XSS)
    class(aceNoMT), intent(inout)                       :: self
    character(*), intent(in)                            :: filePath  !! Path to file with ACE data card
    integer(shortInt), intent(in)                       :: line      !! Line at which the ACE data card begins in the file
    integer(shortInt), dimension(16),intent(out)        :: NXS
    integer(shortInt), dimension(32),intent(out)        :: JXS
    real(defReal),dimension(:), allocatable,intent(out) :: XSS
    character(pathLen)                                  :: localFilePath
    integer(shortInt),parameter                         :: aceFile = 66
    integer(shortInt)                                   :: i

    localFilePath = trim(adjustl(filePath))
    call openToRead(aceFile,localFilePath)

    ! Skip lines
    if (line > 1) then
      do i = 1, line-1
        read(aceFile,*)
      end do
    endif

    ! Read Header information
    read(aceFile,'(A10, F12.6, E12.0, 1X, A10)') self % zzId, &
                                                 self % atomWeight,&
                                                 self % temp, &
                                                 self % date
    read(aceFile,'(A70, A10)') self % comment, &
                               self % MATid

    ! Skip enteries for IZ(I) and AW(I) tabels -> they are legacy empty entery
    do i=1,4
      read(aceFile,*)
    end do

    ! Read NXS, JXS and XSS tables
    read(aceFile, '(8I9)') NXS
    read(aceFile, '(8I9)') JXS
    if (allocated(XSS)) deallocate(XSS)
    allocate(XSS(NXS(1)))
    read(aceFile,*) XSS

    close(aceFile)

  end subroutine readAceLibrary

  !!
  !! Read cross-sections into type variables from NXS,JXS and XSS tabels
  !! This implementation ignores reaction other then total,elastic scattering, total capture and
  !! fission.
  !!
  subroutine readXS(self,NXS,JXS,XSS)
    class(aceNoMT), intent(inout)                       :: self
    integer(shortInt), dimension(16),intent(in)         :: NXS
    integer(shortInt), dimension(32),intent(in)         :: JXS
    real(defReal), dimension(:),intent(in)              :: XSS

    real(defReal), dimension(:), allocatable    :: xsEScatter
    real(defReal), dimension(:), allocatable    :: xsCapture
    real(defReal), dimension(:), allocatable    :: xsFission
    integer(shortInt)                           :: reactionNum
    integer(shortInt)                           :: i,j
    integer(shortInt)                           :: firstIdx, numIdx

    if (allocated(self % energyGrid)) deallocate(self % energyGrid)
    if (associated(self % xsData)) deallocate(self % xsData)


    ! Check if nuclide is fissile (or fissonable) -> has fission cross-sections in ACE
    self % isFissile = ( JXS(21) /= 0 )

    ! Allocate space for the energy grid and xsEnergyPoints
    allocate(self % energyGrid (NXS(3)))
    allocate(self % xsData (NXS(3)))

    ! Allocate temporary storage for cross sections
    allocate(xsEScatter(NXS(3) ))
    allocate(xsCapture (NXS(3) ))
    allocate(xsFission (NXS(3) ))

    ! Ensure that all cross-sections are initially 0.0
    xsEScatter = 0.0
    xsCapture  = 0.0
    xsFission  = 0.0

    ! Load Energy Grid
    self % energyGrid = [XSS(JXS(1):JXS(1)+NXS(3)-1)]

    ! Load cross-sections from ESZ block of ACE data card.
    xsCapture  = [XSS(JXS(1)+2*NXS(3):JXS(1)+3*NXS(3)-1)]  ! Total capture XS
    xsEScatter = [XSS(JXS(1)+3*NXS(3):JXS(1)+4*NXS(3)-1)]  ! Elastic Scattering XS

    ! Load fission data if present
    if (self % isFissile) then
      ! First index on energy grid for which fission data is present
      i = int(XSS(JXS(21)))

      ! Number of consecutive eneteries for data
      j = int(XSS(JXS(21)+1))

      xsFission(i:i+j-1) = [XSS(JXS(21)+2:JXS(21)+j+1)]

    end if

    ! Move cross section data into xsEnergyPointType
    call self % xsData % load (xsEScatter,xsCapture,xsFission)

    ! Read emission data
    self % emissionData(captureIdx)   = emissionFromACE(NXS,JXS,XSS,N_gamma)
    self % emissionData(escatterIdx)  = emissionFromACE(NXS,JXS,XSS,N_N_elastic)

    if (self % isFissile) then
      self % emissionData(fissionIdx) = emissionFromACE(NXS,JXS,XSS,N_fission)

    else
      ! Put placeholder emission type equivalent to capture for non=present fission
      self % emissionData(fissionIdx) = emissionFromACE(NXS,JXS,XSS,N_gamma)

    end if

  end subroutine readXS


end module aceNoMT_class



