module aceNoMT_class

  use numPrecision
  use genericProcedures ,   only : openToRead
  use RNG_class,            only : RNG
  use emissionFromACE_func, only : emissionFromACE
  use emissionENDF_class,   only : emissionENDF, emissionENDF_ptr


  implicit none
  private

  type, public :: aceNoMT
    !! Class that stores isotopic reaction cross-sections read from ACE data card.
    private
    real(defReal)      :: atomWeight    !! Atomic weight ratio A/m_neutron
    real(defReal)      :: temp          !! Temperature of nucleide [MeV]
    character(zzIdLen) :: zzId          !! Isotope ZZid ie. ZZZAAA.nnC
    character(10)      :: date          !! Date the data were processed
    character(70)      :: comment       !! Quick! Guess what is it!
    character(10)      :: MATid         !! MAT indentifier (see p. F-9 Appendic F MCNP 4 manual)

    logical(defBool)   :: isFissile = .false.

    integer(shortInt),dimension(:), allocatable     :: xsMT         !! MT numbers of cross-sections in xs
    type(emissionENDF_ptr),dimension(:),allocatable :: emissionData
    real(defReal),dimension(:),allocatable          :: energyGrid   !! Energy grid for xs data
    real(defReal),dimension(:,:),allocatable        :: xs           !! Microscpoic cross-sections table order -> xs(reactionMTnumber,energyPoint)
  contains
    procedure :: init

    procedure,private :: readAceLibrary
    procedure,private :: readXS
  end type aceNoMT

contains

  subroutine init(self,filePath,line)
    !! Load reaction cross-section data from ACE file.
    class(aceNoMT), intent(inout)           :: self
    character(*), intent(in)                :: filePath  !! Path to file with ACE data card
    integer(shortInt), intent(in)           :: line      !! Line at which the ACE data card begins in the file
    integer(shortInt), dimension(16)        :: NXS
    integer(shortInt), dimension(32)        :: JXS
    real(defReal),dimension(:), allocatable :: XSS

    call self % readAceLibrary(filePath, line, NXS, JXS, XSS)
    call self % readXS(NXS,JXS,XSS)

  end subroutine

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

  end subroutine



  subroutine readXS(self,NXS,JXS,XSS)
    !! Read cross-sections into type variables from NXS,JXS and XSS tabels
    !! This implementation ignores reaction other then total,elastic scattering, total capture and
    !! fission.
    class(aceNoMT), intent(inout)                       :: self
    integer(shortInt), dimension(16),intent(in)         :: NXS
    integer(shortInt), dimension(32),intent(in)         :: JXS
    real(defReal), dimension(:),intent(in)              :: XSS

    integer(shortInt)      :: reactionNum
    integer(shortInt)      :: i,j

    if (allocated(self % energyGrid)) deallocate(self % energyGrid)
    if (allocated(self % xsMT)) deallocate(self % xsMT)
    if (allocated(self % xs)) deallocate(self % xs)

    ! Check if isotope is fissile (or fissonable) -> has fission cross-sections in ACE
    if (JXS(21) == 0 ) then
      reactionNum = 3
    else
      self % isFissile = .true.
      reactionNum = 4
    end if

    allocate(self % energyGrid (NXS(3)))
    allocate(self % xs(reactionNum, NXS(3)))

    allocate(self % xsMT(reactionNum))
    allocate(self % emissionData(2:reactionNum))

    ! Ensure that all cross-sections are initially 0.0
    self % xs = 0.0

    self % energyGrid = [XSS(JXS(1):JXS(1)+NXS(3)-1)]

    ! Load cross-sections from ESZ block of ACE data card.
    ! NOTICE: Elastic scattering is as second entery despite the fact it is 3rd in ESZ

    self % xs(1,:) = [XSS(JXS(1)+NXS(3):JXS(1)+2*NXS(3)-1)]    ! Total XS
    self % xs(2,:) = [XSS(JXS(1)+3*NXS(3):JXS(1)+4*NXS(3)-1)]  ! Elastic Scattering XS
    self % xs(3,:) = [XSS(JXS(1)+2*NXS(3):JXS(1)+3*NXS(3)-1)]  ! Total capture XS

    ! Assign MT numbers
    self % xsMT(1:3) = [ 1, 2, 101]

    ! Load fission data if present
    if (self % isFissile) then
      i=int(XSS(JXS(21)))   ! First index on energy grid for which fission data is present
      j=int(XSS(JXS(21)+1)) ! Number of consecutive eneteries for data
      self % xs(4,i:i+j-1) = [XSS(JXS(21)+2:JXS(21)+j+1)]
      self % xsMT(4) = 18
    end if

    ! Read emission data
    do i=2,size(self % xsMT)
      self % emissionData(i) = emissionFromACE(NXS,JXS,XSS,self % xsMT(i))
    end do

  end subroutine


end module aceNoMT_class



