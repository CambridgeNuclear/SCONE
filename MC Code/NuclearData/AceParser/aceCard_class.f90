module aceCard_class

  use numPrecision
  use genericProcedures, only : fatalError, openToRead, isInteger, linFind,&
                                targetNotFound, searchError

  implicit none
  private

  interface real2Int
    module procedure real2Int_scalar
    module procedure real2Int_array
  end interface

  integer(shortInt), parameter :: unInit = -17

  type,public :: MTreaction
    integer(shortInt) :: MT        = 0        ! ENDF MT number of reaction
    real(defReal)     :: Q         = -999.0   ! Q-value for reaction
    integer(shortInt) :: TY        = -999     ! Neutron release for reaction +ve integer. 19 means fission
    logical(defBool)  :: CMframe   = .true.   ! Is data for reaction given in CM frame
    integer(shortInt) :: IE        = unINIT   ! First energy grid index
    integer(shortInt) :: N_xs      = unINIT   ! Number of consecutive XS data points
    integer(shortInt) :: XSp       = unINIT   ! Location of cross sections on XSS table
    integer(shortInt) :: LOCB      = unINIT   ! Raw ACE LOCB parameter from MCNP manual
    integer(shortInt) :: ANDp      = unINIT   ! Location of angular distribution data in XSS table
    integer(shortInt) :: DLWp      = unINIT   ! Location of energy distribution data in XSS table
    logical(defBool)  :: isotropic = .false.  ! Indicates isotropic scattering at all energies
    logical(defBool)  :: law44     = .false.  ! Correlated angle-energy scattering with ENDF LAW 44
  contains
    procedure :: print => print_MTreaction
  end type MTreaction


  type, public :: aceCard
    !private
    integer(shortInt)  :: head = 0    ! Current read location on a card

    character(zzIdLen) :: ZAID = ''   ! 10 character name ZZZAAA.nnC
    real(defReal)      :: AW   = -ONE ! Atomic weight ratio. Atomic weight divided by the neutron mass
    real(defReal)      :: TZ   = -ONE ! Temperature at thich data were processed [MeV]
    character(10)      :: HD   = ' '  ! 10 character date when data wre processed
    character(70)      :: HK   = ' '  ! 70 character comment
    character(10)      :: HM   = ' '  ! 10 character MAT indentifier

    type(MTreaction),dimension(:),allocatable :: MTdata

    ! Fission related data
    logical(defBool)  :: isFissile = .false.    ! Flag is true if JXS(21) /= 0
    integer(shortInt) :: fissIE    = unINIT     ! First energy grid index for fission XSs
    integer(shortInt) :: fissNE    = unINIT     ! Number of fission XS points
    integer(shortInt) :: fissXSp   = unINIT     ! Location of first fission xs point on XSS table
    integer(shortInt) :: promptNUp = unINIT     ! Location of prompt NU data
    integer(shortInt) :: delayNUp  = unINIT     ! Location of delay NU data
    integer(shortInt) :: totalNUp  = unINIT     ! Location of total NU data


    integer(shortInt),dimension(16)        :: NXS
    integer(shortInt),dimension(32)        :: JXS
    real(defReal),dimension(:),allocatable :: XSS

  contains
    procedure :: ESZblock

    procedure :: firstIdxMT
    procedure :: numXsPointsMT
    procedure :: xsMT
    ! Return TY for MT
    ! Return Q for MT
    ! Return Ref frame for MT
    ! Return LOCB for MT
    ! Set head to angular for MT + Escater
    ! Set head to energy for MT

    procedure :: readInt
    procedure :: readReal
    procedure :: readIntArray
    procedure :: readRealArray

    procedure :: setMTdata
    procedure :: setFissionData
    procedure :: readFromFile
    procedure :: print => print_aceCard

    procedure,private :: getMTidx

  end type aceCard

contains

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

    if ( isInteger(r) ) then
      i = int(r,shortInt)
    else
      call fatalError(Where,' Tried to read real into integer variable')
    end if

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
  !! Returns data from ESZ block
  !! request string specifies which array is returned:
  !!  'energyGrid'   -> energy grid
  !!  'totalXS'      -> total cross-section
  !!  'absorbtionXS' -> total absorbtion cross-section (without fission)
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

      case('absorbtionXS')
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

      end select

      ! Return data
      data = self % XSS(ptr : ptr + N-1)

  end subroutine ESZblock

  !!
  !! For a given MT number returns first index of its XS at the energy grid
  !!
  function firstIdxMT(self, MT) result(IE)
    class(aceCard),intent(in)    :: self
    integer(shortInt),intent(in) :: MT
    integer(shortInt)            :: IE
    integer(shortInt)            :: idx

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
        self % MTdata(i) % TY = abs(self % MTdata(i) % TY)
        self % MTdata(i) % CMframe = .true.

      else
        self % MTdata(i) % CMframe = .false.

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

        case(-1) ! Correleted scattering with LAW 44
          self % MTdata(i) % law44 = .true.

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
    self % MTdata(1:NMTs) % DLWp = real2Int( self % XSS(ptr : ptr+NMTs-1),Here )

    ! Relate energy distributions pointers to XSS table
    self % MTdata(1:NMTs) % DLWp = self % MTdata(1:NMTs) % DLWp + self % JXS(9) - 1
  end subroutine setMTdata

  !!
  !! Loads data related to fission
  !!
  subroutine setFissionData(self)
    class(aceCard), intent(inout) :: self
    integer(shortInt)             :: ptr
    logical(defBool)              :: hasNoNuBlock
    integer(shortInt)             :: KNU          ! Pointer to prompt/total Nu data
    character(100), parameter :: Here ='setFissionData (aceCard_class.f90)'

    ! Check is the nuclide is fissile
    self % isFissile = (self % JXS(21) /= 0)

    ! Read data releted to FIS block
    if(self % isFissile) then
      ptr = self % JXS(21)
      self % fissIE  = real2Int( self % XSS(ptr),Here)
      self % fissNE  = real2Int( self % XSS(ptr+1),Here)
      self % fissXSp = self % JXS(21) + 2

    end if

    ! Check if NU data is provided
    hasNoNuBlock = (self % JXS(2) == 0)

    ! NU data can be provided even when nuclide is not marked fissle
    if(self % isFissile .and. hasNoNuBlock) then
      call fatalError(Here,'Nuclide is fissile but no NU data is provided. WTF?')

    else if (hasNoNuBlock) then ! Leve the subroutine and do not process the NU data
      return

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
  !! Read ACE card from dile in provided filePath that beggins at provided lineNum
  !!
  subroutine readFromFile(self,filePath,lineNum)
    class(aceCard), intent(inout)  :: self
    character(*), intent(in)       :: filePath
    integer(shortInt), intent(in)  :: lineNum
    integer(shortInt)              :: aceFile
    character(pathLen)             :: localFilePath
    integer(shortInt)              :: i
    character(100)                 :: debug

    ! Copy filepath and make shure it is left adjusted
    localFilePath = trim(adjustl(filePath))

    ! Open file to read data
    call openToRead(aceFile,localFilePath)

    ! Skip lines
    if (lineNum > 1) then
      do i = 1, lineNum-1
        read(aceFile,*)
      end do
    endif

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
    if (allocated(self % XSS)) deallocate(self % XSS)
    allocate(self % XSS(self % NXS(1)))
    read(aceFile,*) self % XSS

    ! Close input file
    close(aceFile)

    call self % setMTdata()
    call self % setFissionData()
  end subroutine

  !!
  !! Print contents to the screen
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
    print *, self % isFissile, self % fissIE, self % fissNE, self % fissXSp, self % promptNUp, &
             self % delayNup, self % totalNup

  end subroutine print_aceCard

  !!
  !! Print contents of MT reaction data to screen
  !!
  subroutine print_MTreaction(self)
    class(MTreaction),intent(in) :: self

    print *, self % MT, self % Q, self % TY, self % CMframe, self % IE, self % N_xs, self % XSp,&
             self % LOCB, self % ANDp, self % DLWp, self % isotropic, self % law44

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
