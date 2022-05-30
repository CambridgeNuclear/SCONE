module aceSabCard_class

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
  !! ACE DATA CARD PARSER FOR A SINGLE NUCLIDE (e.g. 92238.03c)
  !!
  !!
  !! Please note that header of the ACE DATA CARD is stored in public components of type(aceCard)
  !!
  !! Procedures for MT accept elastic scattering which has MT = N_N_elastic (from endfConstants)
  !!
  type, public, extends(dataDeck) :: aceSabCard
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

    ! RAW ACE TABLES *** PUBLIC in DEBUG * WILL BE PRIVATE
    integer(shortInt),dimension(16),public        :: NXS
    integer(shortInt),dimension(32),public        :: JXS
    real(defReal),dimension(:),allocatable,public :: XSS

    ! Flags
    logical(defBool)    :: hasElasticXs
    logical(defBool)    :: elIsCoherent


    ! Sampling parameters
    integer(shortInt)   :: inelOutgoingE
    integer(shortInt)   :: elOutgoingMu
    real(defReal),dimension(:),allocatable     :: inelCDF

  contains
    ! Superclass procedures
    procedure :: myType

    ! Output cross sections and energy info
    procedure :: ESZ_elastic
    procedure :: ESZ_inelastic
    procedure :: elasticEnergies
    procedure :: inelasticEnergies
    procedure :: isCoherent
    procedure :: isInelContinuous
    procedure :: getCDF
    procedure :: inelOutE
    procedure :: inelOutMu
    procedure :: elOutMu
    procedure :: hasElastic
    procedure :: hasITCA

    ! Procedures to set head
    procedure :: setToElastic
    procedure :: setToInelastic
    procedure :: setToElasticOut
    procedure :: setToInelasticOut

    ! Procedures related to reading under head and head status
    !
    procedure :: readInt             ! Read int under head and advance
    procedure :: readReal            ! Read real under head and advance
    procedure :: readIntArray        ! Read N ints beginning under head and advance by N
    procedure :: readRealArray       ! Read N ints beginning under head and advance by N

    procedure :: advanceHead         ! Move head by an integer +ve forward -ve backward
    procedure :: setRelativeTo

    ! Initialisation procedure
    procedure :: readFromFile

    ! Private procedures
    procedure,private :: setInelastic
    procedure,private :: setElastic

  end type aceSabCard

contains
  !!
  !! Return String with type name
  !!
  !! See dataDeck_inter for details
  !!
  pure function myType(self) result(type)
    class(aceSabCard), intent(in) :: self
    character(:),allocatable      :: type

    type = 'aceSabCard_class'

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
  !! Returns data from ESZ block
  !!
  !! Hopefully will not blow up the stack
  !!
  !! Following request can be delivered (CASE SENSITIVE):
  !!  'energyGrid'   -> energy grid
  !!  'Pvalues'      -> table entry. It has different meaning if scattering is coherent or incoherent
  !!
  !! Args:
  !!   request [in] -> character with the XS request
  !!
  !! Error:
  !!   fatalError if request is unrecognised
  !!
  function ESZ_elastic(self, request) result(array)
    class(aceSabCard), intent(in)                     :: self
    character(*), intent(in)                          :: request
    real(defReal),dimension(self % elasticEnergies()) :: array
    integer(shortInt)                                 :: ptr, N
    character(100), parameter :: Here = 'ESZ_elastic (aceSabCard_class.f90)'

    ! Obtain size of the XS grid
    N = self % elasticEnergies()

    select case(request)
      case('energyGrid')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(4) + 1

      case('Pvalues')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(4) + N + 1

      case default
        call fatalError(Here,'Unrecognised request string: |' // request //'|')
        ptr = 0
      end select

      ! Return data
      array = self % XSS(ptr : ptr + N-1)

  end function ESZ_elastic

  !!
  !! Returns data from ESZ block
  !!
  !! Hopefully will not blow up the stack
  !!
  !! Following request can be delivered (CASE SENSITIVE):
  !!  'energyGrid'   -> energy grid
  !!  'inelasticXS'      -> inelastic cross section
  !!
  !! Args:
  !!   request [in] -> character with the XS request
  !!
  !! Error:
  !!   fatalError if request is unrecognised
  !!
  function ESZ_inelastic(self, request) result(array)
    class(aceSabCard), intent(in)                       :: self
    character(*), intent(in)                            :: request
    real(defReal),dimension(self % inelasticEnergies()) :: array
    integer(shortInt)                                   :: ptr, N
    character(100), parameter :: Here = 'ESZ_inelastic (aceSabCard_class.f90)'

    ! Obtain size of the XS grid
    N = self % inelasticEnergies()

    select case(request)
      case('energyGrid')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + 1

      case('inelasticXS')
        ! Set pointer to approperiate place in XSS
        ptr = self % JXS(1) + N + 1

      case default
        call fatalError(Here,'Unrecognised request string: |' // request //'|')
        ptr = 0
    end select

    ! Return data
    array = self % XSS(ptr : ptr + N-1)

  end function ESZ_inelastic


  !!
  !! Return number of elastic scattering ingoing energies groups
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of ingoing elastic scattering energies
  !!
  pure function elasticEnergies(self) result(N)
    class(aceSabCard), intent(in) :: self
    integer(shortInt)             :: N

    N = self % XSS(self % JXS(4))

  end function elasticEnergies

  !!
  !! Return number of inelastic scattering ingoing energies groups
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of ingoing inelastic scattering energies
  !!
  pure function inelasticEnergies(self) result(N)
    class(aceSabCard), intent(in) :: self
    integer(shortInt)             :: N

    N = self % XSS(self % JXS(1))

  end function inelasticEnergies

  !!
  !! Sets read head to beginning of elastic scattering energies
  !! Returns error is nuclide doesn't have data for elastic scattering
  !!
  subroutine setToElastic(self)
    class(aceSabCard), intent(inout) :: self
    character(100),parameter :: Here ='setToElastic (aceSabCard_class.f90)'

    if(self % hasElasticXs) then
      self % head = self % JXS(4) + 1

    else
      call fatalError(Here,'Elastic cross sections are not present!')

    end if

  end subroutine setToElastic

  !!
  !! Sets read head to beginning of inelastic scattering energies
  !!
  subroutine setToInelastic(self)
    class(aceSabCard), intent(inout) :: self

    self % head = self % JXS(1) + 1

  end subroutine setToInelastic

  !!
  !! Sets read head to beginning of elastic scattering energies
  !! Returns error is nuclide doesn't have data for elastic scattering
  !!
  subroutine setToElasticOut(self)
    class(aceSabCard), intent(inout) :: self
    character(100),parameter :: Here ='setToElastic (aceSabCard_class.f90)'

    self % head = self % JXS(6)

  end subroutine setToElasticOut

  !!
  !! Sets read head to beginning of inelastic scattering outgoing data
  !!
  subroutine setToInelasticOut(self)
    class(aceSabCard), intent(inout) :: self

    self % head = self % JXS(3)

  end subroutine setToInelasticOut

  !!
  !!
  !!
  pure  function isCoherent(self) result(isIt)
    class(aceSabCard), intent(in) :: self
    logical(defBool)              :: isIt

    isIt = self % elIsCoherent

  end function isCoherent

  !!
  !!
  !!
  pure  function hasElastic(self) result(hasIt)
    class(aceSabCard), intent(in) :: self
    logical(defBool)              :: hasIt

    hasIt = self % hasElasticXs

  end function hasElastic

  !!
  !!
  !!
  pure  function hasITCA(self) result(doesIt)
    class(aceSabCard), intent(in) :: self
    logical(defBool)              :: doesIt

    doesIt = (self % JXS(4) /= 0 .and. self % NXS(6) /= -1)

  end function hasITCA

  !!
  !!
  !!
  pure  function isInelContinuous(self) result(isIt)
    class(aceSabCard), intent(in) :: self
    logical(defBool)              :: isIt

    isIt = .false.
    if (self % NXS(7) == 2) isIt = .true.

  end function isInelContinuous

  !!
  !!
  !!
  pure  function getCDF(self) result(CDF)
    class(aceSabCard), intent(in) :: self
    real(defReal),dimension(self % inelOutgoingE + 1)    :: CDF

    CDF = self % inelCDF

  end function getCDF

  !!
  !!
  !!
  pure  function inelOutE(self) result(N)
    class(aceSabCard), intent(in) :: self
    integer(shortInt)             :: N

    N = self % inelOutgoingE

  end function inelOutE

  !!
  !!
  !!
  pure  function inelOutMu(self) result(N)
    class(aceSabCard), intent(in) :: self
    integer(shortInt)             :: N

    N = self % NXS(3)

  end function inelOutMu

  !!
  !!
  !!
  pure  function elOutMu(self) result(N)
    class(aceSabCard), intent(in) :: self
    integer(shortInt)             :: N

    N = self % elOutgoingMu

  end function elOutMu


  !!
  !! Read single integer and advance read head
  !! Return error if under head is not integer
  !!
  function readInt(self) result(i)
    class(aceSabCard), intent(inout) :: self
    integer(shortInt)                :: i
    character(100),parameter :: Here ='readInt (aceSabCard_class.f90)'

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
    class(aceSabCard), intent(inout) :: self
    real(defReal)                    :: r
    character(100),parameter :: Here ='readReal (aceSabCard_class.f90)'

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
    class(aceSabCard), intent(inout)  :: self
    integer(shortInt), intent(in)     :: N
    integer(shortInt),dimension(N)    :: i
    integer(shortInt)                 :: ptr
    character(100),parameter :: Here ='readIntArray (aceSabCard_class.f90)'

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
    class(aceSabCard), intent(inout)  :: self
    integer(shortInt), intent(in)     :: N
    real(defReal),dimension(N)        :: r
    integer(shortInt)                 :: ptr
    character(100),parameter :: Here ='readRealArray (aceSabCard_class.f90)'

    ! Check head status
    if(self % head <= 0) call fatalError(Here,'Reading with unset (not +ve) head')

    ! Read value
    ptr = self % head
    r = self % XSS(ptr : ptr + N-1)

    ! Advance head
    self % head = self % head + N

  end function readRealArray

  !!
  !! Advances head by N
  !! Moves forward for +ve N
  !! Moves backwards for -ve N
  !!
  subroutine advanceHead(self,N)
    class(aceSabCard), intent(inout) :: self
    integer(shortInt), intent(in)    :: N

    self % head = self % head + N

  end subroutine advanceHead

  !!
  !! Sets head to position given by Root and and offset
  !!
  !! Head = root + offset
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
    class(aceSabCard), intent(inout) :: self
    integer(shortInt), intent(in)    :: root
    integer(shortInt), intent(in)    :: offset
    integer(shortInt)                :: pos
    character(100), parameter :: Here = 'setRelativeTo (aceSabCard_class.f90)'

    ! Validate root
    if (root < 1 .or. root > self % NXS(1)) then
      call fatalError(Here,trim(self % ZAID)//' root is outside bounds: '//numToChar(root))
    end if

    pos = root + offset

    ! Validate position
    if(pos < 1 .or. pos > self % NXS(1)) then
      call fatalError(Here,trim(self % ZAID)//' position reached outside the bounds')
    end if

    self % head = pos

  end subroutine setRelativeTo

  !!
  !! Load data for every MT reaction type
  !! NOTE: in ACE format reactions with no secondary neutrons are at the end of MT numbers list.
  !!       As the result it is safe to read data with MTdata(1:NMTs).
  !!
  subroutine setInelastic(self)
    class(aceSabCard),intent(inout) :: self
    integer(shortInt)               :: i
    character(100), parameter :: Here ='setInelastic (aceSabCard_class.f90)'

    self % inelOutgoingE = self % NXS(4)

    if (allocated(self % inelCDF)) deallocate(self % inelCDF)
    allocate(self % inelCDF(self % inelOutgoingE + 1))
    ! Initialise the CDF to one everywhere
    self % inelCDF = ZERO

    select case(self % NXS(7))
      case(0)
        ! Fill up the vector with a uniform distribution
        do i = 2,(self % inelOutgoingE + 1)
          self % inelCDF(i) = self % inelCDF(i-1) + ONE
        end do

      case(1)
        ! Fill up the vector with a skewed distribution
        self % inelCDF(2) = ONE
        self % inelCDF(3) = 5.0_defReal
        do i = 4,(self % inelOutgoingE - 1)
          self % inelCDF(i) = self % inelCDF(i-1) + 10.0_defReal
        end do
        self % inelCDF(self % inelOutgoingE) = self % inelCDF(self % inelOutgoingE - 1) + 4.0_defReal
        self % inelCDF(self % inelOutgoingE + 1) = self % inelCDF(self % inelOutgoingE) + ONE

      case(2)

        self % inelCDF(self % inelOutgoingE + 1) = ONE

      case default
        call fatalError(Here, 'Flag for inelastic outgoing energy/angle should be 0, 1 or 2')

      end select

    ! Normalise the distribution and check that the final value is one
    self % inelCDF = self % inelCDF/self % inelCDF(self % inelOutgoingE + 1)
    if (self % inelCDF(self % inelOutgoingE + 1) /= ONE) call fatalError(Here,'CDF does not add up')

  end subroutine setInelastic

  !!
  !! Loads data related to fission
  !!
  subroutine setElastic(self)
    class(aceSabCard), intent(inout) :: self
    character(100), parameter :: Here ='setElastic (aceSabCard_class.f90)'

    self % hasElasticXs = (self % JXS(4) /= 0)
    self % elIsCoherent = (self % NXS(5) == 4)

    if (self % hasElasticXs) then
      self % elOutgoingMu = self % NXS(6) + 1
    end if

  end subroutine setElastic

  !!
  !! Read ACE card from dile in provided filePath that beggins at provided lineNum
  !!
  subroutine readFromFile(self,filePath,lineNum)
    class(aceSabCard), intent(inout)  :: self
    character(*), intent(in)          :: filePath
    integer(shortInt), intent(in)     :: lineNum
    integer(shortInt)                 :: aceFile = 8
    character(pathLen)                :: localFilePath
    integer(shortInt)                 :: i, xssLen
    character(13)                     :: skip
    character(100),parameter :: Here ='readFromFile (aceSabCard_class.f90)'

    ! Copy filepath and make shure it is left adjusted
    localFilePath = trim(adjustl(filePath))

    ! Open file to read data
    call openToRead(aceFile,localFilePath)

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

    call self % setInelastic()
    call self % setElastic()

  end subroutine readFromFile

end module aceSabCard_class
