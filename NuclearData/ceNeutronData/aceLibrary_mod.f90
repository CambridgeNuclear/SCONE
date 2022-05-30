!!
!! This module is a factory that produces aceCards for a given nuclide
!!
!! Library format:
!!   '!' line comment indicators
!!
module aceLibrary_mod

  use numPrecision
  use iso_fortran_env,   only : IOSTAT_END
  use genericProcedures, only : fatalError, replaceChar, charToInt
  use charLib_func,      only : splitChar
  use charMap_class,     only : charMap
  use aceCard_class,     only : aceCard
  use aceSabCard_class,  only : aceSabCard

  implicit none
  private

  !! Module parameters
  integer(shortInt), parameter         :: UNDEF   = 0
  integer(shortInt), parameter         :: ACE_CE  = 1
  integer(shortInt), parameter         :: ACE_SAB = 2
  integer(shortInt), parameter         :: MAX_COL = 900
  character(*),parameter               :: READ_FMT = '(A900)'
  character(1),      parameter         :: COMMENT_TOKEN = '!'

  character(1), parameter :: TAB   = char(9)
  character(1), parameter :: SPACE = ' '

  !!
  !! Storage for data associated with each library entry
  !!
  type, private :: item
    character(nameLen) :: ZAID
    integer(shortInt)  :: firstLine = 1
    integer(shortInt)  :: type = UNDEF
    character(pathLen) :: path
  end type

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Module Members
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  type(item),dimension(:),allocatable :: entry
  type(charMap)                       :: map
  character(pathLen)                  :: libFile

  public :: load
  public :: new_neutronACE
  public :: new_moderACE
  public :: kill

contains

  !!
  !! Load library from file
  !!
  !! Args:
  !!   path [in] -> Path to library file. Should be absolute path for safety
  !!
  !! Erros:
  !!   fatalError or Frotran intrinsic read  error if there are ill-formated lines in
  !!     provided library file
  !!   Fortran intrinsic error if the requested file does not exist
  !!
  subroutine load(path)
    character(*), intent(in) :: path
    integer(shortInt)        :: library
    integer(shortInt)        :: i, errorCode, libLen, readStat, last
    character(99)            :: errorMsg
    character(MAX_COL)       :: buffor
    character(100),parameter :: Here = 'load (aceLibrary_mod.f90)'

    ! Clean
    call kill()

    libFile = path

    ! Open file for reading
    open(newunit = library, &
         file = path, &
         status='old', &
         action='read', &
         iostat = errorCode, &
         iomsg = errorMsg)

    if(errorCode /= 0) call fatalError(Here, "File Error: "//trim(adjustl(errorMsg)))

    ! Find number of entries in the library
    libLen = 0
    do
      read(unit = library, fmt=READ_FMT, iostat=readStat) buffor
      if(readStat == IOSTAT_END) exit

      ! Preform line preprocessing
      call preprocessLine(buffor)
      if(len_trim(buffor) /= 0) libLen = libLen + 1
    end do
    rewind(library)

    ! Allocate space for the library
    allocate(entry(libLen))
    call map % init(libLen)

    ! Load Library information
    i = 1
    do
      read(unit = library, fmt=READ_FMT, iostat=readStat) buffor
      if(readStat == IOSTAT_END) exit

      ! Preform line preprocessing
      call preprocessLine(buffor)

      ! Read if line is not empty
      if( len_trim(buffor) /= 0) then
        associate ( bounds => splitChar(buffor, ' ') )
          if (size(bounds,2) < 3) call fatalError(Here, 'Ill formatted line: '//trim(buffor))
          ! Read Content
          entry(i) % ZAID      = buffor(bounds(1,1):bounds(2,1))
          entry(i) % firstLine = charToInt(buffor(bounds(1,2):bounds(2,2)))
          entry(i) % path      = buffor(bounds(1,3):len_trim(buffor))
          i = i +1
        end associate
      end if
    end do

    ! Detect type of each entry and store in map for quick access
    do i=1,size(entry)
      last = len_trim(entry(i) % ZAID)

      ! Set entry type
      select case(entry(i) % ZAID(last:last))
        case ('c')
          entry(i) % TYPE = ACE_CE

        case ('t')
          entry(i) % TYPE = ACE_SAB

        case default
          call fatalError(Here,'Unrecognised ACE CARD type: '// entry(i) % ZAID(last:last))

      end select

      ! Remove type information from the entry
      entry(i) % ZAID = entry(i) % ZAID(1:last-1)

      ! Add to map
      call map % add(entry(i) % ZAID, i)

    end do

    ! Clean up
    close(library)

  end subroutine load

  !!
  !! Load new CE Neutron XSs data for a given ZAID
  !!
  !! Note ZAID is provided without MCNP suffix
  !!   1001.03 -> will work
  !!   1001.03c -> will NOT work
  !!   1001.03d -> will NOT work
  !!
  !! Args:
  !!   ACE [inout] -> ACE Card which will store the data
  !!   ZAID [in]   -> Requested nuclide ZAID of the form ZZAAA.TT
  !!
  !! Errors:
  !!   fatalError if ZAID is not present in the library
  !!   fatalError if ZAID points to data which is not for CE Neutrons, that is does not have a
  !!     MCNP 'c' suffix e.g. {1001.03c}
  !!
  subroutine new_neutronACE(ACE, ZAID)
    class(aceCard), intent(inout)  :: ACE
    character(nameLen), intent(in) :: ZAID
    integer(shortInt)              :: idx
    integer(shortInt), parameter     :: NOT_FOUND = -1
    character(100), parameter :: Here = 'new_neutronACE (aceLibrary_mod.f90)'

    ! Find index of the requested ZAID identifier
    idx = map % getOrDefault(ZAID, NOT_FOUND)
    if(idx == NOT_FOUND) then
      call fatalError(Here, trim(ZAID) //" was not found in ACE library from: "//trim(libFile))
    end if

    ! Verify that type is correct
    if( entry(idx) % type /= ACE_CE) then
      call fatalError(Here,trim(ZAID)//" is not a ACE data with CE XSs.")
    end if

    ! Load data
    call ACE % readFromFile(entry(idx) % path, entry(idx) % firstLine)

  end subroutine new_neutronACE

  !!
  !! Load new thermal scattering neutron data for a given file name
  !!
  !! Note: file name is provided without suffix
  !!   h-h2o.42  -> will work
  !!   h-h2o.42t -> will NOT work
  !!
  !! Args:
  !!   ACE [inout] -> ACE Sab Card which will store the data
  !!   file [in]   -> Requested thermal scattering file
  !!
  !! Errors:
  !!   fatalError if file is not present in the library
  !!   fatalError if file points to data which is not for Sab data
  !!
  subroutine new_moderACE(ACE, file)
    class(aceSabCard), intent(inout) :: ACE
    character(nameLen), intent(in)   :: file
    integer(shortInt)                :: idx
    integer(shortInt), parameter     :: NOT_FOUND = -1
    character(100), parameter :: Here = 'new_moderACE (aceLibrary_mod.f90)'

    ! Find index of the requested ZAID identifier
    idx = map % getOrDefault(file, NOT_FOUND)
    if (idx == NOT_FOUND) then
      call fatalError(Here, trim(file) //" was not found in ACE library from: "//trim(libFile))
    end if

    ! Verify that type is correct
    if( entry(idx) % type /= ACE_SAB) then
      call fatalError(Here,trim(file)//" is not a ACE data with CE XSs.")
    end if

    ! Load data
    call ACE % readFromFile(entry(idx) % path, entry(idx) % firstLine)

  end subroutine new_moderACE


  !!
  !! Returns module to uninitialised state
  !!
  subroutine kill()

    if(allocated(entry)) deallocate(entry)
    call map % kill()
    libFile = ''

  end subroutine kill

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Utility Functions
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Perform preprocessing of a line
  !!   - changes all TABS to SPACES
  !!   - removes any content after (and including) COMMENT_TOKEN
  !!
  !! Args:
  !!   line [in] -> character that contain a line
  !!
  subroutine preprocessLine(line)
    character(*), intent(inout) :: line
    integer(shortInt)           :: i

    do i=1,len(line)
      select case(line(i:i))
        case(TAB)
          line(i:i) = SPACE

        case(COMMENT_TOKEN)
          line(i:len(line)) = ''
          return

      end select
    end do

  end subroutine preprocessLine


end module aceLibrary_mod
