!!
!! This module is a factory that produces aceCards for a given nuclide
!!
!! Library format:
!!   '//' or '!' line comment indicators
module aceLibrary_mod

  use numPrecision
  use iso_fortran_env,   only : IOSTAT_END
  use genericProcedures, only : fatalError, replaceChar, charToInt
  use charLib_func,      only : splitChar
  use charMap_class,     only : charMap

  implicit none
  private

  !! Module parameters
  integer(shortInt), parameter         :: UNDEF  = 0
  integer(shortInt), parameter         :: ACE_CE = 1
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

contains

  !!
  !! Load library from file
  !!
  subroutine load(path)
    character(*), intent(in) :: path
    integer(shortInt)        :: library
    integer(shortInt)        :: i, errorCode, libLen, readStat, last
    character(99)            :: errorMsg
    character(MAX_COL)       :: buffor
    character(100),parameter :: Here = 'load (aceLibrary_mod.f90)'

    libFile = path

    ! Open file for reading
    open(newunit = library, &
         file = path, &
         status='old', &
         action='read', &
         iostat = errorCode, &
         iomsg = errorMsg)

    if(errorCode /= 0) call fatalError(Here, errorMsg)

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
        case default
          call fatalError(Here,'Unrecognised ACE CARD type: '// entry(i) % ZAID(last:last))
      end select

      ! Remove type information from the entry
      entry(i) % ZAID = entry(i) % ZAID(1:last-1)

      ! Add to map
      call map % add(entry(i) % ZAID, i)

    end do
  end subroutine load

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
