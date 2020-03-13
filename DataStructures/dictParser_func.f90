module dictParser_func

  use numPrecision
  use genericProcedures,  only : fatalError, numToChar, replaceChar
  use charTape_class,     only : charTape
  use dictionary_class,   only : dictionary

  implicit none
  private

  !! Public procedures
  public :: fileToDict


  integer(shortInt), parameter :: MAX_COLUMN = 1000
  character(2),dimension(2),parameter :: cmtSigns    = ['! ','//']

contains

  !!
  !! Reads contents of a file into dictionary
  !!
  !! Follows the SCONE dictionary grammar
  !!  TODO: Define grammar and tell here where one can find it
  !!
  !! Args:
  !!   dict [inout] -> dictionary that will be filled with contents of a file. Can be both
  !!                   initialised or uninitialised
  !!   filePath [in] -> Path to the file that is to be read
  !!
  !! Errors:
  !!   fatalError if file under filePath does not exist
  !!
  subroutine fileToDict(dict, filePath)
    class(dictionary), intent(inout) :: dict
    character(*), intent(in)         :: filePath
    type(charTape)                   :: file
    integer(shortInt)                :: unit, stat, pos, i
    character(100)                   :: errorMsg
    character(:),allocatable         :: formatStr
    character(MAX_COLUMN)            :: buffer
    character(100),parameter :: Here = 'fileToDict (dictParser_func.f90)'

    ! Open file and read its contents to a charTape
    open ( newunit=unit, file=filePath, status="old", action="read", iostat=stat, iomsg=errorMsg)

    ! Process possible errors
    if (stat > 0) call fatalError(Here, errorMsg )

    ! Create format string for reading
    formatStr = '(A'//numToChar(MAX_COLUMN)//')'

    ! Read file contents
    stat = 0
    do
      read(unit=unit, fmt=formatStr, iostat=stat) buffer
      if(is_iostat_end(stat)) exit

      ! Remove all character after comment
      ! Loop over all comment signs
      do i=1,size(cmtSigns)
        pos = index(buffer, trim(cmtSigns(i)))

        ! If there is no comment pos = 0.
        ! If there is pos > 0. Then comment is removed
        if (pos > 0) buffer = buffer(1:pos-1)
      end do

      ! Replace all tabs with a space
      ! char(9) is TAB
      call replaceChar(buffer,char(9),' ')

      ! Append the file
      call file % append(trim(adjustl(buffer)))
      call file % append(' ')
    end do

    ! Close file
    close(unit)

    ! Append with dictionary terminator
    call file % append('}')

  end subroutine fileToDict



end module dictParser_func
