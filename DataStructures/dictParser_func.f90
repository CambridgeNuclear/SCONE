module dictParser_func

  use numPrecision
  use genericProcedures,  only : fatalError, numToChar, replaceChar
  use charTape_class,     only : charTape
  use dictionary_class,   only : dictionary

  implicit none
  private

  !! Public procedures
  public :: fileToDict
  public :: charToDict


  ! Parameters
  integer(shortInt), parameter :: MAX_COLUMN = 1000
  character(2),dimension(2),parameter :: cmtSigns    = ['! ','//']
  integer(shortInt), parameter :: CONV_INT = 1, CONV_REAL = 2, CONV_CHAR = 3, CONV_UDEF = 0

  !!
  !! Psuedo dynamic type for reading entries
  !!
  !! Can read INT, REAL and CHARACTER
  !! Which one was read is indicated by type in {CONV_INT, CONV_REAL, CONV_CHAR}
  !!
  !! Public Members:
  !!   i -> Value of the integer if type == CONV_INT
  !!   r -> Value of the real if type == CONV_REAL
  !!   c -> Contents of character if type == CONV_CHAR
  !!   type -> Current type of the reader
  !!
  !! Interface:
  !!   convert -> read contents of the pathLen-long character into
  !!              i,r or c & set approperiate type
  !!
  type, private :: reader
    integer(shortInt)  :: i = 0
    real(defReal)      :: r = ZERO
    character(pathLen) :: c = ''
    integer(shortInt)  :: type = CONV_UDEF
  contains
    procedure :: convert => convert_reader
  end type reader


contains

  !!
  !! Reads contents of a file into dictionary
  !!
  !! Follows the SCONE dictionary grammar see Documentation for more detailed explenation.
  !!
  !! Informal definition of the grammar follows here. In general leading and trailing spaces
  !! are allowed for each <component>. When at least one space is REQUIRED it is marked with ' '
  !!
  !! <DICTIONARY> ::= <ITEM>+                           ! Dictionary consists of 1 or more ITEMS
  !! <ITEM>       ::= <ENTRY>|<KEYWORD>{<DICTIONARY>}
  !! <ENTRY>      ::= <KEYWORD>' '<CONTENT>;            ! At least 1 space in-between <keyword> & <content>
  !! <CONTENT>    ::= <SINGLE_CONTENT>|(<LIST_CONTENT>)
  !! <SINGLE_CONTENT> ::= Int|Real|Word
  !! <LIST_CONTENT> ::= <IntList>|<realList>|<charList>
  !! <IntList> ::= Space Separated Integers e.g.: 1 2   3
  !! <realList> ::= Space separated Numbers e.g.: 1.0 2   3.0E+2
  !! <charList> ::= Space separated not-numbers e.g.: 1.0\ word       char
  !! <KEYWORD>  ::= Character string with len_trim < nameLen constant
  !!
  !! Important facts to note:
  !!  1) Keywords are Case Sensitive keYwORD /= keyword
  !!  2) No words can contain whitespaces: kkkk -> OK; kk kk -> NOT OK
  !!  3) Format is free form but maximum column is set as a parameter to 1000
  !!  4) Only 32bit Integers can be read
  !!
  !! Args:
  !!   dict [inout] -> dictionary that will be filled with contents of a file. Can be both
  !!                   initialised or uninitialised
  !!   filePath [in] -> Path to the file that is to be read
  !!
  !! Errors:
  !!   fatalError if file under filePath does not exist
  !!   fatalError if file contains extra } symbols
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

    ! Reinitialise dictionary
    call dict % kill()
    call dict % init(5)

    ! Parse
    pos = 1
    call parseDict(dict, pos, file)

    ! Make sure there are no leftovers in the file
    ! Remember that pos will be 1 over last '}'
    if(pos-1 /= file % length()) then
      call fatalError(Here,"Not entire file was read. It means that there must be an &
                           &extra '}' bracket somwhere. ")
    end if

  end subroutine fileToDict

  !!
  !! Reads contents of a char string into dictionary
  !!
  !! Follows the SCONE dictionary grammar. See documentation and fileToDict for more details.
  !! String is interpreted as a single line! Even when it contains NEWLINE Characters!
  !! As a result no comments are allowed!
  !!
  !! Args:
  !!   dict [inout] -> dictionary that will be filled with contents of a string. Can be both
  !!                   initialised or uninitialised
  !!   data [in]    -> Character String that contains the data
  !!
  !! Errors:
  !!   fatalError if 'data' contains comment tokens ('//' & '!')
  !!   fatalError if 'data' contains extra } tokens
  !!
  subroutine charToDict(dict, data)
    class(dictionary), intent(inout) :: dict
    character(*), intent(in)         :: data
    character(len(data))             :: loc_data
    type(charTape)                   :: file
    integer(shortInt)                :: i, pos
    character(100),parameter :: Here = 'charToDict (dictParser_func.f90)'

    ! Check that data string has no comment
    do i=1,size(cmtSigns)
      if(index(data, trim(cmtSigns(i))) /= 0) then
        call fatalError(Here, 'Detected line comment sign: '//trim(cmtSigns(i))//&
                              ' Comments are nor allowes in dictionary made of character string')
      end if
    end do

    ! Remove Blank characters that can cause confusion
    loc_data = data
    call replaceChar(loc_data,char(9),' ')        ! Change tabs to single space
    call replaceChar(loc_data,new_line(' '), ' ') ! Replace New Line with space

    ! Create document charTape
    call file % append(loc_data)
    ! parceDict does not like '}}' ending so add an extra space
    call file % append(' }' )

    ! Reinitialise dictionary
    call dict % kill()
    call dict % init(5)

    ! Parse
    pos = 1
    call parseDict(dict, pos, file)

    ! Make sure there are no leftovers in the file
    ! Remember that pos will be 1 over last '}'
    if(pos-1 /= file % length()) then
      call fatalError(Here,"Not entire string was read. It means that there must be an &
                           &extra '}' bracket somwhere. ")
    end if

  end subroutine charToDict



  !!
  !! Parse contents of a tape starting from pos into dictionary
  !!
  !! Finish parsing when dict end symbol '}' is encountered
  !!
  !! This is recursive subroutine. Goes down a level whenever new subdictionary is encountered
  !!
  !! Args:
  !!   dict [inout] -> Initialised dictionary to load content into
  !!   pos [inout]  -> Starting position. Set to ending position on return
  !!   tape [in]    -> charTape with the contents to parse
  !!
  !! Errors:
  !!   fatalError if was given a tape without ';' '{' or '}' characters
  !!   fatalError if invalid syntax is encountered
  !!
  recursive subroutine parseDict(dict, pos, tape)
    class(dictionary), intent(inout) :: dict
    integer(shortInt), intent(inout) :: pos
    type(charTape), intent(in)       :: tape
    integer(shortInt)                :: nextSymbol, fin
    type(dictionary)                 :: tempDict
    character(nameLen)               :: name
    character(nameLen + 5)           :: buffer
    character(1)                     :: token
    character(100), parameter :: Here = 'parseDict (dictParser_func.f90)'

    ! Find next location & token
    do while (pos < tape % length())
      nextSymbol = tape % scanFrom(pos,';{}')
      if(nextSymbol == 0) call fatalError(Here,"There are tokens ';{}' in the charTape!'")
      fin = pos + nextSymbol - 1
      token = tape % get(fin)

      select case(token)
        case(';') ! Found an entry
          call readEntry(dict, pos, fin-1, tape)
          pos = fin + 1

        case('{') ! Beginning of subdictionary
          ! Get name and see if it fits into nameLen characters
          buffer = trim(adjustl(tape % segment(pos, fin - 1)))
          if( len_trim(buffer) <= nameLen) then
            name = trim(buffer)
          else
            call fatalError(Here, 'Subdictionary name: '//trim(buffer)//'... must fit into: '//&
                                  numToChar(nameLen)//' characters')
          end if

          ! Check that name does not contain any blanks
          if(scan(trim(name),' ') /= 0) then
            call fatalError(Here, 'Subdictionary name: '//trim(buffer)//' cannot contain spaces.')
          end if

          ! Prepare Subdictionary
          call tempDict % init(5)

          fin = fin + 1
          call parseDict(tempDict, fin, tape)
          call dict % store(name, tempDict)
          pos = fin
          ! Kill dictionary before the next use
          call tempDict % kill()

        case('}') ! End of dictionary
          ! Check that only blanks are present
          ! Must protect againt brackets with no content e.g. {}
          if (pos /= fin) then
            if (0 /= verify(tape % segment(pos, fin-1),' ')) then
              call fatalError(Here, "There is unparsed content. Missing a token ';' or '{...}' ")
            end if
          end if
          pos = fin + 1
          return

        case default
          call fatalError(Here,'token: '//token//" is not one of the tokens ';{}'. WTF!")

      end select
    end do

    ! Reach the end of file without terminal character
    call fatalError(Here, "End of file was reached without finding a termination symbol for&
                         & a subdictionary. It means that '}' must be missing somwhere.")


  end subroutine parseDict

  !!
  !! Read dictionary entry begining at 'start' and ending at 'end'
  !!
  !! Args:
  !!   dict [inout] -> Initialise dictionary to load content into
  !!   start [in]   -> Starting position on the tape
  !!   end [in]     -> Ending position on the tape
  !!   tape [in]    -> charTape with file contents
  !!
  !! Errors:
  !!   fatalError if start <= end
  !!   fatalError if parsing fails due to wrong input
  !!
  subroutine readEntry(dict, start, end, tape)
    class(dictionary), intent(inout) :: dict
    integer(shortInt), intent(in)    :: start
    integer(shortInt), intent(in)    :: end
    type(charTape), intent(in)       :: tape
    integer(shortInt)                :: p1, l1, l2, temp
    character(pathLen)               :: buffer
    character(nameLen)               :: name
    type(reader)                     :: converter
    character(100), parameter :: Here = 'readEntry (dictParser_func.f90)'

    ! Catch Invalid Entry
    if (end <= start) then
      call fatalError(Here,'End: '//numToChar(end)//' <= Start: '//numToChar(start)//&
                           ". Perhaps there is ;; or ; ; in the dictionary?")
    end if

    ! Separate name from content
    p1 = start
    buffer = readWord(p1, end, tape)
    if(len_trim(buffer) == 0) then
      call fatalError(Here,'There is an entry composed of blanks only. &
                           &Must have two ; with only spaces in-between. ' )

    elseif(len_trim(buffer) > nameLen) then
      call fatalError(Here,'Entry name: '//trim(buffer)//' is too long. Must fit into: '//&
                            numToChar(nameLen)//' characters.')
    end if
    name = trim(buffer)

    ! Verify that content exists
    if (p1 == end) call fatalError(Here,'Entry: '//trim(name)//' has no content!')

    ! Check if content is a list
    l1 = tape % scanFrom(p1, '(', end)
    l2 = tape % scanFrom(p1, ')', end)

    if(l1 == 0 .and. l2 == 0) then ! Read single entry
      buffer = readWord( p1, end, tape)
      ! Convert
      call converter % convert(buffer)

      ! Verify that there are no multiple entries
      if (p1 /= end) then
        buffer = readWord( p1, end, tape)
        if (len_trim(buffer) > 0) then
          call fatalError(Here,'Entry: '//trim(name)//' has unexpected content:'//trim(buffer))
        end if
      end if

      ! Read Entry
      select case(converter % type)
        case(CONV_INT)
          call dict % store(name, converter % i)

        case(CONV_REAL)
          call dict % store(name, converter % r)

        case(CONV_CHAR)
          call dict % store(name, converter % c)

        case default
          call fatalError(Here, 'Unrecognised entry type! WTF?  ')
      end select


    elseif(l1 == 0 .neqv. l2 == 0) then ! Only single bracket is present
      call fatalError(Here,'Entry: '//trim(name)//' has only one bracket. Two are needed for a list')

    elseif (l2 < l1) then ! Brackets are in wrong order
      call fatalError(Here,'Entry: '//trim(name)//' has wrong order of brackets: ) ( ')

    else ! Try to read a list
      l1 = p1 + l1 - 1
      l2 = p1 + l2 - 1
      call readList(dict, name, l1 + 1, l2 - 1, tape )

      ! Verify that there are no multipe entries
      ! After List
      if(l2 /= end) then
        l2 = l2 + 1
        buffer = readWord(l2, end, tape)
        if (len_trim(buffer) > 0) then
          call fatalError(Here,'Entry: '//trim(name)//' has unexpected content after list: '//trim(buffer))
        end if
      end if

      ! Before List
      if (l1 /= start) then
        temp = p1
        buffer = readWord(temp, l1 - 1, tape)
        if (len_trim(buffer) > 0) then
          call fatalError(Here,'Entry: '//trim(name)//' has unexpected content before list: '//trim(buffer))
        end if
      end if

    end if
  end subroutine readEntry

  !!
  !! Read List Content of an entry
  !!
  !! Valid lists are:
  !!   intList  -> contains only integers e.g. name (1 2 3);
  !!   realList -> contains only integers and reals e.g. name (1.0 2 3 1.0E-70);
  !!   charList -> contains only words (no numbers) e.g.:
  !!               name (1.0/ char word); -> IS OK!
  !!               name (1.0 char word);  -> IS NOT ALLOWED
  !!
  !! Args:
  !!   dict [inout] -> Initialise dictionary to load content into
  !!   name [in]    -> Name of the list entry
  !!   start [in]   -> Starting position on the tape
  !!   end [in]     -> Ending position on the tape
  !!   tape [in]    -> charTape with file contents
  !!
  !! Errors:
  !!   fatalError if input syntax is incorrect
  !!   fatalError if if list is not one of the valid types
  !!
  subroutine readList(dict, name, start, end, tape)
    class(dictionary), intent(inout) :: dict
    character(nameLen), intent(in)   :: name
    integer(shortInt), intent(in)    :: start
    integer(shortInt), intent(in)    :: end
    type(charTape), intent(in)       :: tape
    type(reader)                     :: converter
    integer(shortInt)                :: p1
    character(pathLen)               :: buffer
    integer(shortInt)                :: listType, N, i
    integer(shortInt), dimension(:), allocatable  :: temp_int
    real(defReal), dimension(:), allocatable      :: temp_real
    character(nameLen), dimension(:), allocatable :: temp_char
    character(100), parameter :: Here = 'readList (dictParser_func.f90)'

    if (end < start) then
      call fatalError(Here,'In: '//trim(name)//' End: '//numToChar(end)//' < Start: '//&
                           numToChar(start)//'. It means that () is present. Empty &
                           &lists are not allowed.' )
    end if

    ! Identify type of list
    p1 = start
    buffer = readWord(p1, end, tape)
    if (len_trim(buffer) == 0) call fatalError(Here,'In: '//trim(name)//' Empty lists are not allowed')
    call converter % convert(buffer)

    listType = converter % type
    N = 1
    do while (p1 < end)
      buffer = readWord(p1, end, tape)
      if (len_trim(buffer) == 0) exit  ! Empty string -> no more entries
      call converter % convert(buffer)

      ! Update state
      N = N +1
      select case(listType)
        case(CONV_INT)
          if(converter % type == CONV_REAL) then
            listType = CONV_REAL
          elseif(converter % type == CONV_CHAR) then
            call fatalError(Here,'In: '//trim(name)//" Mixed Numbers/Character lists are not allowed")
          end if

        case(CONV_REAL)
          if(converter % type == CONV_CHAR) then
            call fatalError(Here, 'In: '//trim(name)//" Mixed Numbers/Character lists are not allowed" )
          end if

        case(CONV_CHAR)
          if(converter % type /= CONV_CHAR) then
            call fatalError(Here, 'In: '//trim(name)//" Mixed Numbers/Character lists are not allowed" )
          end if

        case default
          call fatalError(Here,'In: '//trim(name)//' Unknown list type: '//numtoChar(listType)//'. WTF!')
      end select
    end do

    ! Read The list
    select case(listType)
      case(CONV_INT)
        allocate(temp_int(N))
        p1 = start
        do i=1,N
          call converter % convert( readWord(p1, end, tape))
          temp_int(i) = converter % i
          if(converter % type /= CONV_INT) then
            call fatalError(Here,'In: '//trim(name)//' WTF?  Not-Int in IntList!')
          end if
        end do
        call dict % store(name, temp_int)

      case(CONV_REAL)
        allocate(temp_real(N))
        p1 = start
        do i=1,N
          call converter % convert( readWord(p1, end, tape))
          if (converter % type == CONV_REAL) then
            temp_real(i) = converter % r
          elseif (converter % type == CONV_INT) then
            temp_real(i) = real(converter % i, defReal)
          else
            call fatalError(Here,'In: '//trim(name)//' WTF? Char in RealList!')
          end if
        end do
        call dict % store(name, temp_real)

      case(CONV_CHAR)
        allocate(temp_char(N))
        p1 = start
        do i=1,N
          call converter % convert( readWord(p1, end, tape))
          if (converter % type /= CONV_CHAR) call fatalError(Here ,'WTF? Not-char in CharList!')
          if(len_trim(converter % c) <= nameLen ) then
            temp_char(i) = trim(converter % c)
          else
            call fatalError(Here,'In: '//trim(name)//' CharList Entry: '//converter % c//' is to &
                                &long. Maximum length of characters in a charList is:'//&
                                 numToChar(nameLen))
          end if
        end do
        call dict % store(name, temp_char)

    end select
  end subroutine readList

  !!
  !! Read Word
  !!
  !! Given interval in a tape
  !!  -> skip blanks
  !!  -> Read Non Blank characters
  !!  -> Return when reaching new blank or interval end
  !!
  !! Args:
  !!   pos [inout] -> beginning of an interval. Set to first blank after word
  !!                  on exit or end if segment does not end with blank
  !!   end [in] -> end of an intervale
  !!   tape [in] -> charTape with data
  !!
  !! Returns:
  !!   Word without leading blanks. E.g.:
  !!   '    word  ' -> returns 'word' pos==9
  !!   ' a  bcd'    -> returns 'a' pos==3
  !!   '    end'    -> returns 'end' pos==7
  !!   '       '    -> returns '' pos==7
  !!   '1'          -> returns '1' pos==1
  !!   Return is a character(pathLen) with left-adjusted word
  !!
  !! Errors:
  !!   Returns blanks if end <= pos
  !!   Return blanks if segment contain only blanks
  !!
  function readWord(pos, end, tape) result(word)
    integer(shortInt), intent(inout) :: pos
    integer(shortInt), intent(in)    :: end
    type(charTape), intent(in)       :: tape
    character(pathLen)               :: word
    character(1)                     :: test
    integer(shortInt)                :: posEnd

    word = ''
    ! Process obvious errors
    if (end == pos .and. tape % get(pos) /= ' ') then
      word = tape % get(pos)
    elseif(end <= pos) then
      return
    end if

    ! Skip blanks
    test = ' '
    pos = pos - 1
    do while(test == ' ' .and. pos < end)
      pos = pos + 1
      test = tape % get(pos)
    end do

    ! Skip non-blanks
    posEnd = pos
    do while (test /= ' ' .and. posEnd < end)
      posEnd = posEnd + 1
      test = tape % get(posEnd)
    end do

    ! Update Position and returns string
    if (posEnd == end .and. test /= ' ') then ! Segment ends without a blank
      word = adjustl(tape % segment(pos, posEnd))
    elseif (pos == end) then ! Segment contains only blanks
      word = ''
    else ! Normal case there is word and ends on blank
       word = adjustl(tape % segment(pos, posEnd-1))
    end if
    pos = posEnd

  end function readWord

  !!
  !! Interpret contents of buffer
  !!
  !! Sets self % type to {CONV_INT, CONV_FLOAT, CONV_CHAR}
  !! loads value into approperite member {i, f, c}
  !!
  !! To identify entries uses 'read' Fortran data-transfer with following formats:
  !!   INT  -> '(I80)'
  !!   REAL -> '(ES.80.80)'
  !!   CHAR -> If both INT & REAL return error during reading
  !! The above assumes that pathLen=80
  !!
  !! Args:
  !!   buffer [in] -> pathLen long character with content to read
  !!
  !! Errors:
  !!   None
  !!
  subroutine convert_reader(self, buffer)
    class(reader), intent(inout)   :: self
    character(pathLen), intent(in) :: buffer
    integer(shortInt)              :: readErr
    character(100)                 :: fmt
    character(100),parameter :: Here = 'convert_reader (dictParser_func.f90)'

    ! Try to read data as a INTEGER and check if this gives an error
    fmt = "(I"//numToChar(pathLen)//")"
    read(unit = buffer, fmt = fmt , iostat = readErr) self % i

    if (readErr == 0) then
      self % type = CONV_INT
      return
    end if

    ! Try to read data as a REAL and check if this gives an error
    fmt = "(ES"//numToChar(pathLen)//"."//numToChar(pathLen)//")"
    read(unit = buffer, fmt = fmt , iostat = readErr) self % r

    if (readErr == 0) then
      self % type = CONV_REAL
      return
    end if

    ! Interpret as a character
    self % c = buffer
    self % type = CONV_CHAR

  end subroutine convert_reader

end module dictParser_func
