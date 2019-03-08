module IOdictionary_class

  use numPrecision
  use genericProcedures,  only : fatalError, arrayConcat, symbolBalance, indexOfNext, compressBlanks, &
                                 replaceChar
  use dictionary_class,   only : dictionary, charLen

  implicit none
  private
  !  Extension of dictionary class, which supports reading the dictionary from a text file.
  !  Format of the dictionary file follows the basic rules of an OpenFOAM dictionary, but not
  !  all features found at the website (https://cfd.direct/openfoam/user-guide/basic-file-format/)
  !  are supported. In principle following entries are supported:
  !
  !  KEYWORD PAIR :
  !  <keyword> <dataEntry>;
  !  ==>  pair of keyword and a Character, integer or real. Maximum size of keyword field is
  !       parameter nameLen in numPrecision module. Max size of data entry field i controled
  !       by parameter charLen.
  !
  !  1-D ARRAY (LIST) :
  !  <keyword> (<dataEntry1> <dataEntry2> <dataEntry3> ... );
  !  ==> keyword followed by "()" brackets enclosing a colections of space separated data entries.
  !      Array of type real, integer or character will be loaded depending on type of data Entries
  !      provided:
  !      -> if list contains at least one character entry, character array will be loaded
  !      -> if list contains numbers and at least one real, real array will be loaded
  !      -> if list is composed only of integers, integer array will be loaded
  !
  !  NESTED DICTIONARY:
  !  <keyword> { <keyword1> <dataEntry1>; <keyword2> (<dataEntry1> ...) ; <keyword3> {...} ... }
  !  ==> keyword fillowed by"{}" baracket enclosing colestions of KEYWORD PAIRS, LISTS or NESTED
  !      DICTIONARIES. Arbitrary number of dictionaries can be nested.
  !
  !  Notes:
  !  - Keywords are case sensitive keYwORD /= keyword
  !  - When entering characters (words) no "" or '' are needed :  NOT "word" BUT word
  !  - Input is generally free form but it cannot exceed column 500
  !
  ! TODO:
  !  -> create procedures for writing dictionary into a file for verification of input
  !

  ! PARAMTERES
  ! Note : maxColumn be an integer with 11 digits only. Should not be a problem for 32 bit.
  integer(shortInt),parameter         :: maxColumn   = 500
  integer(shortInt),parameter         :: numInt      = 1
  integer(shortInt),parameter         :: numReal     = 2
  integer(shortInt),parameter         :: word        = 3
  character(2),dimension(2),parameter :: cmtSigns    = ['! ','//']
  character(1),parameter              :: startOfDict = '{'
  character(1),parameter              :: startOfList = '('
  character(1),parameter              :: endOfEnt    = ';'
  character(1),dimension(3),parameter :: delim = [endOfEnt, startOfList, startOfDict]
  character(1),parameter              :: endOfList   = ')'
  character(1),parameter              :: endOfDict   = '}'


  type, public,extends(dictionary) :: IOdictionary
    private

  contains
    procedure :: initFrom
    final     :: final_IOdictionary
  end type IOdictionary

contains

  subroutine initFrom(self,filename)
    class(IOdictionary), intent(inout)             :: self
    character(*),intent(in)                        :: filename
    character(maxColumn),dimension(:),allocatable  :: file
    character(:),allocatable                       :: tape
    character(1)                                   :: tab
    character(100),parameter                       :: Here='initFrom (IOdictionary_class.f90)'

    file = readFile(filename)
    call removeComments(file)

    ! Crate a single long character (tape) from a file
    tape = arrayConcat(file)
    deallocate(file)

    ! Create tab character and replace any tabs with blanks
    !*** WORKS ON LINUX. BEWERE DIFFRENT TABS ON DIFFRENT OS
    tab = char(9)
    call replaceChar(tape,tab,' ')

    call compressBlanks(tape)

    call checkBrackets(tape,'(',')',Here,filename)
    call checkBrackets(tape,'{','}',Here,filename)

    call readDictionary(tape,self)

  end subroutine initFrom



  recursive subroutine readDictionary(tape,dict)
    character(*), intent(in)            :: tape
    class(dictionary), intent(inout)    :: dict
    character(len(tape))                :: localTape
    integer(shortInt)                   :: start, pos, listEnd, tapeEnd, dictEnd
    character(100),parameter            :: Here='readDictionary (IOdictionary_class.f90)'

    ! Kill dictionary
    call dict % kill()
    call dict % init(20,5)

    ! Remove leading blanks from input
    localTape = adjustl(tape)
    tapeEnd   = len(localTape)

    ! Find next delimiter ; or ( or { and read
    start = 1
    do while (len_trim(localTape) > 0)
      pos   = indexOfNext(delim,localTape(1:))
      if (pos == 1) then
        print *, localTape
        call fatalError(here,'Repeated delimiters or delim. at the beginning of a file  : ' // &
                              localTape(pos:len(localTape)))
      else if (pos == 0) then
        print *,localTape
        call fatalError(here,'There are no ";" delimiters at the end of entries in this dictionary')

      else if (localTape(pos:pos) == endOfEnt)  then
        ! Found KEYWORD PAIR
        call readEntry(dict,localTape(1:pos-1))
        start = pos +1

      else if (localTape(pos:pos) == startOfList) then
        ! Found LIST
        listEnd = index(localTape(1:),endOfEnt)
        if (listEnd == 0 ) call fatalError(Here, 'No ; after list : ' // localTape)
        call readList(dict,localTape(1:listEnd-1))
        start = listEnd +1

      else if (localTape(pos:pos) == startOfDict) then
        ! Found DICTIONARY
        dictEnd = findDictEnd(localTape(pos:)) + pos
        if (dictEnd == 0) call fatalError(Here, '} to end nested dictionary not found.')
        call readNestedDict(dict,localTape(1:dictEnd))
        start = dictEnd +1

      else
        print *, localTape, pos
        call fatalError(here,'Impossible error. Tape reading stoped at undefined delimiter')

      end if
      ! "Cut" the read entry from the beginning of the tape
      localTape = adjustl(localTape(start:))

    end do

  end subroutine readDictionary



  function findDictEnd(tape) result( pos )
    !! Function to find end of top level dictionary, which ignores '}' from nested dictionaries
    !! Tape need to start with '{'. Goes linearly through the string and tracks nesting level.
    character(*),intent(in)  :: tape
    integer(shortInt)        :: pos
    integer(shortInt)        :: level
    character(100),parameter :: Here='findDictEnd (IOdictionary_class.f90)'

    if (tape(1:1) /= startOfDict ) then
      call fatalError(Here,'Input tape does not start with '//startOfDict// ' : '//tape)
    end if

    level = 0
    do pos=1,len(tape)
      if (tape(pos:pos) == startOfDict) level = level -1
      if (tape(pos:pos) == endOfDict)   level = level +1
      if (level == 0) exit

    end do

  end function findDictEnd



  recursive subroutine readNestedDict(dict,entry)
    !! Reads keyword of the DICTIONARY entry and calls readDictionary to read its contents.
    class(dictionary), intent(inout) :: dict
    character(*), intent(in)         :: entry
    type(dictionary)                 :: localDict
    character(:),allocatable         :: tempKeyword
    character(:),allocatable         :: tempData
    character(nameLen)               :: keyword
    integer(shortInt)                :: lastC
    character(100),parameter         :: Here='readNestedDict (IOdictionary_class.f90)'

    call readWord(entry,tempKeyword,tempData)
    lastC = len_trim(tempData)

    ! Perform checks
    if (len(trim(tempKeyword)) > nameLen) then
      call fatalError(Here,'Keyword name is longer than allowed {nameLen-parameter}: ' // entry)

    else if (tempData(lastC:lastC) /= endOfDict ) then
    print *, entry
      call fatalError(Here,'Garbage detected after the end bracket of the dict. : ' // entry)

    else if (tempData(1:1) /= startOfDict ) then
      call fatalError(Here,'Garbage before the dictionary bracket : ' // entry)

    end if

    ! Set keyword
    keyword = adjustl(tempKeyword)

    ! Remove brackets and repeted blanks. Align to the left.
    tempData(1:1)         = ' '
    tempData(lastC:lastC) = ' '
    call compressBlanks(tempData)
    tempData = adjustl(tempData)

    ! Read and load the nested dictionary into higher level dictionary
    call readDictionary(tempData,localDict)
    call dict % store(keyword,localDict)

  end subroutine readNestedDict



  subroutine readEntry(dict,entry)
    !! Read and load KEYWORD PAIR
    class(dictionary), intent(inout)  :: dict
    character(*), intent(in)          :: entry
    integer(shortInt)                 :: dataType, dummyInt
    real(defReal)                     :: dummyReal
    character(nameLen)                :: keyword
    character(:),allocatable          :: tempData, tempKeyword
    character(charLen)                :: data
    character(100),parameter          :: Here='readEntry (IOdictionary_class.f90)'
    character(100)                    :: format

    call readWord(entry,tempKeyword,tempData)

    ! Preform error checks
    if (index(trim(tempData)," ") /= 0) then
       call fatalError(Here,'There are multiple data entries in : ' // entry)

    else if (len(trim(tempData)) > charLen) then
      call fatalError(Here,'The entery is longer than allowed {nameLen-parameter}: ' // entry)

    else if (len(trim(tempKeyword)) > nameLen) then
      call fatalError(Here,'Keyword name is longer than allowed {nameLen-parameter}: ' // entry)

    end if

    ! Store temporary data
    data    = adjustl(tempData)
    keyword = adjustl(tempKeyword)

    ! Identify and load data into dictionary after creating approperiate read format
    dataType = identifyData(data)

    select case (dataType)
      case(numInt)
        write(format,"(A2,I10.1,A2)") "(I",charLen,")"
        read(data,format) dummyInt
        call dict % store(keyword,dummyInt)

      case(numReal)
        write(format,"(A2,I10.1,A1,I10.1,A2)") "(ES",charLen,".",charLen,")"
        read(data,format) dummyReal
        call dict % store(keyword,dummyReal)

      case(word)
        call dict % store(keyword,data)

      case default
        call fatalError(Here,'Unrecognised data type.')

    end select

  end subroutine readEntry



  subroutine readWord(tape,word,tail)
    !! Function that reads first ' ' delimited entry of a tape and returns the rest of the tape as
    !! tail.
    character(*),intent(in)                :: tape
    character(len(tape))                   :: localTape
    character(:), allocatable,intent(out)  :: word
    character(:), allocatable,intent(out)  :: tail
    integer(shortInt)                      :: blankPos
    character(100),parameter               :: Here='readWord (IOdictionary_class.f90)'

    ! Preprocess tape and find position of first ' '
    localTape = adjustl(tape)
    call compressBlanks(localTape)
    blankPos = index(trim(localTape)," ")

    ! Cut tape into keyword and data after performing checks
    if(blankPos == 0 ) then
      call fatalError(Here,'Missing data for keyword: ' // localTape)

    else if (blankPos == 1) then
      call fatalError(Here,'Impossible error. Local tape begins with blank : ' // localTape)

    else
      word = localTape(1:blankPos-1)
      tail = trim( localTape(blankPos+1:len(localTape)) )

    end if

  end subroutine readWord



  subroutine readList(dict,list)
    class(dictionary), intent(inout)            :: dict
    character(*), intent(in)                    :: list
    character(:),allocatable                    :: tempKeyword, tempData
    character(nameLen)                          :: keyword
    integer(shortInt)                           :: lastC, datSize, i, blankPos, startPos
    character(100),parameter                    :: Here='readList (IOdictionary_class.f90)'
    integer(shortInt),dimension(:),allocatable  :: dataTypes
    character(charLen),dimension(:),allocatable :: listEntries
    logical(defBool)                            :: isCharacter, isReal, isInteger


    ! Read keyword and find last non-blank character index
    call readWord(list, tempKeyword, tempData)
    lastC = len_trim(tempData)

    ! Perform checks
    if (len(trim(tempKeyword)) > nameLen) then
      call fatalError(Here,'Keyword name is longer than allowed {nameLen-parameter}: ' // list)

    else if (tempData(lastC:lastC) /= endOfList ) then
      call fatalError(Here,'Garbage detected after the end bracket of the list : ' // list)

    else if (tempData(1:1) /= startOfList) then
      call fatalError(Here,'Garbage detected after keyword in a list : ' // list)

    end if

    ! Set keyword
    keyword = adjustl(tempKeyword)

    ! Remove all brackets and repeted blanks. Align to the left.
    ! Brackets of the nested (illigal) list will be removed as well,
    ! but reading should have failed before on delimiter serach it
    ! such nested list exists.
    call replaceChar(tempData,startOfList,' ')
    call replaceChar(tempData,endOfList,' ')
    call compressBlanks(tempData)
    tempData = adjustl(tempData)

    ! Find length of the list
    datSize = listLen(tempData)

    allocate( dataTypes(datSize)  )
    allocate( listEntries(datSize))

    startPos = 1

    ! Separate entries into individual char array components. It is crucial that tempData does not
    ! start with a blank " "! Should not becouse of the "adjustl" ealier.
    do i=1,datSize
      blankPos = index(tempData(startPos:)," ")

      if (blankPos -1 > charLen) then
        call fatalError(Here,'The list entery is longer than allowed {charLen-parameter}: ' // list)

      end if

      listEntries(i) = tempData(startPos:startPos+blankPos-1)
      startPos = startPos + blankPos
    end do

    ! Indentity list entries
    do i=1,datSize
      dataTypes(i) = identifyData( listEntries(i))
    end do

    ! Select array type
    isCharacter = all( dataTypes == word )
    isReal    = all((dataTypes == numReal) .or. (dataTypes == numInt)) .and. all( dataTypes /= word)
    isInteger = all( dataTypes == numInt)

    ! Load list depending on the type
    if (isCharacter) then
      call dict % store(keyword,listEntries)

    elseif(isInteger) then
      call readIntArray()   ! Internal subroutine

    elseif(isReal .and. (.not.isInteger)) then
      call readRealArray()  ! Internal subroutine

    else
      call fatalError(Here,'Mixed type list:'// tempData )

    end if

    contains

      subroutine readIntArray()
        !! A MACRO for reading a integer array to help with code clarity
        integer(shortInt), dimension(:), allocatable :: tempInt
        character(100)                               :: format

        allocate(tempInt(datSize))
        write(format,"(A2,I10.1,A2)") "(I",charLen,")"

        do i=1,datSize
          read(listEntries(i),format) tempInt(i)
        end do

        call dict % store(keyword,tempInt)

      end subroutine readIntArray


      subroutine readRealArray()
        !! A MACRO for reading real Array to help with code clarity
        real(defReal), dimension(:), allocatable :: tempReal
        character(100)                           :: formatInt
        character(100)                           :: formatReal
        integer(shortInt)                        :: tempInt


        allocate(tempReal(datSize))
        write(formatInt,"(A2,I10.1,A2)") "(I",charLen,")"
        write(formatReal,"(A2,I10.1,A1,I10.1,A2)") "(ES",charLen,".",charLen,")"

        do i=1,datSize
          select case (dataTypes(i))
            case(numInt)
              ! If integer i.e. 7 is read as a real there will be no
              ! read error but the read number will be wrong i.e. 7.0E-30
              read(listEntries(i),formatInt) tempInt
              tempReal(i) = real(tempInt,defReal)

            case(numReal)
              read(listEntries(i),formatReal) tempReal(i)

            case default
              call fatalError(Here,'Unrecoginised data type')

          end select
        end do

        call dict % store(keyword, tempReal)

      end subroutine readRealArray

  end subroutine readList



  recursive function listLen(data) result(length)
    !! Finds length of blank deliited list of entries. "data" must NOT begin with a blank.
    character(*),intent(in)  :: data
    integer(shortInt)        :: length
    integer(shortInt)        :: pos, end
    character(100),parameter :: Here='listLen (IOdictionary_class.f90)'

    if (data(1:1) == " ") call fatalError(Here, 'Input "data" begins with a blank')

    pos = index(trim(data)," ")
    end = len(data)

    if (pos == 0) then
      if(len_trim(data) == 0) then
        length = 0
      else
        length = 1
      end if
    else
      length = 1 + listLen(data(pos+1:end))
    end if

  end function listLen



  function identifyData(data) result(dataType)
    !! Identify type of the data stored as a character
    character(*), intent(in) :: data
    integer(shortInt)        :: dataType
    integer(shortInt),save   :: noError        = 0
    logical(defBool),save    :: firstTime      = .true.
    integer(shortInt)        :: i, readErr
    real(defReal)            :: r
    character(5)             :: dummy
    character(100)           :: format
    character(100),parameter :: Here='identifyData (IOdictionary_class.f90)'

    ! On first execution generate no error code for local system (should be 0).
    ! This is propably redundant but better safe then sorry.
    if (firstTime) then
      dummy='17.01'
      read (unit = dummy, fmt ="(F5.5)" , iostat = noError )  r
      firstTime = .false.

    end if

    ! Try to read data as a INTEGER and check if this gives an error
    write(format,"(A2,I10.1,A2)") "(I",charLen,")"
    read(unit = data, fmt = format , iostat = readErr) i

    if (readErr == noError) then
      dataType = numInt
      return
    end if

    ! Try to read data as a REAL and check if this gives an error
    write(format,"(A2,I10.1,A1,I10.1,A2)") "(ES",charLen,".",charLen,")"
    read(unit = data, fmt = format , iostat = readErr) r

    if (readErr == noError) then
      dataType = numReal
      return
    end if

    ! Try to read data as a CHARACTER and check if this gives an error
    write(format,"(A2,I10.1,A2)") "(A",charLen,")"
    read(unit = data, fmt = format , iostat = readErr) dummy

    if (readErr == noError) then
      dataType = word
      return
    end if

    ! Identification failed for some reason. Give error
    call fatalError(Here,'Procedure failed to identify data entery type. Sorry...')

    ! AVoid compier warning
    dataType = huge(word)
  end function identifyData



  subroutine checkBrackets(tape,leftB,rightB,Where,filePath)
    !! Checks consistency in number and sequence of brackets [catches ')(']
    character(*), intent(in) :: tape
    character(*), intent(in) :: Where
    character(*), intent(in) :: filePath
    character(1), intent(in) :: leftB
    character(1), intent(in) :: rightB
    integer(shortInt)        :: balance

    balance = symbolBalance(tape,leftB,rightB)

    if (balance > 0) then
      call fatalError(Where,'Missing '//rightB//' bracket in: ' // trim(filePath))

    elseif (balance < 0) then
      call fatalError(Where,'Missing '//leftB//' bracket in: ' // trim(filePath))
    end if

  end subroutine checkBrackets



  subroutine removeComments(file)
    !! Removes line comments
    character(*),dimension(:),intent(inout) :: file
    integer(shortInt)                       :: i,j, pos

    ! Loop over all lines
    do i=1,size(file)

      ! Loop over all comment signs
      do j=1,size(cmtSigns)

        pos = index(file(i), trim(cmtSigns(j)))

        ! If there is no comment pos = 0.
        ! If there is pos > 0. Then comment is removed
        if (pos > 0) then
          file(i)=file(i)(1:pos-1)
        end if
      end do
    end do

  end subroutine removeComments



  function readFile(filename) result(lines)
    !! Reads file into a character array for further processing
    character(*), intent(in)                      :: filename
    character(maxColumn),allocatable,dimension(:) :: lines
    integer(shortInt)                             :: errorStat
    character(100)                                :: errorMsg
    integer(shortInt),parameter                   :: fileUnit = 72
    integer(shortInt)                             :: fileLen,i
    integer(shortInt)                             :: readStat=0
    character(maxColumn)                          :: line, pastLine
    character(100)                                :: format

    ! Open the given file
    open ( unit   = fileUnit,   &
           file   = fileName,      &
           status = "old",     &
           action = "read",    &
           iostat = errorStat, &
           iomsg  = errorMsg)

    if (errorStat > 0) call fatalError('initFrom (IOdictionary.f90)', errorMsg )

    ! Find number of lines in the input file
    fileLen=0
    do
      write(format,"(A1,I10.1,A2)") "(",maxColumn,"A)"

      read(unit = fileUnit, fmt=format,iostat=readStat) line

      if(readStat == -1 .and. line==pastLine ) exit ! See character equality to indicate double read of last line
      fileLen = fileLen + 1
      pastLine = line
    end do
    rewind(fileUnit)

    ! Read file into character array

    allocate( lines(fileLen))

    do i=1,fileLen
      write(format,"(A1,I10.1,A2)") "(",maxColumn,"A)"
      read(unit = fileUnit, fmt=format,iostat=readStat) lines(i)
    end do

    close(fileUnit)

  end function readFile

  !!
  !! Finalisation subroutine for IOdictionary
  !!
  subroutine final_IOdictionary(self)
    type(IOdictionary), intent(inout) :: self

    call self % kill()

  end subroutine

    
end module IOdictionary_class
