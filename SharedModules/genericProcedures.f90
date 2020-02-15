module genericProcedures
  ! Intrinsic fortran Modules
  use iso_fortran_env, only : compiler_version

  use numPrecision
  use endfConstants
  use universalVariables

  implicit none

  interface swap
    module procedure swap_shortInt
    module procedure swap_defReal
    module procedure swap_char_nameLen
    module procedure swap_defReal_defReal
  end interface

  interface quickSort
    module procedure quickSort_shortInt
    module procedure quickSort_defReal
    module procedure quickSort_defReal_defReal
  end interface

  interface hasDuplicates
    module procedure hasDuplicates_shortInt
    module procedure hasDuplicates_defReal
  end interface

  interface removeDuplicates
    module procedure removeDuplicates_Char
    module procedure removeDuplicates_shortInt
  end interface removeDuplicates

  interface linFind
    module procedure linFind_Char
    module procedure linFind_defReal
    module procedure linFind_defReal_withTOL
    module procedure linFind_shortInt
  end interface

  interface findDuplicates
    module procedure findDuplicates_Char
  end interface

  interface binarySearch
    module procedure binaryFloorIdxClosed_Real
  end interface

  interface endfInterpolate
    module procedure RealReal_endf_interpolate
  end interface

  interface interpolate
    module procedure RealReal_linlin_elemental_interpolate
  end interface

  interface isSorted
    module procedure isSorted_defReal
    module procedure isSorted_shortInt
  end interface

  interface isDescending
    module procedure isDescending_defReal
    module procedure isDescending_shortInt
  end interface


  interface numToChar
    module procedure numToChar_shortInt
    module procedure numToChar_longInt
    module procedure numToChar_defReal
    module procedure numToChar_shortIntArray
    module procedure numToChar_defRealArray
  end interface

  contains

  pure function binaryFloorIdxClosed_Real(array,value) result(idx)
    !! Performes binary search of an real sorted array and returns index of the largest element
    !! smaller-or-equal to the requested value. For the value equalt to the largest element
    !! array(size(array)) it returns size(array)-1. For the value equal to the smallest element
    !! it returns 1. It returns -ve index in case of an error. Specific value is defined as a
    !! paramether. Following errors can happen
    !!   valueOutsideArray -> larger or smaller then array bounds
    !!   tooManyIter       -> algorithm did not convarged in required number of iterations
    real(defReal),dimension(:),intent(in) :: array
    real(defReal),intent(in)              :: value
    integer(shortInt)                     :: idx
    integer(shortInt)                     :: bottom, top, i

    ! Find Top and Bottom Index Array
    bottom = 1
    top = size(array)

    ! Check if the element is in array bounds
    if ( value < array(bottom) .or. value >array(top)) then
      idx = valueOutsideArray
      return
    end if

    do i = 1,70
      !Calculate mid point
      idx = (top + bottom)/2

      ! Termination condition
      if (bottom == idx) return

      ! Binary Step
      if (array(idx) <= value ) then
        bottom = idx
      else
        top = idx
      end if
    end do

    ! Failed to end in 70 steps
    idx = tooManyIter

  end function binaryFloorIdxClosed_Real

  !!
  !! Performes linear search of an real sorted array and returns index of the largest
  !! element smaller-or-equal to the requested value. For the value equal to the largest element
  !! array(size(array)) it returns size(array)-1. For the value equal to the smallest element
  !! it returns 1. It returns -ve index in case of an error. Specific value is defined as a
  !! paramether. Following errors can happen
  !!   valueOutsideArray -> larger or smaller then array bounds
  !!
  pure function linearFloorIdxClosed_Real(array,value) result (idx)
    real(defReal),dimension(:),intent(in) :: array
    real(defReal),intent(in)              :: value
    integer(shortInt)                     :: idx

    if (value > array(size(array)) .or. value < array(1)) then
      idx = valueOutsideArray
      return
    end if

    do idx = size(array)-1,1,-1
      if( array(idx) <= value) return
    end do

  end function linearFloorIdxClosed_Real

  !!
  !! Performes linear search of an integer sorted array and returns index of the largest element,
  !! which is smaller-or-equal to the requested value. Returns errors for emelents smaller and larger
  !! than the bounds of the array. For the value equal to the smallest element it returns 1 and
  !! for the value equal to the largest element it returns an error.
  !!
  function linearFloorIdxClosed_shortInt(array,value) result(idx)
    integer(shortInt),dimension(:),intent(in) :: Array
    integer(shortInt),intent(in)              :: Value
    integer(shortInt)                         :: idx
    character(100),parameter                  :: Here='linearFloorIdxClosed_shortInt (genericProcedures.f90)'

    ! Check if the value is above the bounds of an array
    if ( Value >= array(size(array))) call fatalError(Here,'Value is above upper bound of the array')

    do idx=size(array),1,-1
      if ( array(idx) <= value ) return
    end do

    call fatalError(Here,'Value is below lower bound of the array')

  end function linearFloorIdxClosed_shortInt

  !!
  !! Performes linear search of an integer sorted array and returns index of the smallest element,
  !! which is greater-or-equal to the requested value. Returns errors for elements larger than
  !! the upper bound of the array. Returns 1 for values smaller or equal to the lower bound of the
  !! array. Following errors can happen:
  !!   valueOutsideArray -> larger then the upper bound of array
  !!
  pure function linearCeilingIdxOpen_shortInt(array,value) result(idx)
    integer(shortInt),dimension(:),intent(in) :: Array
    integer(shortInt),intent(in)              :: Value
    integer(shortInt)                         :: idx
    character(100),parameter                  :: Here='linearCeilingIdxOpen_shortInt (genericProcedures.f90)'

    do idx=1,size(array)
      if ( array(idx) >= value ) return
    end do

    ! Value is larger than the upper bound of the array
    idx = valueOutsideArray

  end function linearCeilingIdxOpen_shortInt

  !!
  !! Subroutine that checks whether there was an error during search and returns approperiate
  !! message.
  !!
  subroutine searchError(idx,Here)
    integer(shortInt),intent(in)  :: idx
    character(*),intent(in)       :: Here

    if (idx < 0) then
      select case (idx)
        case (valueOutsideArray)
          call fatalError(Here,'The requested value was outide the array bounds')

        case (tooManyIter)
          call fatalError(Here,'Search did not terminate in hardcoded number of iterations')

        case (targetNotFound)
          call fatalError(Here,'Search failed to find target in the array')

        case default
          call fatalError(Here,'Search returned unknown error flag (negative index)')

      end select
    end if
  end subroutine searchError

  !!
  !! Displays error message and stops execution
  !!
  subroutine fatalError(Where,Why)
    character(*), intent(in)    :: Why, Where
    character(100)              :: Line, locWhy, locWhere
    character(20)               :: format

    Line = repeat('<>',50)
    format = '(A100)'
    locWhere = adjustR(where)
    locWhy = adjustR(why)

    print format, Line
    print format, 'Fatal Error has occured in:'
    print format, locWhere
    print *
    print format, 'For the following reason:'
    print format, locWhy
    print *
    print format, Line
    stop

  end subroutine fatalError

  !!
  !! Open "File" for reading under with "unitNum" reference
  !!
  subroutine openToRead(unitNum,File)
    integer(shortInt), intent(in)  :: unitNum
    character(*), intent(in)       :: File
    integer(shortInt)              :: errorStat
    character(99)                  :: errorMsg

        open ( unit   = unitNum,   &
               file   = File,      &
               status = "old",     &
               action = "read",    &
               iostat = errorStat, &
               iomsg  = errorMsg)

        if (errorStat > 0) call fatalError('openToRead subroutine (genericProcedures.f90)', &
                                           errorMsg )
  end subroutine openToRead


  !!
  !! Removes duplicates from an unsorted intArray
  !! Does not preserve the order
  !!
  pure function removeDuplicates_shortInt(intArray) result(out)
    integer(shortInt), dimension(:), intent(in) :: intArray
    integer(shortInt), dimension(:),allocatable :: out
    integer(shortInt),dimension(size(intArray)) :: array
    integer(shortInt)                           :: i, j, N, test, idx

    N = size(intArray)

    ! If provided array is of size 0, return size 0 array
    if (N == 0 ) then
      allocate(out(0))
      return
    end if

    array(1) = intArray(1)    ! Copy first element
    j = 2                     ! Set new empty index
    do i=2,N
      test = intArray(i)
      idx = linFind(intArray(1:i-1),test)  ! Check if element is present before i in the array
      if  (idx == targetNotFound ) then    ! If it isn't copy it to array without duplicates
        array(j) = test
        j = j + 1

      end if
    end do

    out = array(1:j-1)  ! Allocate output

  end function removeDuplicates_shortInt

  !!
  !! Function that removes duplicates from input character array. It returns array of equal or
  !! smaller size. Unfortunatly Fortran requires output character to have specified length. Length
  !! of 100 is hardcoded at the moment. Function returns fatal error if input characters are of
  !! length greater then 100
  !!
  function removeDuplicates_Char(charArray) result(out)
    character(nameLen),dimension(:),intent(in)   :: charArray
    character(nameLen),dimension(:),allocatable  :: out
    logical(defBool),dimension(:),allocatable    :: unique
    integer(shortInt)                            :: i,j

    if (len(charArray) > len(out)) call fatalError('removeDuplicates_Char (genericProcedures.f90)',&
                                                   'Maximum length of input character is 100 ')
    if (size(charArray) == 1) then
      out = charArray
    else
      allocate(unique(size(charArray)))
      unique = .true.
      ! For every element search if it matches any previous element. Change uniqe to false if it is
      ! repeted.
        do i = 1,size(charArray)
          search: &
          do j = 1, i-1
            if( trim(charArray(i)) == trim(charArray(j)) ) then
              unique(i) = .false.
              exit search
            end if
          end do search
        end do
        ! Select elements from charArray for which unique == .true.
        out = pack(charArray, unique)
    end if
  end function removeDuplicates_Char

  !!
  !! Function that finds duplicates in array of characters. Returns array that contains repeted
  !! element. Unfortunatly Fortran requires output character to have specified length. Length
  !! of 100 is hardcoded at the moment. Function returns fatal error if input characters are of
  !! length greater then 100
  !!
  function findDuplicates_Char(charArray) result(out)
    character(nameLen),dimension(:),intent(in)  :: charArray
    character(nameLen),dimension(:),allocatable :: out
    logical(defBool),dimension(:),allocatable   :: unique
    integer(shortInt)                           :: i,j

    if (len(charArray) > len(out)) call fatalError('removeDuplicates_Char (genericProcedures.f90)',&
                                                   'Maximum length of input character is 100 ')
    if (size(charArray) == 1) then
      out = charArray
    else
      allocate(unique(size(charArray)))
      unique = .true.
      ! For every element search if it matches any previous element. Change uniqe to false if it is
      ! repeted.
        do i = 1,size(charArray)
          search: &
          do j = 1, i-1
            if( trim(charArray(i)) == trim(charArray(j)) ) then
              unique(i) = .false.
              exit search
            end if
          end do search
        end do
        ! Select elements from charArray for which unique == .true.
        out = pack(charArray, .not.unique)
        out = removeDuplicates(out)
    end if
  end function findDuplicates_Char

  !!
  !! Searches linearly for the occurance of target in charArray. Removes left blanks.
  !! Following Errors can occur:
  !! targetNotFound -> target is not present in the array
  !!
  function linFind_Char(charArray,target) result(idx)
    character(*),dimension(:),intent(in) :: charArray
    character(*),intent(in)              :: target
    integer(shortInt)                    :: idx

    do idx=1,size(charArray)
     ! if( trim(charArray(idx)) == trim(target) ) return
      if( adjustl(charArray(idx)) == adjustl(target)) return
    end do
    idx = targetNotFound

  end function linFind_Char

  !!
  !! Searches linearly for the occurance of target in defRealArray. Following Errors can occur
  !! valueOutsideArray -> target is not present in the array
  !!
  function linFind_defReal(defRealArray,target) result (idx)
    real(defReal), dimension(:), intent(in)  :: defRealArray
    real(defReal), intent(in)                :: target
    integer(shortInt)                        :: idx

    do idx=1,size(defRealArray)
      if (defRealArray(idx) == target) return
    end do
    idx = targetNotFound
  end function linFind_defReal

  !!
  !! Searches linearly for the occurance of target in defRealArray.
  !!
  !! Uses relative tolerance of TOL
  !!
  !! Args:
  !!   defRealArray [in] -> array to search
  !!   target [in]       -> traget to search for
  !!   tol [in]          -> relative tolerance (+ve)
  !!
  function linFind_defReal_withTOL(defRealArray, target, tol) result (idx)
    real(defReal), dimension(:), intent(in)  :: defRealArray
    real(defReal), intent(in)                :: target
    real(defReal), intent(in)                :: tol
    integer(shortInt)                        :: idx
    real(defReal)                            :: inv_T

    inv_T = ONE/target

    do idx=1,size(defRealArray)
      if (abs(defRealArray(idx) * inv_T - ONE) < tol) return
    end do
    idx = targetNotFound
  end function linFind_defReal_withTOL

  !!
  !! Searches linearly for the occurance of target in shortIntArray. Following Errors can occur
  !! valueOutsideArray -> target is not present in the array
  !!
  pure function linFind_shortInt(shortIntArray,target) result (idx)
    integer(shortInt), dimension(:), intent(in)  :: shortIntArray
    integer(shortInt), intent(in)                :: target
    integer(shortInt)                            :: idx

    do idx=1,size(shortIntArray)
      if (shortIntArray(idx) == target) return
    end do
    idx = targetNotFound
  end function linFind_shortInt

  !!
  !! Concatenate strings from an array into a single long character (tape). Asjusts left and trims
  !! elements of char Array. Adds a blank at the end of a line
  !!
  function arrayConcat(charArray) result(out)
    character(*),dimension(:),intent(in)       :: charArray
    character(:),allocatable                   :: out
    integer(shortInt)                          :: elementLen, trimLen , i

    ! Find total length of elements of charArray after adjusting left and trimming
    trimLen=0
    do i=1,size(charArray)
      elementLen =  len( trim( adjustl( charArray(i))))
      trimLen = trimLen + elementLen
    end do

    ! Allocate output array
    i = trimLen+size(charArray)
    allocate(character(i) :: out)

    ! Make shure output tape is empty
    out = ''

    ! Write elements of the array to output tape
    do i=1,size(charArray)
      out = out // trim( adjustl(charArray(i))) // ' '
    end do

  end function arrayConcat

  !!
  !! Function that searches counts all occurences of a "symbol" in a "string"
  !!
  function countSymbol(string,symbol) result(num)
    character(*),intent(in)  :: string
    character(1),intent(in)  :: symbol
    integer(shortInt)        :: num
    integer(shortInt)        :: start, end, pos

    start = 1
    end   = len(string)
    num   = 0

    do
      pos = index(string(start:end),symbol)
      if (pos == 0) return ! No more occurences of symbol in "string"
      num = num + 1
      start = start + pos ! pos is returned relative to start
    end do

  end function countSymbol

  !!
  !! Goes through the string and adds +1 to balance for each leftS and -1 for each rightS. It
  !! terminates and returns -1 when balance becomes -ve.
  !!
  function symbolBalance(str,leftS,rightS) result (balance)
    character(*),intent(in)       :: str
    character(1),intent(in)       :: leftS
    character(1),intent(in)       :: rightS
    integer(shortInt)             :: balance
    integer(shortInt)             :: i

    balance = 0
    do i=1,len(str)

      if(str(i:i) == leftS) then
        balance = balance + 1

      elseif(str(i:i) == rightS) then
        balance = balance -1

      end if

      if (balance < 0) return

    end do

  end function symbolBalance

  !!
  !! Finds index of next character in signs in the string
  !!
  function indexOfNext(signs,string) result (idx)
    character(1),dimension(:),intent(in)     :: signs
    character(*), intent(in)                 :: string
    integer(shortInt)                        :: idx
    integer(shortInt),dimension(size(signs)) :: temp_idx
    character(100),parameter                 :: here='indexOf (genericProcedures.f90)'

    temp_idx = index(string,signs)


    idx = minval(temp_idx,temp_idx > 0)

    if (idx == huge(temp_idx)) idx = 0

  end function indexOfNext

  !!
  !! Replaces all repeated blanks in string with a single blank
  !!
  subroutine compressBlanks(string)
    character(*), intent(inout)     :: string
    character(len(string))          :: stringCopy
    integer(shortInt)               :: i, j
    logical(defBool)                :: lastBlank


    lastBlank = .false.
    stringCopy = ''
    j = 1

    do i=1,len(string)
      if (lastBlank) then
        if (string(i:i) /= " ") then
          lastBlank = .false.
          stringCopy(j:j) = string(i:i)
          j=j+1
        end if

      else
        stringCopy(j:j) = string(i:i)
        j=j+1
        if (string(i:i) == " ") then
          lastBlank = .true.
        endif

      end if
    end do
    string = stringCopy

  end subroutine compressBlanks

  !!
  !! Replaces all symbols "oldS" with "newS" in string
  !!
  subroutine replaceChar(string,oldS,newS)
    character(*), intent(inout) :: string
    character(1), intent(in)    :: oldS
    character(1), intent(in)    :: newS
    integer(shortInt)           :: i

    do i=1,len(string)
      if(string(i:i) == oldS) string(i:i) = newS
    end do

  end subroutine replaceChar

  !!
  !! Compares strings for equality. Ignores leading blanks.
  !!
  elemental function charCmp(char1, char2) result(areEqual)
    character(*), intent(in)  :: char1
    character(*), intent(in)  :: char2
    logical(defBool)          :: areEqual

    areEqual = (trim(adjustl(char1)) == trim(adjustl(char2)))

  end function charCmp

  !!
  !! Perform linear interpolation between defReals
  !!
  elemental function RealReal_linlin_elemental_interpolate(xMin,xMax,yMin,yMax,x) result(y)
    real(defReal), intent(in) :: xMin, xMax, yMin, yMax, x
    real(defReal)             :: y
    real(defReal)             :: interFactor

    interFactor = (x-xMin)/(xMax-xMin)
    y = yMax * interFactor + (1-interFactor)*yMin

  end function RealReal_linlin_elemental_interpolate

  !!
  !! Perform one of ENDF defined interpolation types
  !!
  function RealReal_endf_interpolate(xMin,xMax,yMin,yMax,x,endfNum) result(y)
    real(defReal), intent(in)     :: xMin, xMax, yMin, yMax, x
    integer(shortInt), intent(in) :: endfNum
    real(defReal)                 :: y
    character(100),parameter      :: Here='RealReal_endf_interpolate (genericProcedures.f90)'

    select case (endfNum) ! Naming Convention for ENDF interpolation (inY-inX) i.e. log-lin => logarithmic in y; linear in x
      case (histogramInterpolation)
        y = yMin

      case (linLinInterpolation)
        y = interpolate(xMin,xMax,yMin,yMax,x)

      case (linLogInterpolation)
        y = interpolate(log(xMin),log(xMax),yMin,yMax,log(x))

      case (logLinInterpolation)
        y = interpolate(xMin,xMax,log(yMin),log(yMax),x)
        y = exp(y)

      case (loglogInterpolation)
        y = interpolate(log(xMin),log(xMax),log(yMin),log(yMax),log(x))
        y = exp(y)

      case (chargedParticleInterpolation)
        ! Not implemented
        call fatalError(Here, 'ENDF interpolation law for charged Particles is not implemented')

      case default
        call fatalError(Here, 'Unknown ENDF interpolation number')
    end select

  end function RealReal_endf_interpolate

  !!
  !! Checks if the provided float is an integer. It may not be rebust and requires further
  !! testing. It uses the fact that ceiling and floor of an integer are the same.
  !!
  elemental function isInteger(float) result (isIt)
    real(defReal), intent(in) :: float
    logical(defBool)          :: isIt

    ! It is Fortran 2008 Standard Compliant. Really!
    isIt = (floor(float,longInt) == ceiling(float,longInt))

  end function isInteger

  !!
  !! Function that check if the array is sorted in ascending order (a(i) >= a(i-1) for all i).
  !!
  function isSorted_defReal(array) result (isIt)
    real(defReal),dimension(:),intent(in) :: array
    logical(defBool)                      :: isIt
    integer(shortInt)                     :: i

    do i=2,size(array)
      if (array(i) < array(i-1)) then
        isIt = .false.
        return
      end if
    end do

    isIt = .true.

  end function isSorted_defReal

  !!
  !! Function that check if the array is sorted in ascending order (a(i) >= a(i-1) for all i).
  !!
  function isSorted_shortInt(array) result (isIt)
    integer(shortInt),dimension(:),intent(in) :: array
    logical(defBool)                          :: isIt
    integer(shortInt)                         :: i

    do i=2,size(array)
      if (array(i) < array(i-1)) then
        isIt = .false.
        return
      end if
    end do

    isIt = .true.

  end function isSorted_shortInt

  !!
  !! Function that check if the array is sorted in descending order (a(i) <= a(i-1) for all i).
  !!
  function isDescending_defReal(array) result (isIt)
    real(defReal),dimension(:),intent(in) :: array
    logical(defBool)                      :: isIt
    integer(shortInt)                     :: i

    do i=2,size(array)
      if (array(i) > array(i-1)) then
        isIt = .false.
        return
      end if
    end do

    isIt = .true.

  end function isDescending_defReal

  !!
  !! Function that check if the array is sorted in descending order (a(i) <= a(i-1) for all i).
  !!
  function isDescending_shortInt(array) result (isIt)
    integer(shortInt),dimension(:),intent(in) :: array
    logical(defBool)                          :: isIt
    integer(shortInt)                         :: i

    do i=2,size(array)
      if (array(i) > array(i-1)) then
        isIt = .false.
        return
      end if
    end do

    isIt = .true.

  end function isDescending_shortInt

  !!
  !! Convert shortInt to character
  !! TODO: tempChar should have a parametrised length - need to come up with a smart way of doing it!
  !!
  function numToChar_shortInt(x) result(c)
    integer(shortInt),intent(in) :: x
    character(:), allocatable    :: c
    character(40)                :: tempChar

    write(tempChar,'(I0)') x
    c = trim(tempChar)

  end function numToChar_shortInt

  !!
  !! Convert shortInt array to character
  !!
  function numToChar_shortIntArray(x) result(c)
    integer(shortInt), dimension(:) ,intent(in) :: x
    character(:), allocatable                   :: c
    character(40)                               :: tempChar
    integer(shortInt)                           :: i

    c = ''
    if(size(x) == 0) then
      c = ''
    else
      write(tempChar,'(I0)') x(1)
      c = trim(tempChar)
    end if

    do i=2,size(x)
      write(tempChar,'(I0)') x(i)
      c = c//' '//trim(tempChar)
    end do

  end function numToChar_shortIntArray

  !!
  !! Convert longInt to character
  !! TODO: tempChar should have a parametrised length - need to come up with a smart way of doing it!
  !!
  function numToChar_longInt(x) result(c)
    integer(longInt),intent(in) :: x
    character(:), allocatable   :: c
    character(40)               :: tempChar

    write(tempChar,'(I0)') x
    c = trim(tempChar)

  end function numToChar_longInt

  !!
  !! Convert defReal to character
  !! TODO: tempChar should have a parametrised length - need to come up with a smart way of doing it!
  !!
  function numToChar_defReal(x) result(c)
    real(defReal),intent(in)  :: x
    character(:), allocatable :: c
    character(40)             :: tempChar

    write(tempChar,*) x
    c = trim(tempChar)

  end function numToChar_defReal

  !!
  !! Convert defReal array to character
  !!
  function numToChar_defRealArray(x) result(c)
    real(defReal), dimension(:),intent(in)  :: x
    character(:), allocatable               :: c
    character(40)                           :: tempChar
    integer(shortInt)                       :: i

    c = ''
    if(size(x) == 0) then
      c = ''
    else
      write(tempChar,*) x(1)
      c = trim(tempChar)
    end if

    do i=2,size(x)
      write(tempChar,*) x(i)
      c = c//' '//trim(tempChar)
    end do

  end function numToChar_defRealArray

  !!
  !! Convert character to an Integer
  !!
  !! Looks at first 20 places in the character and if they contain a valid integer
  !! it performs conversion.
  !!
  !! E.G
  !! '2' -> 2 OK!
  !! '-003' -> -3 OK!
  !! '  2E+03' -> FAIL!
  !! ' 7 is Swell' -> FAIL!
  !! '7                   A' -> 7 OK!
  !!
  !! Args:
  !!   str [in]      -> character to convert
  !!   error [inout] -> Optional. Set to .TRUE. if conversion has failed
  !! Result:
  !!   An integer contained within first 20 fields of str.
  !! Error:
  !!   Fortran Intrinsic error upon failed conversion Unless error is present.
  !!   If error argument is present it is set to .TRUE. if conversion failed by any reason and
  !!   the value of i becomes undefined
  !!
  function charToInt(str, error) result(i)
    character(*), intent(in)                 :: str
    logical(defBool),intent(inout), optional :: error
    integer(shortInt)                        :: i
    integer(shortInt)                        :: state

    if(present(error)) then
      read(str,'(I20)', IOSTAT=state) i
      error = state /= 0
    else
      read(str,'(I20)') i
    end if

  end function charToInt

  !!
  !! Subroutine takes a normilised direction vector dir and rotates it by cosine of a polar angle
  !! mu and azimuthal angle phi (in radians).
  !! Procedure will produce incorrect results WITHOUT error message if dir is not normalised
  !!
  pure function rotateVector(dir, mu, phi) result(newDir)
    real(defReal), dimension(3), intent(in)    :: dir
    real(defReal), intent(in)                  :: mu
    real(defReal), intent(in)                  :: phi
    real(defReal), dimension(3)                :: newDir
    real(defReal)                              :: u,v,w
    real(defReal)                              :: sinPol, cosPol, A, B

    ! Precalculate cosine and sine of polar angle
    sinPol = sin(phi)
    cosPol = cos(phi)

    ! Load Cartesian components of direction.
    u = dir(1)
    v = dir(2)
    w = dir(3)

    ! Perform standard roatation. Note that indexes are parametrised
    A = sqrt(max(ZERO, ONE - mu*mu))
    B = sqrt(max(ZERO, ONE - w*w  ))


    if ( B > 1E-8) then
      newDir(1) = mu * u + A * (u*w*cosPol - v * sinPol) / B
      newDir(2) = mu * v + A * (v*w*cosPol + u * sinPol) / B
      newDir(3) = mu * w - A * B * cosPol

    else
      B = sqrt(max(ZERO, ONE - v*v))
      newDir(1) = mu * u + A *(u*v*cosPol + w * sinPol) / B
      newDir(2) = mu * v - A * B * cosPol
      newDir(3) = mu * w + A* (v*w*cosPol - u * sinPol) /B

    end if

    newDir = newDir / norm2(newDir)

  end function rotateVector

  !!
  !! Dot product for 3D vector
  !!
  function dotProduct(a,b) result(x)
    real(defReal),dimension(3), intent(in) :: a,b
    real(defReal)                          :: x

    x = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

  end function dotProduct

  !!
  !! Cross product for 3D vectors
  !!
  function crossProduct(a,b) result(c)
    real(defReal),dimension(3),intent(in) :: a,b
    real(defReal),dimension(3)            :: c

    c = [a(3)*b(2) - a(2)*b(3), &
         a(1)*b(3) - a(3)*b(1), &
         a(2)*b(1) - a(1)*b(2) ]

  end function crossProduct

  !!
  !! Returns true if array contains duplicates
  !!
  pure function hasDuplicates_shortInt(array) result(doesIt)
    integer(shortInt), dimension(:), intent(in) :: array
    logical(defBool)                            :: doesIt
    integer(shortInt), dimension(size(array))   :: arrayCopy
    integer(shortInt)                           :: i

    ! Copy and sort array
    arrayCopy = array
    call quickSort(arrayCopy)

    ! Search through the array looking for duplicates
    doesIt = .false.
    do i=2,size(array)
      doesIt = doesIt .or. arrayCopy(i) == arrayCopy(i-1)

    end do

  end function hasDuplicates_shortInt

  !!
  !! Returns tue if array contains duplicates
  !!
  pure function hasDuplicates_defReal(array) result(doesIt)
    real(defReal), dimension(:), intent(in) :: array
    logical(defBool)                        :: doesIt
    real(defReal), dimension(size(array))   :: arrayCopy
    integer(shortInt)                       :: i

    ! Copy and sort array
    arrayCopy = array
    call quickSort(arrayCopy)

    ! Search through the array looking for duplicates
    doesIt = .false.
    do i=2,size(array)
      doesIt = doesIt .or. arrayCopy(i) == arrayCopy(i-1)

    end do

  end function hasDuplicates_defReal

  !!
  !! Quicksort for integer array
  !!
  recursive pure subroutine quickSort_shortInt(array)
    integer(shortInt), dimension(:), intent(inout) :: array
    integer(shortInt)                              :: pivot
    integer(shortInt)                              :: i, maxSmall

    if (size(array) > 1 ) then
      ! Set a pivot to the rightmost element
      pivot = size(array)

      ! Move all elements <= pivot to the LHS of the pivot
      ! Find position of the pivot in the array at the end (maxSmall)
      maxSmall = 0
      do i=1,size(array)

        if( array(i) <= array(pivot)) then
          maxSmall = maxSmall + 1
          call swap(array(i),array(maxSmall))
        end if
      end do

      ! Recursivly sort the sub arrays divided by the pivot
      call quickSort(array(1:maxSmall-1))
      call quickSort(array(maxSmall+1:size(array)))
    end if

  end subroutine quickSort_shortInt

  !!
  !! Quicksort for real array
  !!
  recursive pure subroutine quickSort_defReal(array)
    real(defReal), dimension(:), intent(inout) :: array
    integer(shortInt)                          :: i, maxSmall, pivot

    if (size(array) > 1 ) then
      ! Set a pivot to the rightmost element
      pivot = size(array)

      ! Move all elements <= pivot to the LHS of the pivot
      ! Find position of the pivot in the array at the end (maxSmall)
      maxSmall = 0
      do i=1,size(array)

        if( array(i) <= array(pivot)) then
          maxSmall = maxSmall + 1
          call swap(array(i),array(maxSmall))
        end if
      end do

      ! Recursivly sort the sub arrays divided by the pivot
      call quickSort(array(1:maxSmall-1))
      call quickSort(array(maxSmall+1:size(array)))
    end if

  end subroutine quickSort_defReal

  !!
  !! Quicksort for array1 and array2 of reals by array1
  !! If size(array2) > size(array1) ignores extra elements (does not sort them)
  !! If size(array1) > size(array2) behaviour is undefined.
  !!
  recursive subroutine quickSort_defReal_defReal(array1, array2)
    real(defReal), dimension(:), intent(inout) :: array1
    real(defReal), dimension(:), intent(inout) :: array2
    integer(shortInt)                          :: i, maxSmall, pivot
    character(100),parameter :: Here = 'quickSort_defReal_defReal (genericProcdures.f90)'

    if(size(array1) /= size(array2)) then
      call fatalError(Here,'Arrays have diffrent size!')
    end if

    if (size(array1) > 1 ) then
      ! Set a pivot to the rightmost element
      pivot = size(array1)

      ! Move all elements <= pivot to the LHS of the pivot
      ! Find position of the pivot in the array1 at the end (maxSmall)
      maxSmall = 0
      do i=1,size(array1)

        if( array1(i) <= array1(pivot)) then
          maxSmall = maxSmall + 1
          call swap(array1(i), array2(i), array1(maxSmall), array2(maxSmall))
        end if
      end do

      ! Recursivly sort the sub arrays divided by the pivot
      call quickSort(array1(1:maxSmall-1), array2(1:maxSmall-1))
      call quickSort(array1(maxSmall+1:size(array1)), array2(maxSmall+1:size(array1)))
    end if

  end subroutine quickSort_defReal_defReal

  !!
  !! Swap two integers
  !! Use XOR to avoid extra space (which is completly unnecessary)
  !!
  elemental subroutine swap_shortInt(i1,i2)
    integer(shortInt), intent(inout) :: i1
    integer(shortInt), intent(inout) :: i2
    integer(shortInt)                :: temp

    temp = i1
    i1 = i2
    i2 = temp

  end subroutine swap_shortInt

  !!
  !! Swap to reals
  !!
  elemental subroutine swap_defReal(r1,r2)
    real(defReal), intent(inout) :: r1
    real(defReal), intent(inout) :: r2
    real(defReal)                :: temp

    temp = r1
    r1 = r2
    r2 = temp

  end subroutine swap_defReal

  !!
  !! Swap two pair of reals
  !!
  elemental subroutine swap_defReal_defReal(r1_1, r1_2, r2_1, r2_2)
    real(defReal), intent(inout) :: r1_1
    real(defReal), intent(inout) :: r1_2
    real(defReal), intent(inout) :: r2_1
    real(defReal), intent(inout) :: r2_2
    real(defReal)                :: temp1, temp2

    ! Load first pair into temps
    temp1 = r1_1
    temp2 = r1_2

    ! Assign values of 2nd pair to 1st pair
    r1_1 = r2_1
    r1_2 = r2_2

    ! Assign values of 1st pair to 2nd pair
    r2_1 = temp1
    r2_2 = temp2

  end subroutine swap_defReal_defReal

  !!
  !! Swap character of length nameLen
  !!
  elemental subroutine swap_char_nameLen(c1,c2)
    character(nameLen), intent(inout) :: c1
    character(nameLen), intent(inout) :: c2
    character(nameLen)                :: temp

    temp = c1
    c1 = c2
    c2 = temp

  end subroutine swap_char_nameLen

  pure function charCshift(string,N) result(shifted)
    character(*), intent(in)      :: string
    integer(shortInt), intent(in) :: N
    character(len(string))        :: shifted
    integer(shortInt)             :: i, i_swap

    ! Loop over character and copy character with an offset
    do i=1,len(string)
      i_swap = modulo(i+N-1,len(string))+1
      shifted(i_swap:i_swap) = string(i:i)
    end do

  end function charCShift

  !!
  !! Convert Particle Type to string
  !!
  !! Args:
  !!   type [in] -> particle type
  !!
  !! Result:
  !!   Allocatable String that describes what particle is this
  !!
  !! Errors:
  !!   For unknown type prints "Unknown <int>" where int is number in type
  !!
  function printParticleType(type) result(str)
    integer(shortInt), intent(in) :: type
    character(:), allocatable     :: str

    select case(type)
      case(P_NEUTRON_CE)
        str = "CE Neutron"

      case(P_NEUTRON_MG)
        str = "MG Neutron"

      case default
        str = "Unknown "// numToChar(type)
    end select
  end function printParticleType


  !!
  !! Prints Scone ACII Header
  !!
  subroutine printStart()
    print *, repeat(" ><((((*> ",10)
    print *, ''
    print * ,"        _____ __________  _   ________  "
    print * ,"       / ___// ____/ __ \/ | / / ____/  "
    print * ,"       \__ \/ /   / / / /  |/ / __/     "
    print * ,"      ___/ / /___/ /_/ / /|  / /___     "
    print * ,"     /____/\____/\____/_/ |_/_____/     "
    print * , ''
    print * , ''
    print * , "Compiler Info :   ", compiler_version()
    print *
    print *, repeat(" <*((((>< ",10)

    ! TODO: Add extra info like date & time

  end subroutine printStart

  !!
  !! Prints line of fishes swiming right with an offset
  !!
  subroutine printFishLineR(offset)
    integer(shortInt),intent(in) :: offset
    integer(shortInt)            :: offset_L
    character(100), parameter    :: line = repeat(" ><((((*> ",10)
    character(100),dimension(10), parameter :: lines = [ &
    " ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*> " ,&
    "  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>" ,&
    ">  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*" ,&
    "*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((" ,&
    "(*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><(((" ,&
    "((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((" ,&
    "(((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><(" ,&
    "((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><" ,&
    "<((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  >" ,&
    "><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  ><((((*>  " ]

    offset_L = modulo(offset,10)

    print *, lines(offset_L+1)

  end subroutine  printFishLineR

end module genericProcedures
