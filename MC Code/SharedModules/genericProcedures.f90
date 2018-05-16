module genericProcedures

  use numPrecision
  use endfConstants

  implicit none

  interface removeDuplicates
    module procedure removeDuplicates_Char
  end interface removeDuplicates

  interface linFind
    module procedure linFind_Char
    module procedure linFind_defReal
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

  interface numToChar
    module procedure numToChar_shortInt
    module procedure numToChar_longInt
    module procedure numToChar_defReal
  end interface


  integer(shortInt), parameter :: valueOutsideArray = -1,&
                                  tooManyIter       = -2,&
                                  targetNotFound    = -3

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


  pure function linearFloorIdxClosed_Real(array,value) result (idx)
    !! Performes linear search of an real sorted array and returns index of the largest
    !! element smaller-or-equal to the requested value. For the value equal to the largest element
    !! array(size(array)) it returns size(array)-1. For the value equal to the smallest element
    !! it returns 1. It returns -ve index in case of an error. Specific value is defined as a
    !! paramether. Following errors can happen
    !!   valueOutsideArray -> larger or smaller then array bounds
    real(defReal),dimension(:),intent(in) :: array
    real(defReal),intent(in)              :: value
    integer(shortInt)                     :: idx
    integer(shortInt)                     :: i

    if (value > array(size(array)) .or. value < array(1)) then
      idx = valueOutsideArray
      return
    end if

    do idx = size(array)-1,1,-1
      if( array(idx) <= value) return
    end do

  end function linearFloorIdxClosed_Real


  function linearFloorIdxClosed_shortInt(array,value) result(idx)
    !! Performes linear search of an integer sorted array and returns index of the largest element,
    !! which is smaller-or-equal to the requested value. Returns errors for emelents smaller and larger
    !! than the bounds of the array. For the value equal to the smallest element it returns 1 and
    !! for the value equal to the largest element it returns an error.
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

  pure function linearCeilingIdxOpen_shortInt(array,value) result(idx)
    !! Performes linear search of an integer sorted array and returns index of the smallest element,
    !! which is greater-or-equal to the requested value. Returns errors for elements larger than
    !! the upper bound of the array. Returns 1 for values smaller or equal to the lower bound of the
    !! array. Following errors can happen:
    !!   valueOutsideArray -> larger then the upper bound of array
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

  subroutine searchError(idx,Here)
    !! Subroutine that checks whether there was an error during search and returns approperiate
    !! message.
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

  end subroutine

  subroutine fatalError(Where,Why)
    character(*), intent(in)    :: Why, Where
    character(100)              :: Line, locWhy, locWhere
    character(20)               :: format
    integer(shortInt)           :: i

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

  subroutine openToRead(unitNum,File)
    integer(kind=shortInt), intent(in)    :: unitNum
    character(len=*), intent(in)          :: File
    integer(kind=shortInt)                :: errorStat
    character(len=99)                     :: errorMsg

        open ( unit   = unitNum,   &
               file   = File,      &
               status = "old",     &
               action = "read",    &
               iostat = errorStat, &
               iomsg  = errorMsg)

        !errorMsg=adjustR(errorMsg)

        if (errorStat > 0) call fatalError('openToRead subroutine (genericProcedures.f90)', &
                                           errorMsg )
  end subroutine openToRead

  function removeDuplicates_Char(charArray) result(out)
    !! Function that removes duplicates from input character array. It returns array of equal or
    !! smaller size. Unfortunatly Fortran requires output character to have specified length. Length
    !! of 100 is hardcoded at the moment. Function returns fatal error if input characters are of
    !! length greater then 100
    character(len=*),dimension(:),intent(in)       :: charArray
    character(len=100),dimension(:),allocatable    :: out
    logical(kind=defBool),dimension(:),allocatable :: unique
    integer(kind=shortInt)                         :: i,j

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
  end function

  function findDuplicates_Char(charArray) result(out)
    !! Function that finds duplicates in array of characters. Returns array that contains repeted
    !! element. Unfortunatly Fortran requires output character to have specified length. Length
    !! of 100 is hardcoded at the moment. Function returns fatal error if input characters are of
    !! length greater then 100
    character(len=*),dimension(:),intent(in)       :: charArray
    character(len=100),dimension(:),allocatable    :: out
    logical(kind=defBool),dimension(:),allocatable :: unique
    integer(kind=shortInt)                         :: i,j

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
  end function



  function linFind_Char(charArray,target) result(idx)
    !! Searches linearly for the occurance of target in charArray. Removes left blanks.
    !! Following Errors can occur:
    !! targetNotFound -> target is not present in the array
    character(*),dimension(:),intent(in) :: charArray
    character(*),intent(in)              :: target
    integer(shortInt)                    :: idx

    do idx=1,size(charArray)
     ! if( trim(charArray(idx)) == trim(target) ) return
      if( adjustl(charArray(idx)) == adjustl(target)) return
    end do
    idx = targetNotFound

  end function



  function linFind_defReal(defRealArray,target) result (idx)
    !! Searches linearly for the occurance of target in defRealArray. Following Errors can occur
    !! valueOutsideArray -> target is not present in the array
    real(defReal), dimension(:), intent(in)  :: defRealArray
    real(defReal), intent(in)                :: target
    integer(shortInt)                        :: idx

    do idx=1,size(defRealArray)
      if (defRealArray(idx) == target) return
    end do
    idx = targetNotFound
  end function linFind_defReal


  function linFind_shortInt(shortIntArray,target) result (idx)
    !! Searches linearly for the occurance of target in shortIntArray. Following Errors can occur
    !! valueOutsideArray -> target is not present in the array
    integer(shortInt), dimension(:), intent(in)  :: shortIntArray
    integer(shortInt), intent(in)                :: target
    integer(shortInt)                            :: idx

    do idx=1,size(shortIntArray)
      if (shortIntArray(idx) == target) return
    end do
    idx = targetNotFound
  end function linFind_shortInt

!  function arrayConcat(charArray) result(out)
!    !! Concatenate strings from an array into a single long character. Trims elements of char Array
!    !! and ads on blank between them for separation.
!    character(*),dimension(:),intent(in)       :: charArray
!    character(:),allocatable                   :: out
!    integer(shortInt)                          :: trimLen , i
!
!    ! Find total trim length of elements of charArray
!    trimLen=0
!    do i=1,size(charArray)
!      trimLen = trimLen + len(trim(charArray(i)))
!    end do
!
!    i = trimLen + size(charArray)
!    allocate(character(i) :: out)
!
!    out = ''
!
!    do i=1,size(charArray)
!      out = out // trim(charArray(i)) // ' '
!    end do
!
!  end function arrayConcat



  function arrayConcat(charArray) result(out)
    !! Concatenate strings from an array into a single long character (tape). Asjusts left and trims
    !! elements of char Array. Adds a blank at the end of a line
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



  function countSymbol(string,symbol) result(num)
    ! Function that searches counts all occurences of a "symbol" in a "string"
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


  function symbolBalance(str,leftS,rightS) result (balance)
    !! Goes through the string and adds +1 to balance for each leftS and -1 for each rightS. It
    !! terminates and returns -1 when balance becomes -ve.
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



  function indexOfNext(signs,string) result (idx)
    character(1),dimension(:),intent(in)     :: signs
    character(*), intent(in)                 :: string
    integer(shortInt)                        :: idx
    integer(shortInt)                        :: i
    integer(shortInt),dimension(size(signs)) :: temp_idx
    character(100),parameter                 :: here='indexOf (genericProcedures.f90)'

    temp_idx = index(string,signs)


    idx = minval(temp_idx,temp_idx > 0)

    if (idx == huge(temp_idx)) idx = 0

  end function indexOfNext


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



  elemental function RealReal_linlin_elemental_interpolate(xMin,xMax,yMin,yMax,x) result(y)
    real(defReal), intent(in) :: xMin, xMax, yMin, yMax, x
    real(defReal)             :: y
    real(defReal)             :: interFactor

    interFactor = (x-xMin)/(xMax-xMin)
    y = yMax * interFactor + (1-interFactor)*yMin

  end function RealReal_linlin_elemental_interpolate



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



  elemental function isInteger(float) result (isIt)
    !! Checks if the provided float is an integer. It may not be rebust and requires further
    !! testing. It uses the fact that ceiling and floor of an integer are the same.

    real(defReal), intent(in) :: float
    logical(defBool)          :: isIt
    integer(longInt)          :: tempI
    real(defReal)             :: a

    isIt = (floor(float,longInt) == ceiling(float,longInt))

  end function isInteger



  function isSorted_defReal(array) result (isIt)
    !! Function that check if the array is sorted in ascending order (a(i) >= a(i-1) for all i).
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



  function isSorted_shortInt(array) result (isIt)
    !! Function that check if the array is sorted in ascending order (a(i) >= a(i-1) for all i).
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

  !! Convert shortInt to character
  !! TODO: tempChar should have a parametrised length - need to come up with a smart way of doing it!
  function numToChar_shortInt(x) result(c)
    integer(shortInt)         :: x
    character(:), allocatable :: c
    character(40)             :: tempChar

    write(tempChar,'(I0)') x
    c = trim(tempChar)

  end function numToChar_shortInt

  !! Convert longInt to character
  !! TODO: tempChar should have a parametrised length - need to come up with a smart way of doing it!
  function numToChar_longInt(x) result(c)
    integer(longInt)          :: x
    character(:), allocatable :: c
    character(40)             :: tempChar

    write(tempChar,'(I0)') x
    c = trim(tempChar)

  end function numToChar_longInt

  !! Convert defReal to character
  !! TODO: tempChar should have a parametrised length - need to come up with a smart way of doing it!
  function numToChar_defReal(x) result(c)
    real(defReal)             :: x
    character(:), allocatable :: c
    character(40)             :: tempChar

    write(tempChar,'(F0.0)') x
    c = trim(tempChar)

  end function numToChar_defReal


  !!
  !! Subroutine takes a normilised direction vector dir and rotates it by cosine of a polar angle
  !! mu and azimuthal angle phi (in radians).
  !! Procedure will produce incorrect results without error message if dir is not normalised
  !!
  function rotateVector(dir,mu,phi) result(newDir)
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



end module genericProcedures
