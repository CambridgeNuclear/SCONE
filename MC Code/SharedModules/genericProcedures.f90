module genericProcedures

  use numPrecision
  use endfConstants

  implicit none

  interface removeDuplicates
    module procedure removeDuplicates_Char
  end interface removeDuplicates

  interface linFind
    module procedure linFind_Char
  end interface

  interface findDuplicates
    module procedure findDuplicates_Char
  end interface

  interface binaryFloorIdx
    module procedure binaryFloorIdx_Real
  end interface

  interface linearFloorIdx
    module procedure linearFloorIdx_shortInt
  end interface

  interface endfInterpolate
    module procedure RealReal_endf_interpolate
  end interface

  interface interpolate
    module procedure RealReal_linlin_elemental_interpolate
  end interface

  contains

  function binaryFloorIdx_Real(array,value) result(idx)
    !! Performes binary search of an real sorted array and returns index of the largest element
    !! smaller-or-equal to the requested value. Returns error for elemets smaller and larger then
    !! the bounds of the array. For the value equal to the smallest element it returns 1 and for
    !! the value equal to the largest element it returns size(array).
    real(defReal),dimension(:),intent(in) :: array
    real(defReal),intent(in)              :: value
    integer(shortInt)                     :: idx
    integer(shortInt)                     :: bottom, top, i
    character(100),parameter              :: Here='binaryFloorIdx_Real (genericProcedures.f90)'

    ! Find Top and Bottom Index Array
    bottom = 1
    top = size(array)

    ! Check if the element is in array bounds
    if ( value <= array(bottom)) then
      if (value == array(bottom)) then ! Catch special case where value is equal to lower bound
        idx = bottom
        return
      else
        call fatalError(Here,'Requested Value is smaller than bottom element of array')
      end if
    else if ( value >= array(top)) then
      if (value == array(top)) then     ! Catch special case where value is equal to upper bound
        idx = top
        return
      else
        call fatalError(Here,'Requested Value is larger then top element of array')
      end if
    end if

    ! Perform search

    do i = 1,70
      !Calculate mid point
      idx = (top + bottom)/2

      ! Termination condition
      if (bottom == idx) return

      if (array(idx) < value ) then
        bottom = idx
      else if (array(idx) > value) then
        top = idx
      else ! When array(idx) == value
        return
      end if
    end do

    call fatalError(Here, 'Search loop failed to terminate after 70 iterations')

  end function binaryFloorIdx_Real

  function linearFloorIdx_shortInt(array,value) result(idx)
    !! Performes linear search of an integer sorted array and returns index of the largest element,
    !! which is smaller-or-equal to the requested value. Returns errors for emelents smaller and larger
    !! than the bounds of the array. For the value equal to the smallest element it returns 1 and
    !! for the value equal to the largest element it returns size(array)
    integer(shortInt),dimension(:),intent(in) :: Array
    integer(shortInt),intent(in)              :: Value
    integer(shortInt)                         :: idx
    character(100),parameter                  :: Here='linearFloorIdx_Int (genericProcedures.f90)'

    ! Check if the value is above the bounds of an array
    if ( Value > array(size(array))) call fatalError(Here,'Value is above upper bound of the array')

    do idx=size(array),1,-1
      if ( array(idx) <= value ) return
    end do

    call fatalError(Here,'Value is below lower bound of the array')

  end function linearFloorIdx_shortInt



  subroutine fatalError(Where,Why)
    character(*), intent(in)    :: Why, Where
    character(100)              :: Line, locWhy, locWhere
    character(20)               :: format
    integer(shortInt)           :: i

    Line = repeat('*',100)
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

  function linFind_Char(charArray,target) result(index)
    !! Searches linearly for the occurance of target in charArray. Returns index of -1 if target
    !! is not found. The index assumes that array begins at 1 (i.e. charArray(1:N)). If array begins
    !! with diffrent index (i.e. A(-5:N))the returned value needs to be approperiatly translated.
    character(len=*),dimension(:),intent(in) :: charArray
    character(len=*),intent(in)              :: target
    integer(kind=shortInt)                   :: index

    do index=1,size(charArray)
      if( trim(charArray(index)) == trim(target) ) return
    end do
    index = -1
  end function

  function arrayConcat(charArray) result(out)
    !! Concatenate strings from an array into a single long character. Trims elements of char Array
    !! and ads on blank between them for separation.
    character(*),dimension(:),intent(in)       :: charArray
    character(:),allocatable                   :: out
    integer(shortInt)                          :: trimLen , i

    ! Find total trim length of elements of charArray
    trimLen=0
    do i=1,size(charArray)
      trimLen = trimLen + len(trim(charArray(i)))
    end do

    allocate(character(trimLen+size(charArray)):: out)
    out = ''
    do i=1,size(charArray)
      out = out // trim(charArray(i)) // ' '
    end do
  end function arrayConcat

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

end module genericProcedures
