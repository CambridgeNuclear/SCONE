module charLib_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dynArray_class,    only : dynIntArray

  implicit none
  private

  public :: splitChar


contains

  !!
  !! Split Character
  !!
  !! Loops once over the string and saves flip locations between region of only delimiters
  !! and not-delimiters.
  !!
  !! Args:
  !!   str [in]   -> input string
  !!   delim [in] -> 1 character delimiter
  !!
  !! Result:
  !!   2xN array with start and end index of separated substring
  !!
  function splitChar(str, delim) result (subLoc)
    character(*), intent(in)                       :: str
    character(1), intent(in)                       :: delim
    integer(shortInt), dimension(:,:), allocatable :: subLoc
    logical(defBool)                               :: inContent, notDelim
    integer(shortInt)                              :: i
    type(dynIntArray)                              :: flips
    character(100), parameter :: Here ='splitChar ( charLib_func.f90)'

    inContent = .false.

    do i =1,len(str)
      notDelim = str(i:i) /= delim

      ! Detect flip from sequence of delimiters to sequence of not delimiters
      if (inContent .neqv. notDelim) then
        ! Save flip location
        ! Need to take care to save i-1 if we are flippping from content to delimiters
        if( inContent) then
          call flips % add(i-1)
        else
          call flips % add(i)
        end if

        ! Flip state
        inContent = .not. inContent

      end if
    end do

    ! We need to save final location if we ended in content
    if(inContent) call flips % add(i-1)

    ! Transform flips to result array
    if( flips % getSize() == 0) then ! Empty string or only delimiters
      allocate( subLoc(2,1))
      subLoc = 0

    else if( mod(flips % getSize(), 2) /= 0) then ! Weird faliure
      call fatalError(Here,' Failed to obtain even number of flips. WTF?')

    else
      subLoc = reshape(flips % expose(), [2, flips % getSize()/2])

    end if
  end function splitChar


end module charLib_func
