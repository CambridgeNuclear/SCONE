module genericProcedures
  use numPrecision
  implicit none

  interface removeDuplicates
    module procedure removeDuplicates_Char
  end interface removeDuplicates

  contains

  function binarySearch(Array,Value) result(index)
    real(kind=defReal),dimension(:),intent(in) :: Array
    real(kind=defReal),intent(in)              :: Value
    integer(kind=shortInt)                     :: index
    integer(kind=shortInt)                     :: Bottom, Top, i

    ! Find Top and Bottom Index Array
    Bottom = 1
    Top = size(Array)

    ! Perform Search
    Search: do i=1,10000
      if(Bottom +1 >= Top) then
        index = Bottom
        return
      end if

      index = (Bottom + Top)/2;

      if ( Array(index) <= Value ) then
        Bottom = index
      else
        Top = index
      end if
    end do Search

    ! Give Error if search failed to end
    call fatalError('Binary Search Function (genericProcedures.f03)',&
                    'Search loop failed to terminate after 10000 iterations')

  end function binarySearch


  subroutine fatalError(Where,Why)
    character(len=*), intent(in)    :: Why, Where
    character(len=100)              :: Line, locWhy, locWhere
    character(len=20)               :: format
    integer(kind=shortInt)           :: i

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

end module genericProcedures
