module genericProcedures
  use numPrecision
  implicit none

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

        if (errorStat > 0) call fatalError('openToRead subroutine (genericProcedures.f03)', &
                                           errorMsg )
  end subroutine openToRead




end module genericProcedures
