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
    character(len=100)              :: Line
    character(len=20)               :: format

    Line(1:100) = '*'
    format = 'A100'

    print format, Line
    print format, 'Fatal Error has occured in:'
    print format, Where
    print format, 'For the following reason:'
    print format, Why
    print format, Line
    stop
  end subroutine fatalError


end module genericProcedures
