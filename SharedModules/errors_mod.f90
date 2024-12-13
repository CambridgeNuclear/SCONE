!!
!! Collection of functions to handle errors and warnings
!!
!! It declared as a singleton-like module with a state since it will contain
!! a warning log at a later date.
!!
module errors_mod
  use iso_fortran_env, only : error_unit

  use numPrecision
  use universalVariables, only : MAX_COL

  implicit none

contains

  !!
  !! Kill the execution and print the error message
  !!
  !! A pretty(ish) formatted error is printed to STDOUT. A backtrace should
  !! be produced as well (at least with gfortran but may be different for
  !! other compilers)
  !!
  !! Args:
  !!   where [in] -> Location of the error
  !!   why [in] -> Error message
  !!
  !! Note:
  !!   Error message will be broken at MAX_COL column. Breaking can happen
  !!   mid-word.
  !!
  subroutine fatalError(where, why)
    character(*), intent(in) :: where
    character(*), intent(in) :: why
    integer(shortInt)        :: end, start

    ! Upper frame
    write(error_unit, *) repeat('<>', MAX_COL / 2)
    write(error_unit, *) 'Fatal has occurred in:'
    write(error_unit, *) where, new_line('')
    write(error_unit, *) 'Because:'

    ! Do simple message folding
    ! Break on a whitespace if possible
    ! Does not recognise TABS, and other whitespace characters
    start = 1
    do while (start < len(why))
      end = min(start + MAX_COL, len(why))

      if (end < len(why)) then ! A line break is required
        end = start + scan(why(start: end), " ", back = .true.) - 1

        ! Should not happen often but be safe in case a string with no whitespace
        if (end <= start) end = min(start + MAX_COL, len(why))
      end if

      write(error_unit, '( " ",A)') why(start: end)
      start = end + 1
    end do

    ! Lower frame
    write(error_unit, *) repeat('<>', MAX_COL / 2)

    ! Terminate with backtrace
    ! NOTE: We assume MPI implementation will terminate all processes if one of them
    ! returns with an error code.
    error stop

  end subroutine fatalError


end module errors_mod
