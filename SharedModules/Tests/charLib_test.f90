module charLib_test

  use numPrecision
  use charLib_func, only : splitChar
  use pFUnit_mod

  implicit none


contains

  !!
  !! Test char splitting
  !!
@Test
  subroutine testSplitChar()
    character(:), allocatable :: line
    integer(shortInt), dimension(:,:),allocatable :: SE

    ! CASE 1: Standard Use Case
    line = " This is a first attempt at splitting"

    SE = splitChar(line, ' ')

    @assertEqual([2, 7, 10, 12, 18, 26, 29], SE(1,:))
    @assertEqual([5, 8, 10, 16, 24, 27, 37], SE(2,:))

    deallocate(line)
    deallocate(SE)

    ! CASE 2: Long delimiter sequences
    line ="...........12345......123..."
    SE = splitChar(line,'.')

    @assertEqual([12, 23], SE(1,:))
    @assertEqual([16, 25], SE(2,:))
    deallocate(line)
    deallocate(SE)

    ! CASE 3: Empty String
    line = ""
    SE = splitChar(line,'.')
    @assertEqual([0], SE(1,:))
    @assertEqual([0], SE(2,:))
    deallocate(line)
    deallocate(SE)

    ! CASE 4: Delimiter only string
    line ="#################################"
    SE = splitChar(line,'#')
    @assertEqual([0], SE(1,:))
    @assertEqual([0], SE(2,:))
    deallocate(line)
    deallocate(SE)


  end subroutine testSplitChar


    
end module charLib_test
