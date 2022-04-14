module weightWindowsField_test
  use numPrecision
  use pFUnit_mod
  use particle_class,           only : particle
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use weightWindowsField_class, only : weightWindowsField

  implicit none


@testCase
  type, extends(TestCase) :: test_weightWindows
    private
    type(weightWindowsField) :: wwField
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_weightWindows


  !!
  !! Weight Windows Definition
  !!
  character(*), parameter :: DICT_DEF = &
  & "file /home/vr339/SCONE/scone/IntegrationTestFiles/testWW ;"

contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_weightWindows), intent(inout) :: this
    type(dictionary)                         :: dict

    call charToDict(dict, DICT_DEF)

    call this % wwField % init(dict)

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_weightWindows), intent(inout) :: this

    call this % wwField % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test retrieving the weight window values
  !!
@Test
  subroutine testGetValue(this)
    class(test_weightWindows), intent(inout) :: this
    type(particle)                           :: p
    integer(shortInt), dimension(3)          :: bins, EXPECTED_BINS

    p % isMG = .false.
    p % coords % lvl(1) % r = [0.5, 7.0, 0.0]
    p % E = 10.0

    bins = this % wwField % at(p)
    EXPECTED_BINS = [0.4, 1.5, 0.8]

    @assertEqual(EXPECTED_BINS,bins)

  end subroutine testGetValue


end module weightWindowsField_test
