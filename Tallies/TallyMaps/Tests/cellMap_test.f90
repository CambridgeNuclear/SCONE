module cellMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,      only : particleState
  use dictionary_class,    only : dictionary
  use dictParser_func,     only : charToDict
  use outputFile_class,    only : outputFile
  use cellMap_class,       only : cellMap

  implicit none


@testCase
  type, extends(TestCase) :: test_cellMap
    private
    type(cellMap),allocatable :: map_noUndef
    type(cellMap),allocatable :: map_Undef
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_cellMap


  !!
  !! Test parameters
  !!
  integer(shortInt),dimension(*),parameter :: CELL_NAMES = [1, 2, 3, 4, 5]
  integer(shortInt),dimension(*),parameter :: CELL_IN_MAP = [2, 3, 5]

contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_cellMap), intent(inout) :: this
    integer(shortInt)                  :: temp
    integer(shortInt)                  :: i
    type(dictionary)                   :: dict
    type(dictionary)                   :: mapDict1

    ! Initialise dictionaries
    call mapDict1 % init(2)

    ! Build material map definition
    call mapDict1 % store('cells', CELL_IN_MAP)
    allocate(this % map_noUndef, source = cellMap(mapDict1))

    call mapDict1 % store('undefBin','true')
    allocate(this % map_undef, source = cellMap(mapDict1))

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_cellMap), intent(inout) :: this

    call this % map_noUndef % kill()
    call this % map_Undef % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Mapping test without undefined bin
  !!
@Test
  subroutine testMappingNoUndefined(this)
    class(test_cellMap), intent(inout)       :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(5)           :: bins
    integer(shortInt),dimension(5),parameter :: EXPECTED_BINS = [0, 1, 2, 0, 3]

    do i=1,5
      state % cellIdx = i
      bins(i) = this % map_noUndef % map(state)
    end do

    @assertEqual(EXPECTED_BINS,bins)

  end subroutine testMappingNoUndefined


  !!
  !! Mapping test with undefined bin
  !!
@Test
  subroutine testMappingUndefined(this)
    class(test_cellMap), intent(inout)       :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(5)           :: bins
    integer(shortInt),dimension(5),parameter :: EXPECTED_BINS = [4, 1, 2, 4, 3]

    do i=1,5
      state % cellIdx = i
      bins(i) = this % map_undef % map(state)
    end do

    @assertEqual(EXPECTED_BINS,bins)

  end subroutine testMappingUndefined


  !!
  !! Test number of bins inquiry
  !!
@Test
  subroutine testNumberOfBinsInquiry(this)
    class(test_cellMap), intent(inout) :: this

    @assertEqual(3, this % map_noUndef % bins(1),'cellMap without undefined bin')
    @assertEqual(4, this % map_undef % bins(1),'cellMap with undefined bin')
    @assertEqual(3, this % map_noUndef % bins(0), 'Number of all bins')
    @assertEqual(0, this % map_noUndef % bins(2), 'higher dimension')
    @assertEqual(0, this % map_noUndef % bins(-2),'invalid dimension')

  end subroutine testNumberOfBinsInquiry

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequance is without errors
  !!
@Test
  subroutine testPrint(this)
    class(test_cellMap), intent(inout) :: this
    type(outputFile)                   :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_noUndef % print(out)
    @assertTrue(out % isValid(),'For map with no undefined cell bin: ')
    call out % reset()

    call this % map_undef % print(out)
    @assertTrue(out % isValid(),'For map with undefined cell bin: ')
    call out % reset()

  end subroutine testPrint



end module cellMap_test
