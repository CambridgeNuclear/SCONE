module homogMatMap_test
  use numPrecision
  use funit
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use dictParser_func,         only : charToDict
  use outputFile_class,        only : outputFile

  ! May not be ideal but there is a dependance on Global materialMenu
  use materialMenu_mod,        only : mm_init => init, mm_kill => kill

  use homogMatMap_class,       only : homogMatMap

  implicit none


@testCase
  type, extends(TestCase) :: test_homogMatMap
    private
    type(homogMatMap),allocatable :: map_noUndef
    type(homogMatMap),allocatable :: map_Undef
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_homogMatMap


  !!
  !! Test parameters
  !!
  character(*), dimension(*), parameter :: BIN_NAMES   = ['bin1','bin2','bin3']
  character(*), dimension(*), parameter :: MAT_IN_BIN1 = ['mat3','mat4']
  character(*), dimension(*), parameter :: MAT_IN_BIN2 = ['mat1']
  character(*), dimension(*), parameter :: MAT_IN_BIN3 = ['mat5']

  !!
  !! Material Definitions
  !!
  character(*), parameter :: DICT_DEF = &
  " mat1 { temp 17; composition {} } &
  & mat2 { temp 17; composition {} } &
  & mat3 { temp 17; composition {} } &
  & mat4 { temp 17; composition {} } &
  & mat5 { temp 17; composition {} } "



contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_homogMatMap), intent(inout) :: this
    type(dictionary)                  :: dict
    type(dictionary)                  :: mapDict1

    ! Build nuclear data
    call charToDict(dict, DICT_DEF)
    call mm_init(dict)

    ! Initialise dictionaries
    call mapDict1 % init(5)

    ! Build material map definition
    call mapDict1 % store('bins', BIN_NAMES)
    call mapDict1 % store(BIN_NAMES(1), MAT_IN_BIN1)
    call mapDict1 % store(BIN_NAMES(2), MAT_IN_BIN2)
    call mapDict1 % store(BIN_NAMES(3), MAT_IN_BIN3)
    allocate(this % map_noUndef, source = homogMatMap(mapDict1))

    call mapDict1 % store('undefBin','true')
    allocate(this % map_undef, source = homogMatMap(mapDict1))

  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_homogMatMap), intent(inout) :: this

    call mm_kill()
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
    class(test_homogMatMap), intent(inout)   :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(6)           :: bins
    integer(shortInt),dimension(6),parameter :: EXPECTED_BINS = [2, 0, 1, 1, 3, 0]

    do i = 1,6
      state % matIdx = i
      bins(i) = this % map_noUndef % map(state)
    end do

    @assertEqual(EXPECTED_BINS,bins)

  end subroutine testMappingNoUndefined


  !!
  !! Mapping test with undefined bin
  !!
@Test
  subroutine testMappingUndefined(this)
    class(test_homogMatMap), intent(inout)   :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(6)           :: bins
    integer(shortInt),dimension(6),parameter :: EXPECTED_BINS = [2, 4, 1, 1, 3, 4]

    do i = 1,6
      state % matIdx = i
      bins(i) = this % map_undef % map(state)
    end do

    @assertEqual(EXPECTED_BINS,bins)

  end subroutine testMappingUndefined


  !!
  !! Test number of bins inquiry
  !!
@Test
  subroutine testNumberOfBinsInquiry(this)
    class(test_homogMatMap), intent(inout) :: this

    @assertEqual(3, this % map_noUndef % bins(1), 'homogMatMap without undefined bin')
    @assertEqual(4, this % map_undef % bins(1),   'homogMatMap with undefined bin')
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
    class(test_homogMatMap), intent(inout) :: this
    type(outputFile)                       :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_noUndef % print(out)
    @assertTrue(out % isValid(),'For map with no undefined material bin: ')
    call out % reset()

    call this % map_undef % print(out)
    @assertTrue(out % isValid(),'For map with undefined material bin: ')
    call out % reset()

  end subroutine testPrint



end module homogMatMap_test
