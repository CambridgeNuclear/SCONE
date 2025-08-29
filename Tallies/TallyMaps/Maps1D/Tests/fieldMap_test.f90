module fieldMap_test
  use numPrecision
  use funit
  use particle_class,     only : particleState
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use outputFile_class,   only : outputFile
  use fieldMap_class,     only : fieldMap

  implicit none


@testCase
  type, extends(TestCase) :: test_fieldMap
    private
    type(fieldMap),allocatable :: map_cartesian
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_fieldMap

  !!
  !! Test parameters
  !!
  character(*), parameter :: FIELD_DEF = "&
  & type cartesianField; origin (0 0 0); &
  & pitch (3 2 1); materials (all); default 1; &
  & shape (2 2 2); all (1 1 1 1 1 1 1 1);"

contains

  !!
  !! Sets up test_fieldMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_fieldMap), intent(inout) :: this
    type(dictionary)                    :: dict, mapDict1, dictTemp
    character(nameLen)                  :: name

    ! Store surfaces, cells and universes dictionaries
    call charToDict(dictTemp, FIELD_DEF)

    ! Initialise dictionaries
    call mapDict1 % init(2)

    ! Build material map definition
    call mapDict1 % store('field', dictTemp)
    allocate(this % map_cartesian, source = fieldMap(mapDict1))

  end subroutine setUp

  !!
  !! Kills test_fieldMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_fieldMap), intent(inout) :: this

    call this % map_cartesian % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Mapping test without undefined bin
  !!
@Test
  subroutine testMapping(this)
    class(test_fieldMap), intent(inout)      :: this
    type(particleState)                      :: state
    integer(shortInt)                        :: i
    integer(shortInt),dimension(4)           :: bins
    integer(shortInt),dimension(4),parameter :: EXPECTED_BINS = [7, 1, 6, 0]
    real(defReal), dimension(12), parameter :: r = [-2.0, 1.0, 0.5, -1.0, -1.0, -0.5, 1.4, -0.1, 0.1, 10.0, 10.0, 10.0]

    ! Note ordering of cartesian field
    do i = 1, size(EXPECTED_BINS)
      state % r = r((3*(i-1) + 1):(3*i))
      bins(i) = this % map_cartesian % map(state)
    end do

    @assertEqual(EXPECTED_BINS,bins)

  end subroutine testMapping


  !!
  !! Test number of bins inquiry
  !!
@Test
  subroutine testNumberOfBinsInquiry(this)
    class(test_fieldMap), intent(inout) :: this

    @assertEqual(8, this % map_cartesian % bins(1))
    @assertEqual(0, this % map_cartesian % bins(2), 'higher dimension')
    @assertEqual(0, this % map_cartesian % bins(-2),'invalid dimension')

  end subroutine testNumberOfBinsInquiry

  !!
  !! Test correctness of print subroutine
  !! Does not checks that values are correct, but that calls sequance is without errors
  !!
@Test
  subroutine testPrint(this)
    class(test_fieldMap), intent(inout) :: this
    type(outputFile)                   :: out

    call out % init('dummyPrinter', fatalErrors = .false.)

    call this % map_cartesian % print(out)
    @assertTrue(out % isValid())
    call out % reset()

  end subroutine testPrint

end module fieldMap_test
