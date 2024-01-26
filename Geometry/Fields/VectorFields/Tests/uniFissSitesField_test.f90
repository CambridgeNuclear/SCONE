module uniFissSitesField_test
  use numPrecision
  use funit
  use particle_class,           only : particle, particleState
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use geometry_inter,           only : geometry
  use RNG_class,                only : RNG
  use uniFissSitesField_class,  only : uniFissSitesField

  implicit none


@testCase
  type, extends(TestCase) :: test_uniFissSitesField
    private
    type(uniFissSitesField) :: ufsField
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_uniFissSitesField

  !!
  !! Map definition
  !!
  character(*), parameter :: DICT_DEF = &
  " type spaceMap;  axis z;  grid unstruct; &
    &bins (0.0 20.0 40.0 60.0 80.0); "

contains

  !!
  !! Sets up test_weightWindows object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_uniFissSitesField), intent(inout) :: this
    type(dictionary)                             :: dict, dictMap
    class(geometry), pointer                     :: geom
    class(RNG), pointer                          :: rand
    integer(shortInt)                            :: type

    call charToDict(dictMap, DICT_DEF)

    ! Initialise dictionaries
    call dict % init(3)

    ! Build material map definition
    call dict % store('type', 'uniFissSitesField')
    call dict % store('uniformMap', 1)
    call dict % store('map', dictMap)

    call this % ufsField % init(dict)
    call this % ufsField % estimateVol(geom, rand, type)

  end subroutine setUp

  !!
  !! Kills test_weightWindows object
  !!
  subroutine tearDown(this)
    class(test_uniFissSitesField), intent(inout) :: this

    call this % ufsField % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test retrieving the ufs values
  !!
@Test
  subroutine testGetValue(this)
    class(test_uniFissSitesField), intent(inout) :: this
    type(particle)                               :: p
    type(particleState)                          :: state
    real(defReal), dimension(3)                  :: bins, EXPECTED_BINS

    ! Test case in the map
    p % coords % lvl(1) % r = [0.5, 7.0, 50.0]

    bins = this % ufsField % at(p)
    EXPECTED_BINS = [0.25, 0.25, 0.0]
    @assertEqual(EXPECTED_BINS, bins, tolerance=1.0e-6)

    ! Test case outside the map
    p % coords % lvl(1) % r = [0.5, 7.0, 100.0]

    bins = this % ufsField % at(p)
    EXPECTED_BINS = [1.0, 1.0, 1.0]
    @assertEqual(EXPECTED_BINS,bins)

    ! Modify the map by storing fission sites
    state % r   = [0.5, 7.0, 12.0]
    state % wgt = 0.2

    call this % ufsField % storeFS(state)

    state % r   = [0.5, 7.0, 23.2]
    state % wgt = 0.8
    call this % ufsField % storeFS(state)

    call this % ufsField % updateMap()

    ! Test case in the updated map
    p % coords % lvl(1) % r = [0.5, 7.0, 18.1]

    bins = this % ufsField % at(p)
    EXPECTED_BINS = [0.25, 0.06666666667, 0.0]
    @assertEqual(EXPECTED_BINS, bins, tolerance=1.0e-6)

  end subroutine testGetValue


end module uniFissSitesField_test
