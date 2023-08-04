module cellMap_test
  use numPrecision
  use pFUnit_mod
  use universalVariables, only : VOID_MAT
  use particle_class,     only : particleState
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use charMap_class,      only : charMap
  use outputFile_class,   only : outputFile
  use cellMap_class,      only : cellMap
  use geometryReg_mod,    only : gr_kill => kill
  use geometryFactory_func, only : new_geometry
  use materialMenu_mod,     only : mm_nameMap => nameMap

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
  character(*), parameter :: SURF_DEF = "&
  & surf2 {id 8; type zPlane; z0 -1.3;} &
  & surf3 {id 9; type zPlane; z0 0.0;} &
  & surf4 {id 10; type zPlane; z0 1.0;}"

  character(*), parameter :: CELL_DEF = "&
  & cell2 {id 2; type simpleCell; surfaces (-8 10); filltype mat; material fuel;} &
  & cell3 {id 1; type simpleCell; surfaces (-8 -10 9); filltype mat; material void;} &
  & cell4 {id 5; type simpleCell; surfaces (-9 8); filltype mat; material fuel;} &
  & cell5 {id 3; type simpleCell; surfaces (-8); filltype mat; material pecorino;}"

  character(*), parameter :: UNI_DEF = "&
  & root {id 2; type rootUniverse; border 9; fill u<1>;} &
  & uni2 {id 1; type cellUniverse; cells (2 1 5 3);} "

  integer(shortInt),dimension(*),parameter :: CELL_IN_MAP = [2, 5, 3]

contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_cellMap), intent(inout) :: this
    type(dictionary)                   :: dict, mapDict1, dictTemp
    character(nameLen)                 :: name

    call dict % init(6)
    call dict % store('type','geometryStd')
    call dict % store('boundary', [0, 0, 0, 0, 0, 0])
    ! Store graph
    call dictTemp % init(1)
    call dictTemp % store('type','shrunk')
    call dict % store('graph', dictTemp)
    call dictTemp % kill()

    ! Store surfaces, cells and universes dictionaries
    call charToDict(dictTemp, SURF_DEF)
    call dict % store('surfaces', dictTemp)
    call charToDict(dictTemp, CELL_DEF)
    call dict % store('cells', dictTemp)
    call charToDict(dictTemp, UNI_DEF)
    call dict % store('universes', dictTemp)

    ! Initialise material map for materialMenu
    name = 'fuel'
    call mm_nameMap % add(name, 1)
    name = 'pecorino'
    call mm_nameMap % add(name, 2)
    name = 'void'
    call mm_nameMap % add(name, VOID_MAT)

    ! Initialise geometry in geomReg
    name = 'geom'
    call new_geometry(dict, name, silent = .true.)

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
    call gr_kill()

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
    integer(shortInt),dimension(5),parameter :: EXPECTED_BINS = [1, 0, 2, 3, 0]

    do i = 1, size(EXPECTED_BINS)
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
    integer(shortInt),dimension(5),parameter :: EXPECTED_BINS = [1, 4, 2, 3, 4]

    do i = 1, size(EXPECTED_BINS)
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
