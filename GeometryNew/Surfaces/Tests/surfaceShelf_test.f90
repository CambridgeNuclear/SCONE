module surfaceShelf_test

  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict
  use surface_inter,      only : surface
  use surfaceShelf_class, only : surfaceShelf
  use pfUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: SHELF_DEF = &
    " surf1 {id 18; type sphere; origin (1.0 1.0 1.0); radius 17.4;} &
    & surf2 {id 23; type xPlane; x0 3.1;} &
    & surf3 {id 1; type zPlane; z0  0.3;} "

  ! Variables
  type(surfaceShelf) :: shelf


contains

  !!
  !! Build the Shelf
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict

    ! Try killing uninitialised
    call shelf % kill()

    call charToDict(dict, SHELF_DEF)
    call shelf % init(dict)

  end subroutine setUp

  !!
  !! Clean after SHelf
  !!
@After
  subroutine cleanUp()

    call shelf % kill()

  end subroutine cleanUp

  !!
  !! Test
  !!
@Test
  subroutine testShelf()
    class(surface), pointer :: ptr
    integer(shortInt)       :: idx

    ! 1st Surface
    idx = shelf % surfIdx(18)
    ptr => shelf % surfPtr(idx)
    @assertEqual(18, ptr % id())
    @assertEqual('sphere', ptr % myType())
    @assertEqual(18, shelf % surfId(idx))

    ! 2nd Surface
    idx = shelf % surfIdx(1)
    ptr => shelf % surfPtr(idx)
    @assertEqual(1, ptr % id())
    @assertEqual('zPlane', ptr % myType())
    @assertEqual(1, shelf % surfId(idx))

  end subroutine testShelf


end module surfaceShelf_test
