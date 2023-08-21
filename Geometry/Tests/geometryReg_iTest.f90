module geometryReg_iTest

  use numPrecision
  use genericProcedures, only : numToChar
  use dictionary_class,  only : dictionary
  use dictParser_func,   only : charToDict
  use geometry_inter,    only : geometry
  use field_inter,       only : field
  use materialMenu_mod,  only : mm_init => init, mm_kill => kill
  use geometryReg_mod,   only : gr_geomIdx => geomIdx, gr_geomPtr => geomPtr,&
                                gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr, gr_kill => kill
  use geometryFactory_func, only : new_geometry
  use fieldFactory_func,    only : new_field
  use funit
  implicit none

  ! Material definition
  character(*), parameter :: MAT_DEF = &
  " fuel  {temp 8; composition {}} &
  & water {temp 2; composition{}}"

  ! Geometry definition
  character(*), parameter :: GEOM_DEF = &
  " type geometryStd; boundary (0); graph {type shrunk;} &
  & surfaces { surf1 {id 1; type sphere; origin (0.0 0.0 0.0); radius 2.0; }} &
  & cells {} &
  & universes { root {id 1; type rootUniverse; border 1; fill fuel;}}"

  ! Field definition
  character(*), parameter :: FIELD_DEF = &
  "type uniformScalarField; value -1.333;"




contains

  !!
  !! Prepare test enviroment
  !!
@Before
  subroutine set_up()
    type(dictionary) :: dict

    ! Load material definitions
    call charToDict(dict, MAT_DEF)
    call mm_init(dict)

  end subroutine set_up

  !!
  !! Clean test enviroment
  !!
@After
  subroutine clean_up()

    call mm_kill()

  end subroutine clean_up

  !!
  !! Test geometry functionality
  !!
@Test
  subroutine test_geometry()
    type(dictionary)         :: dict
    character(nameLen)       :: name
    integer(shortInt)        :: i
    class(geometry), pointer :: ptr

    ! Get geometry definition
    call charToDict(dict, GEOM_DEF)

    ! Build 10 instances of geometry
    do i = 1, 10
      name = 'geom'//numToChar(i)
      call new_geometry(dict, name, silent=.true.)
    end do

    ! Get some indexes
    name = 'geom7'
    @assertEqual(7, gr_geomIdx(name))

    name = 'geom2'
    @assertEqual(2, gr_geomIdx(name))

    ! Get some pointers -> Check if associated
    ptr => gr_geomPtr(3)
    @assertTrue(associated(ptr))

    ptr => gr_geomPtr(10)
    @assertTrue(associated(ptr))

    ! Kill registry
    call gr_kill()

  end subroutine test_geometry

  !!
  !! Test field functionality
  !!
@Test
  subroutine test_field()
    type(dictionary)      :: dict
    character(nameLen)    :: name
    integer(shortInt)     :: i
    class(field), pointer :: ptr

    ! Get field definition
    call charToDict(dict, FIELD_DEF)

    ! Build 10 instances of fields
    do i = 1, 10
      name = 'field'//numToChar(i)
      call new_field(dict, name)
    end do

    ! Get some indexes
    name = 'field7'
    @assertEqual(7, gr_fieldIdx(name))

    name = 'field2'
    @assertEqual(2, gr_fieldIdx(name))

    ! Get some pointers -> Check if associated
    ptr => gr_fieldPtr(3)
    @assertTrue(associated(ptr))

    ptr => gr_fieldPtr(10)
    @assertTrue(associated(ptr))

    ! Kill registry
    call gr_kill()

  end subroutine test_field


end module geometryReg_iTest
