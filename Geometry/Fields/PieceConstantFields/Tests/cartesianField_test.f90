module cartesianField_test

  use numPrecision
  use dictionary_class,         only : dictionary
  use dictParser_func,          only : charToDict
  use particle_class,           only : particle
  use field_inter,              only : field
  use pieceConstantField_inter, only : pieceConstantField, pieceConstantField_CptrCast
  use cartesianField_class,     only : cartesianField, cartesianField_TptrCast
  use funit

  implicit none
  character(*), parameter :: FIELD_DEF = "type cartesianField; origin (1.0 2.0 1.0); pitch (1.0 2.0 0.5); &
          & materials (all); shape (3 2 3); all (1 2 3 1 2 3 4 5 6 4 5 6 7 8 9 7 8 9 ); default -6.2;"

contains

  !!
  !! Test Cartesian Field
  !!
@Test
  subroutine test_cartesianField()
    type(cartesianField), target       :: fieldT
    class(field), pointer              :: ref
    class(pieceConstantField), pointer :: ptr
    type(cartesianField), pointer      :: ptr2
    type(dictionary)                   :: dict
    type(particle)                     :: p
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Test invalid pointers
    ref => null()

    ptr => pieceConstantField_CptrCast(ref)
    ptr2 => cartesianField_TptrCast(ref)
    @assertFalse(associated(ptr))
    @assertFalse(associated(ptr2))

    ! Test valid pointers
    ref => fieldT

    ptr => pieceConstantField_CptrCast(ref)
    ptr2 => cartesianField_TptrCast(ref)

    @assertTrue(associated(ptr, fieldT))
    @assertTrue(associated(ptr2, fieldT))

    ! Initialise field
    call charToDict(dict, FIELD_DEF)

    call fieldT % init(dict)

    ! Check values at different points
    ! Inside region x = 2, y = 1, z = 1
    call p % point([1.0_defReal, 0.0_defReal, 0.0_defReal])
    call p % teleport([1.0_defReal, 1.0_defReal, 0.4_defReal])
    @assertEqual(8.0_defReal, fieldT % at(p), TOL)
    
    ! Inside region x = 1, y = 2, z = 3
    call p % point([1.0_defReal, 0.0_defReal, 0.0_defReal])
    call p % teleport([-0.4_defReal, 2.2_defReal, 1.3_defReal])
    @assertEqual(1.0_defReal, fieldT % at(p), TOL)
    
    ! Outside the field
    call p % teleport([-8.5_defReal, 4.5_defReal, 0.9_defReal])
    @assertEqual(-6.2_defReal, fieldT % at(p), TOL)

    ! Check distances to the field
    ! Inside, pointing along x
    call p % point([1.0_defReal, 0.0_defReal, 0.0_defReal])
    call p % teleport([1.0_defReal, 3.5_defReal, 0.9_defReal])
    @assertEqual(0.5_defReal, fieldT % distance(p), TOL)

    ! Inside, at an angle
    call p % point([-sqrt(2.0_defReal)/2, sqrt(2.0_defReal)/2, 0.0_defReal])
    call p % teleport([0.9_defReal, 1.0_defReal, 0.4_defReal])
    @assertEqual(0.4_defReal * sqrt(2.0_defReal), fieldT % distance(p), TOL)

    ! Outside, at an angle
    call p % point([-sqrt(2.0_defReal)/2, -sqrt(2.0_defReal)/2, 0.0_defReal])
    call p % teleport([2.0_defReal, 4.5_defReal, 1.0_defReal])
    @assertEqual(1.0_defReal/sqrt(2.0_defReal), fieldT % distance(p), TOL)

    ! Kill
    call fieldT % kill()

  end subroutine test_cartesianField

end module cartesianField_test
