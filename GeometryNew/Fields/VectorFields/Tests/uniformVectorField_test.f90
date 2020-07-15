module uniformVectorField_test

  use numPrecision
  use dictionary_class,         only : dictionary
  use coord_class,              only : coordList
  use field_inter,              only : field
  use vectorField_inter,        only : vectorField, vectorField_CptrCast
  use uniformVectorField_class, only : uniformVectorField, uniformVectorField_TptrCast
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test Uniform Scalar Field
  !!
@Test
  subroutine test_uniformVectorField()
    type(uniformVectorField), target  :: fieldT
    class(field), pointer             :: ref
    class(vectorField), pointer       :: ptr
    type(uniformVectorField), pointer :: ptr2
    type(dictionary)                  :: dict
    type(coordList)                   :: coords
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Test invalid pointers
    ref => null()

    ptr => vectorField_CptrCast(ref)
    ptr2 => uniformVectorField_TptrCast(ref)
    @assertFalse(associated(ptr))
    @assertFalse(associated(ptr2))

    ! Test valid pointers
    ref => fieldT

    ptr => vectorField_CptrCast(ref)
    ptr2 => uniformVectorField_TptrCast(ref)

    @assertTrue(associated(ptr, fieldT))
    @assertTrue(associated(ptr2, fieldT))

    ! Initialise field
    call dict % init(2)
    call dict % store('type', 'uniformVectorField')
    call dict % store('value', [9.6_defReal, -8.0_defReal, 9.7_defReal])

    call fieldT % init(dict)

    ! Check value
    call coords % init([ONE, ZERO, ONE], [ONE, ZERO, ZERO])
    @assertEqual([9.6_defReal, -8.0_defReal, 9.7_defReal], fieldT % at(coords), TOL)

    ! Kill
    call fieldT % kill()

  end subroutine test_uniformVectorField

end module uniformVectorField_test
