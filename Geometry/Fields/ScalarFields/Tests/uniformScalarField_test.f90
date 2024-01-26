module uniformScalarField_test

  use numPrecision
  use dictionary_class,         only : dictionary
  use particle_class,           only : particle
  use field_inter,              only : field
  use scalarField_inter,        only : scalarField, scalarField_CptrCast
  use uniformScalarField_class, only : uniformScalarField, uniformScalarField_TptrCast
  use funit

  implicit none

contains

  !!
  !! Test Uniform Scalar Field
  !!
@Test
  subroutine test_uniformScalarField()
    type(uniformScalarField), target  :: fieldT
    class(field), pointer             :: ref
    class(scalarField), pointer       :: ptr
    type(uniformScalarField), pointer :: ptr2
    type(dictionary)                  :: dict
    type(particle)                    :: p
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Test invalid pointers
    ref => null()

    ptr => scalarField_CptrCast(ref)
    ptr2 => uniformScalarField_TptrCast(ref)
    @assertFalse(associated(ptr))
    @assertFalse(associated(ptr2))

    ! Test valid pointers
    ref => fieldT

    ptr => scalarField_CptrCast(ref)
    ptr2 => uniformScalarField_TptrCast(ref)

    @assertTrue(associated(ptr, fieldT))
    @assertTrue(associated(ptr2, fieldT))

    ! Initialise field
    call dict % init(2)
    call dict % store('type', 'uniformVectorField')
    call dict % store('value', 9.6_defReal)

    call fieldT % init(dict)

    ! Check value
    @assertEqual(9.6_defReal, fieldT % at(p), TOL)

    ! Kill
    call fieldT % kill()

  end subroutine test_uniformScalarField

end module uniformScalarField_test
