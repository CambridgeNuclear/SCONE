module uniformScalarField_class

  use numPrecision
  use dictionary_class,  only : dictionary
  use coord_class,       only : coordList
  use particle_class,    only : particle
  use field_inter,       only : field
  use scalarField_inter, only : scalarField

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: uniformScalarField_TptrCast

  !!
  !! Uniform Scalar Field
  !!
  !! Always returns the same value
  !!
  !! Sample Dictionary Input:
  !!   field { type uniformScalarField; value 3.0;}
  !!
  !! Public Members:
  !!   val -> Value of the field
  !!
  !! Interface:
  !!   scalarField interface
  !!
  type, public, extends(scalarField) :: uniformScalarField
    real(defReal) :: val = ZERO
  contains
    ! Superclass interface
    procedure :: init
    procedure :: kill
    procedure :: at
    procedure :: atP
  end type uniformScalarField

contains

  !!
  !! Initialise from dictionary
  !!
  !! See field_inter for details
  !!
  subroutine init(self, dict)
    class(uniformScalarField), intent(inout) :: self
    class(dictionary), intent(in)            :: dict

    ! Load value
    call dict % get(self % val, 'value')

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(uniformScalarField), intent(inout) :: self

    self % val = ZERO

  end subroutine kill

  !!
  !! Get value of the scalar field at the co-ordinate point
  !!
  !! See scalarField_inter for details
  !!
  function at(self, coords) result(val)
    class(uniformScalarField), intent(in) :: self
    class(coordList), intent(in)          :: coords
    real(defReal)                         :: val

    val = self % val

  end function at
  
  !!
  !! Get value of the scalar field at particle phase space location
  !!
  !! See scalarField_inter for details
  !!
  function atP(self, p) result(val)
    class(uniformScalarField), intent(in) :: self
    class(particle), intent(in)           :: p
    real(defReal)                         :: val

    val = self % val

  end function atP

  !!
  !! Cast field pointer to uniformScalarField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of uniformScalarField
  !!   Pointer to source if source is uniformScalarField type
  !!
  pure function uniformScalarField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(uniformScalarField), pointer :: ptr

    select type (source)
      type is (uniformScalarField)
        ptr => source

      class default
        ptr => null()
    end select

  end function uniformScalarField_TptrCast


end module uniformScalarField_class
