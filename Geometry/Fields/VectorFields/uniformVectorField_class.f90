module uniformVectorField_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use dictionary_class,  only : dictionary
  use coord_class,       only : coordList
  use particle_class,    only : particle
  use field_inter,       only : field
  use vectorField_inter, only : vectorField

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: uniformVectorField_TptrCast

  !!
  !! Uniform Vector Field
  !!
  !! Always returns the same value
  !!
  !! Sample Dictionary Input:
  !!   field { type uniformVectorField; value (3.0 1.0 0.3);}
  !!
  !! Public Members:
  !!   val -> Value of the field
  !!
  !! Interface:
  !!   vectorField interface
  !!
  type, public, extends(vectorField) :: uniformVectorField
    real(defReal), dimension(3) :: val = ZERO
  contains
    ! Superclass interface
    procedure :: init
    procedure :: kill
    procedure :: at
    procedure :: atP
  end type uniformVectorField

contains

  !!
  !! Initialise from dictionary
  !!
  !! See field_inter for details
  !!
  subroutine init(self, dict)
    class(uniformVectorField), intent(inout) :: self
    class(dictionary), intent(in)            :: dict
    real(defReal), dimension(:), allocatable :: temp
    character(100), parameter :: Here = 'init (uniformVectorField_class.f90)'

    ! Load value
    call dict % get(temp, 'value')

    if (size(temp) /= 3) then
      call fatalError(Here, 'Value must have size 3. Has: '//numToChar(size(temp)))
    end if

    self % val = temp

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(uniformVectorField), intent(inout) :: self

    self % val = ZERO

  end subroutine kill

  !!
  !! Get value of the scalar field at the co-ordinate point
  !!
  !! See vectorField_inter for details
  !!
  function at(self, coords) result(val)
    class(uniformVectorField), intent(in) :: self
    class(coordList), intent(in)          :: coords
    real(defReal), dimension(3)           :: val

    val = self % val

  end function at
  
  !!
  !! Get value of the scalar field by particle
  !!
  !! See vectorField_inter for details
  !!
  function atP(self, p) result(val)
    class(uniformVectorField), intent(in) :: self
    class(particle), intent(in)           :: p
    real(defReal), dimension(3)           :: val

    val = self % val

  end function atP

  !!
  !! Cast field pointer to uniformVectorField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of uniformVectorField
  !!   Pointer to source if source is uniformVectorField type
  !!
  pure function uniformVectorField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(uniformVectorField), pointer :: ptr

    select type (source)
      type is (uniformVectorField)
        ptr => source

      class default
        ptr => null()
    end select

  end function uniformVectorField_TptrCast


end module uniformVectorField_class
