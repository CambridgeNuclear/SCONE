module vectorField_inter

  use numPrecision
  use field_inter, only : field
  use coord_class, only : coordList

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: vectorField_CptrCast

  !!
  !! Simple Real Vector Field
  !!
  !! Access to field is via coordList to allow more fancy fields to be defined
  !! (e.g. assign vector to each uniqueID etc.)
  !!
  !! Interface:
  !!   field interface
  !!   at -> Return vector value given position coordinates
  !!
  type, public, abstract, extends(field) :: vectorField
  contains
    procedure(at), deferred :: at
  end type vectorField

  abstract interface

    !!
    !! Get value of the vector field at the co-ordinate point
    !!
    !! Args:
    !!   coords [in] -> Coordinates of the position in the geometry
    !!
    !! Result:
    !!   Size 3 vector of real values.
    !!
    function at(self, coords) result(val)
      import :: vectorField, coordList, defReal
      class(vectorField), intent(in) :: self
      type(coordList), intent(in)    :: coords
      real(defReal), dimension(3)    :: val
    end function at

  end interface

contains

  !!
  !! Cast field pointer to vectorField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of vectorField
  !!   Pointer to source if source is vectorField class
  !!
  pure function vectorField_CptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    class(vectorField), pointer       :: ptr

    select type (source)
    class is (vectorField)
        ptr => source

      class default
        ptr => null()
    end select

  end function vectorField_CptrCast

end module vectorField_inter
