module vectorField_inter

  use numPrecision
  use field_inter,    only : field
  use coord_class,    only : coordList
  use particle_class, only : particle

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
  !! Also allows access by particle until WW and UFS are sorted
  !!
  !! Interface:
  !!   field interface
  !!   at  -> Return vector value given position coordinates
  !!   atP -> Return vector value given particle
  !!
  type, public, abstract, extends(field) :: vectorField
  contains
    procedure(at), deferred  :: at
    procedure(atP), deferred :: atP
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
      class(coordList), intent(in)   :: coords
      real(defReal), dimension(3)    :: val
    end function at

    !!
    !! Get value of the vector field by particle
    !!
    !! Args:
    !!   particle [in] -> Coordinates + other phase space info
    !!
    !! Result:
    !!   Size 3 vector of real values.
    !!
    function atP(self, p) result(val)
      import :: vectorField, particle, defReal
      class(vectorField), intent(in) :: self
      class(particle), intent(in)    :: p
      real(defReal), dimension(3)    :: val
    end function atP

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
