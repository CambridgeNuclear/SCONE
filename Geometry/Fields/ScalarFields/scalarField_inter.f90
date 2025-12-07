module scalarField_inter

  use numPrecision
  use field_inter,    only : field
  use coord_class,    only : coordList
  use particle_class, only : particle

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: scalarField_CptrCast

  !!
  !! Simple Real Scalar Field
  !!
  !! Access to field is via coordList to allow more fancy fields to be defined
  !! (e.g. assign value to each uniqueID etc.)
  !!
  !! Also allows access by particle until WW and UFS are sorted
  !!
  !! Interface:
  !!   field interface
  !!   at  -> Return scalar value given position coordinates
  !!   atP -> Return scalar value given particle
  !!
  type, public, abstract, extends(field) :: scalarField
  contains
    procedure(at), deferred  :: at
    procedure(atP), deferred :: atP
  end type scalarField

  abstract interface

    !!
    !! Get value of the scalar field at the co-ordinate point
    !!
    !! Args:
    !!   coords [in] -> Coordinates of the position in the geometry
    !!
    !! Result:
    !!   Value of the scalar field. Real number.
    !!
    function at(self, coords) result(val)
      import :: scalarField, coordList, defReal
      class(scalarField), intent(in)  :: self
      class(coordList), intent(in)    :: coords
      real(defReal)                   :: val
    end function at
    
    !!
    !! Get value of the scalar field at the particle's phase space point
    !!
    !! Args:
    !!   particle [in] -> Coordinates + other phase space info
    !!
    !! Result:
    !!   Value of the scalar field. Real number.
    !!
    function atP(self, p) result(val)
      import :: scalarField, particle, defReal
      class(scalarField), intent(in)  :: self
      class(particle), intent(in)     :: p
      real(defReal)                   :: val
    end function atP

  end interface

contains

  !!
  !! Cast field pointer to scalarField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of scalarField
  !!   Pointer to source if source is scalarField class
  !!
  pure function scalarField_CptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    class(scalarField), pointer       :: ptr

    select type (source)
    class is (scalarField)
        ptr => source

      class default
        ptr => null()
    end select

  end function scalarField_CptrCast

end module scalarField_inter
