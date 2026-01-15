module pieceConstantField_inter

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use scalarField_inter, only : scalarField
  use coord_class,       only : coordList
  use dictionary_class,  only : dictionary
  use field_inter,       only : field
  use scalarField_inter, only : scalarField

  implicit none
  private

  !!
  !! Public Pointer Cast and clean up
  !!
  public :: pieceConstantField_CptrCast
  public :: kill

  !!
  !! Piecewise constant field. Values of the field are piecewise constant.
  !! Provides a distance calculation to neighbouring elements, allowing surface
  !! tracking to be used.
  !!
  !! Interface:
  !!   field interface
  !!   distance    -> Return distance to next element of the field
  !!   setValue    -> Sets value of field at a given index
  !!   getMaxValue -> Returns the maximum field value
  !!
  type, public, abstract, extends(scalarField) :: pieceConstantField
    real(defReal), dimension(:), allocatable   :: val
    integer(shortInt)                          :: N = 0
  contains
    procedure(distance), deferred :: distance
    
    procedure :: setValue
    procedure :: getMaxValue
    
    procedure :: kill
  end type pieceConstantField

  abstract interface

    !!
    !! Get distance to the next element of the field at the given coordinates
    !!
    !! Args:
    !!   coords [in] -> Coordinates of the position in the geometry
    !!
    !! Result:
    !!   Distance to the next element of the field. Real number.
    !!
    function distance(self, coords) result(d)
      import :: pieceConstantField, coordList, defReal
      class(pieceConstantField), intent(in) :: self
      class(coordList), intent(in)          :: coords
      real(defReal)                         :: d
    end function distance
    
  end interface

contains

  !!
  !! Cast field pointer to pieceConstantField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of pieceConstantField
  !!   Pointer to source if source is pieceConstantField class
  !!
  pure function pieceConstantField_CptrCast(source) result(ptr)
    class(field), pointer, intent(in)  :: source
    class(pieceConstantField), pointer :: ptr

    select type (source)
    class is (pieceConstantField)
        ptr => source

      class default
        ptr => null()
    end select

  end function pieceConstantField_CptrCast
  
  !!
  !! Sets the value of a field at a given index.
  !! May be used to update the field.
  !!
  !! Args:
  !!   val [in] -> Value to input
  !!   idx [in] -> Index of the field to overwrite
  !!
  !! Errors:
  !!  If the values of the field are not allocated
  !!  If the index is out of bounds
  !!
  subroutine setValue(self, val, idx)
    class(pieceConstantField), intent(inout) :: self
    real(defReal), intent(in)                :: val
    integer(shortInt), intent(in)            :: idx
    character(100), parameter                :: Here = 'setValue_idx (pieceConstantField_class.f90)'

    if (.not. allocated(self % val)) call fatalError(Here,&
            'The values of the field have not been allocated')

    if (idx < 1 .or. idx > self % N) call fatalError(Here,&
            'Invalid index to overwrite: '//numToChar(idx)&
            //'. Field has size: '//numToChar(self % N))

    self % val(idx) = val

  end subroutine setValue

  !!
  !! Returns the value of the maximum value of the field.
  !! This may be useful, e.g., in ensuring the majorant is correct.
  !!
  !! Errors:
  !!  If the values of the field are not allocated
  !!
  function getMaxValue(self) result(val)
    class(pieceConstantField), intent(inout) :: self
    real(defReal)                            :: val
    character(100), parameter                :: Here = 'getMaxValue (pieceConstantField_class.f90)'

    if (.not. allocated(self % val)) call fatalError(Here,&
            'The values of the field have not been allocated')

    val = maxval(self % val)

  end function getMaxValue  
  
  !!
  !! Kills elements of the superclass
  !!
  elemental subroutine kill(self)
    class(pieceConstantField), intent(inout) :: self

    if (allocated(self % val)) deallocate(self % val)
    self % N = 0

  end subroutine kill

end module pieceConstantField_inter
