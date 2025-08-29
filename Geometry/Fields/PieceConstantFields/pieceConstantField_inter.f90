module pieceConstantField_inter

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use field_inter,       only : field
  use coord_class,       only : coordList
  use dictionary_class,  only : dictionary
  use particle_class,    only : particle
  use dictParser_func,   only : fileToDict

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: pieceConstantField_CptrCast

  !!
  !! Piecewise constant field. Values of the field are piecewise constant.
  !! Provides a distance calculation to neighbouring elements, allowing surface
  !! tracking to be used.
  !!
  !! Access to field is via coordList to allow more fancy fields to be defined
  !! (e.g. assign value to each uniqueID etc.)
  !!
  !! Interface:
  !!   field interface
  !!   init         -> Initialise from dictionary, which may contain a path
  !!   init_dict    -> Initialises from a dictionary containing full description
  !!   init_file    -> Initialises from a path to a full description
  !!   at           -> Return scalar value given coordlist
  !!   map          -> Returns field index given coordList  
  !!   map_particle -> Returns field index given particle
  !!   distance     -> Return distance to next element of the field
  !!   setValue     -> Sets value of field at a given index
  !!   getMaxValue  -> Returns the maximum field value
  !!   getSize      -> Returns the size of the field
  !!
  type, public, abstract, extends(field) :: pieceConstantField
    real(defReal), dimension(:), allocatable :: val
    integer(shortInt)                        :: N = 0
  contains
    procedure(init_dict), deferred :: init_dict
    procedure(at), deferred        :: at
    procedure(distance), deferred  :: distance
    procedure(map), deferred       :: map
    
    procedure :: init
    procedure :: init_file
    procedure :: map_particle
    procedure :: setValue
    procedure :: getMaxValue
    procedure :: getSize
    
    procedure :: killSuper
  end type pieceConstantField

  abstract interface

    !!
    !! Initialise field from a dictionary
    !!
    !! Args:
    !!   dict [in] -> Dictionary describing field
    !!
    subroutine init_dict(self, dict)
      import :: pieceConstantField, dictionary
      class(pieceConstantField), intent(inout) :: self
      class(dictionary), intent(in)            :: dict
    end subroutine init_dict
    
    !!
    !! Get value of the field at the given coordinates
    !!
    !! Args:
    !!   coords [in] -> Coordinates of the position in the geometry
    !!
    !! Result:
    !!   Value of the scalar field. Real number.
    !!
    function at(self, coords) result(val)
      import :: pieceConstantField, coordList, defReal
      class(pieceConstantField), intent(in) :: self
      class(coordList), intent(in)          :: coords
      real(defReal)                         :: val
    end function at
    
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
    
    !!
    !! Get index of the field at the given coordinates
    !!
    !! Args:
    !!   coords [in] -> Coordinates of the position in the geometry
    !!
    !! Result:
    !!   Index into the field. Integer
    !!
    pure function map(self, coords) result(idx)
      import :: pieceConstantField, coordList, shortInt
      class(pieceConstantField), intent(in) :: self
      class(coordList), intent(in)          :: coords
      integer(shortInt)                     :: idx
    end function map
    
  end interface

contains

  !!
  !! Initialise field from a dictionary. This dictionary may contain
  !! either the details of the field or a path to a further dictionary
  !!
  !! Args:
  !!   dict [in] -> Dictionary describing contents of field
  !!
  subroutine init(self, dict)
    class(pieceConstantField), intent(inout) :: self
    class(dictionary), intent(in)            :: dict
    character(pathLen)                       :: path

    if (dict % isPresent('path')) then
      call dict % get(path, 'path')
      call self % init_file(path)
    else
      call self % init_dict(dict)
    end if

  end subroutine init
  
  !!
  !! Initialise field from a file
  !!
  !! Args:
  !!   file [in] -> Path to a file containing the field definition
  !!
  subroutine init_file(self, path)
    class(pieceConstantField), intent(inout) :: self
    character(pathLen), intent(in)           :: path
    type(dictionary)                         :: dict

    call fileToDict(dict, path)
    call self % init(dict)

  end subroutine init_file

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
  !! Returns the index of the field given a particle.
  !!
  function map_particle(self, p) result(idx)
    class(pieceConstantField), intent(in) :: self
    type(particle), intent(in)            :: p
    integer(shortInt)                     :: idx
    type(coordList)                       :: coords

    coords = p % coords
    idx = self % map(coords)

  end function map_particle

  !!
  !! Returns the value of the maximum value of the field.
  !! This may be useful, e.g., in ensuring the majorant is correct.
  !!
  !! Errors:
  !!  If the values of the field are not allocated
  !!
  function getMaxValue(self) result(val)
    class(pieceConstantField), intent(in) :: self
    real(defReal)                         :: val
    character(100), parameter             :: Here = 'getMaxValue (pieceConstantField_class.f90)'

    if (.not. allocated(self % val)) call fatalError(Here,&
            'The values of the field have not been allocated')

    val = maxval(self % val)

  end function getMaxValue  
  
  !!
  !! Returns the size of the field.
  !!
  function getSize(self) result(N)
    class(pieceConstantField), intent(in) :: self
    integer(shortInt)                     :: N

    N = self % N

  end function getSize
  
  !!
  !! Kills elements of the superclass
  !!
  elemental subroutine killSuper(self)
    class(pieceConstantField), intent(inout) :: self

    if (allocated(self % val)) deallocate(self % val)
    self % N = 0

  end subroutine killSuper

end module pieceConstantField_inter
