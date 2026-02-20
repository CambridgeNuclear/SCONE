!!
!! Registry that contains all defined geometries and fields
!!
!! Stores and manages liftime of all geometries and fields
!! Allows to obtain pointer to a geometry or a field
!! Serves as geometry and field factory
!!
!! It is necessary to load geometries and fields one-by-one.
!!
!! Private Members:
!!   geometries      -> Array with the available geometries
!!   geometryTop     -> Number of defined geometries
!!   geometryNameMap -> Map of geometry names to idx
!!   fields          -> Array with the available fields
!!   fieldTop        -> Number of defined fields
!!   fieldNameMap    -> Map of field names to idx
!!
!! Interface:
!!   addGeom  -> Add new geometry
!!   geomIdx  -> Get index of the geometry
!!   geomPtr  -> Get pointer to a geometry with an index
!!   addField -> Add new field
!!   hasField -> Does the field with a given name exist?
!!   fieldIdx -> Get index of a field
!!   fieldPtr -> Get pointer to a field specified by index
!!   display  -> Display info about defined fields and geometries
!!   kill     -> Return to uninitialised state
!!
module geometryReg_mod

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use dictionary_class,  only : dictionary
  use charMap_class,     only : charMap

  ! Geometry
  use geometry_inter,    only : geometry

  ! Fields
  use field_inter,       only : field

  implicit none
  private

  !!
  !! Small type to hold a single polymorphic geometry
  !!
  type :: geomBox
    character(nameLen)           :: name = ''
    class(geometry), allocatable :: geom
  end type geomBox

  !!
  !! Small type to hold a single polymorphic field
  !!   Kentta ~= Field (Fi.)
  !!
  type :: fieldBox
    character(nameLen)        :: name = ''
    class(field), allocatable :: kentta
  end type fieldBox

  !! Public Interface
  public :: addGeom
  public :: geomIdx
  public :: geomPtr
  public :: geomNum
  public :: addField
  public :: fieldIdx
  public :: fieldPtr
  public :: hasField
  public :: kill

  integer(shortInt), parameter :: START_SIZE = 5
  real(defReal), parameter     :: GROWTH_RATE = 1.6_defReal

  !! Members
  type(geomBox), dimension(:), allocatable, target  :: geometries
  integer(shortInt)                                 :: geometryTop = 0
  type(charMap)                                     :: geometryNameMap

  type(fieldBox), dimension(:), allocatable, target :: fields
  integer(shortInt)                                 :: fieldTop = 0
  type(charMap)                                     :: fieldNameMap

contains

  !!
  !! Add Geometry definition
  !!
  !! Gets information about materials from materialMenu
  !!
  !! Args:
  !!   geom [in] -> Pre-initialised geometry, to be allocated
  !!   name [in] -> Name of the geometry
  !!
  !! Errors:
  !!   fatalError if geometry with the name was already defined
  !!
  subroutine addGeom(geom, name)
    class(geometry), allocatable, intent(inout) :: geom
    character(nameLen), intent(in)              :: name
    integer(shortInt)                           :: idx
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter    :: Here = 'addGeom (geometryReg_mod.f90)'

    ! Get free index
    idx = geometryNameMap % getOrDefault(name, NOT_PRESENT)

    if (idx /= NOT_PRESENT) then
      call fatalError(Here, 'Geometry with name: '//trim(name)//' has already been defined with &
                            &idx: '//numToChar(idx))
    else
      idx = freeGeomIdx()
      geometryTop = geometryTop + 1 ! Increment counter
    end if

    ! Initialise & store index
    call geometryNameMap % add(name, idx)
    geometries(idx) % name = name

    ! Point geometry
    call move_alloc(geom, geometries(idx) % geom)

  end subroutine addGeom

  !!
  !! Get index of geometry given its name
  !!
  !! Args:
  !!   name [in] -> Name of the geoemtry
  !!
  !! Result:
  !!   Index of the geometry with the name
  !!
  !! Errors:
  !!   fatalError if name was not defined
  !!
  function geomIdx(name) result(idx)
    character(nameLen), intent(in) :: name
    integer(shortInt)              :: idx
    integer(shortInt), parameter   :: NOT_PRESENT = -8
    character(*), parameter :: Here = 'geomIdx (geometryReg_mod.f90)'

    idx = geometryNameMap % getOrDefault(name, NOT_PRESENT)
    if (idx == NOT_PRESENT) then
      call fatalError(Here, 'Geometry: '//trim(name)//' was not found.')
    end if

  end function geomIdx

  !!
  !! Get pointer to a geometry by index
  !!
  !! Args:
  !!   idx [in] -> Index of the requested geometry
  !!
  !! Result:
  !!   Pointer to the geometry with the idx
  !!
  !! Errors:
  !!   FatalError if the index is not valid
  !!
  function geomPtr(idx) result(ptr)
    integer(shortInt), intent(in) :: idx
    class(geometry), pointer      :: ptr
    character(*), parameter :: Here = 'geomPtr (geometryReg_mod.f90)'

    ! Check index
    if (idx < 1 .or. idx > geometryTop) then
      call fatalError(Here,'Index: '//numToChar(idx)//' does not correspond to valid geometry. &
                           &Must be 1-'//numToChar(geometryTop))
    end if

    ptr => geometries(idx) % geom

  end function geomPtr

  !!
  !! Get number of geometries stored in the module
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of allocated geometries
  !!
  function geomNum() result(N)
    integer(shortInt) :: N

    ! Check allocation and get size
    if (allocated(geometries)) then
      N = geometryTop
    else
      N = 0
    end if

  end function geomNum

  !!
  !! Add Field definition
  !!
  !! Args:
  !!   kentta [in] -> Pre-initialised field, to be allocated
  !!   name [in]   -> Name of the field
  !!
  !! Errors:
  !!   fatalError if field with the name was already defined
  !!
  subroutine addField(kentta, name)
    class(field), allocatable, intent(inout) :: kentta
    character(nameLen), intent(in)           :: name
    integer(shortInt)                        :: idx
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter    :: Here = 'addField (geometryReg_mod.f90)'

    ! Get free index
    idx = fieldNameMap % getOrDefault(name, NOT_PRESENT)

    if (idx /= NOT_PRESENT) then
      call fatalError(Here, 'Field with name: '//trim(name)//' has already been defined with &
                            &idx: '//numToChar(idx))
    else
      idx = freeFieldIdx()
      fieldTop = fieldTop + 1 ! Increment counter
    end if

    ! Initialise & store index
    call fieldNameMap % add(name, idx)
    fields(idx) % name = name

    ! Point field
    call move_alloc(kentta, fields(idx) % kentta)

  end subroutine addField
  
  !!
  !! Returns whether a field is present
  !!
  !! Args:
  !!   name [in] -> Name of the field
  !!
  !! Result:
  !!   Logical stating whether the field exists
  !!
  pure function hasField(name) result(exists)
    character(nameLen), intent(in) :: name
    logical(defBool)               :: exists
    integer(shortInt), parameter   :: NOT_PRESENT = -8
    integer(shortInt)              :: idx

    idx = fieldNameMap % getOrDefault(name, NOT_PRESENT)
    exists = (idx /= NOT_PRESENT)

  end function hasField

  !!
  !! Get index of a field given its name
  !!
  !! Args:
  !!   name [in] -> Name of the field
  !!
  !! Result:
  !!   Index of the field with the name
  !!
  !! Errors:
  !!   fatalError if name was not defined
  !!
  function fieldIdx(name) result(idx)
    character(nameLen), intent(in) :: name
    integer(shortInt)              :: idx
    integer(shortInt), parameter   :: NOT_PRESENT = -8
    character(*), parameter :: Here = 'fieldIdx (geometryReg_mod.f90)'

    idx = fieldNameMap % getOrDefault(name, NOT_PRESENT)
    if (idx == NOT_PRESENT) then
      call fatalError(Here, 'Field: '//trim(name)//' was not found.')
    end if

  end function fieldIdx

  !!
  !! Get pointer to a field by index
  !!
  !! Args:
  !!   idx [in] -> Index of the requested field
  !!
  !! Result:
  !!   Pointer to the field with the idx
  !!
  !! Errors:
  !!   FatalError if the index is not valid
  !!
  function fieldPtr(idx) result(ptr)
    integer(shortInt), intent(in) :: idx
    class(field), pointer         :: ptr
    character(*), parameter :: Here = 'fieldPtr (geometryReg_mod.f90)'

    ! Check index
    if (idx < 1 .or. idx > fieldTop) then
      call fatalError(Here,'Index: '//numToChar(idx)//' does not correspond to valid field. &
                           &Must be 1-'//numToChar(fieldTop))
    end if

    ptr => fields(idx) % kentta

  end function fieldPtr

  !!
  !! Return to uninitialised state
  !!
  subroutine kill()
    integer(shortInt) :: i

    ! Kill geometries
    call geometryNameMap % kill()
    if (allocated(geometries)) then
      do i = 1, geometryTop
        call geometries(i) % geom % kill()
      end do

      deallocate(geometries)
    end if
    geometryTop = 0

    ! Kill fields
    call fieldNameMap % kill()
    if (allocated(fields)) then
      do i = 1, fieldTop
        call fields(i) % kentta % kill()
      end do

      deallocate(fields)
    end if
    fieldTop = 0

  end subroutine kill

  !!
  !! Get geometry free Idx
  !!
  !! Manages the size and extension of `geometries` array.
  !! Returns a next free index in the array
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   A free index on the `geometries`
  !!   Does NOT increment the geometryTop counter
  !!
  function freeGeomIdx() result(idx)
    integer(shortInt)                        :: idx
    integer(shortInt)                        :: N, i
    type(geomBox), dimension(:), allocatable :: temp

    ! Set index
    idx = geometryTop + 1

    ! Check allocation/extension of geometries array
    if (allocated(geometries)) then
      N = size(geometries)

      if (idx > N) then ! Extend Array
        N = int(GROWTH_RATE * N)

        ! Allocate temporary
        allocate (temp(N))

        ! Load values
        do i = 1, geometryTop
          temp(i) % name = geometries(i) % name
          call move_alloc(geometries(i) % geom, temp(i) % geom)
        end do

        ! Switch allocation
        call move_alloc(temp, geometries)

      end if

    else ! Allocate array
      N = START_SIZE
      allocate (geometries(N))

    end if

  end function freeGeomIdx

  !!
  !! Get field free Idx
  !!
  !! Manages the size and extension of `fields` array.
  !! Returns a next free index in the array
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   A free index on the `fields`
  !!   Does NOT increment the fieldTop counter
  !!
  function freeFieldIdx() result(idx)
    integer(shortInt)                        :: idx
    integer(shortInt)                        :: N, i
    type(fieldBox), dimension(:), allocatable :: temp

    ! Set index
    idx = fieldTop + 1

    ! Check allocation/extension of geometries array
    if (allocated(fields)) then
      N = size(fields)

      if (idx > N) then ! Extend Array
        N = int(GROWTH_RATE * N)

        ! Allocate temporary
        allocate (temp(N))

        ! Move values
        do i = 1, fieldTop
          temp(i) % name = fields(i) % name
          call move_alloc(fields(i) % kentta, temp(i) % kentta)
        end do

        ! Switch allocation
        call move_alloc(temp, fields)

      end if

    else ! Allocate array
      N = START_SIZE
      allocate (fields(N))

    end if

  end function freeFieldIdx

end module geometryReg_mod
