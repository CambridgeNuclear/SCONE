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
!!   fieldIdx -> Get index of a field
!!   fieldPtr -> Get pointer to a field specified by index
!!   display  -> Display info abou defined fields and geometries
!!   kill     -> Return to uninitialised state
!!
module geometryReg_mod

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use dictionary_class,  only : dictionary
  use charMap_class,     only : charMap

  ! Meterial interface
  use materialMenu_mod,  only : mm_nameMap => nameMap

  ! Geometry
  use geometry_inter,    only : geometry
  use geometryStd_class, only : geometryStd
  use geometryGrid_class, only : geometryGrid

  ! Fields
  use field_inter,              only : field
  use uniformScalarField_class, only : uniformScalarField
  use uniformVectorField_class, only : uniformVectorField

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
  public :: addField
  public :: fieldIdx
  public :: fieldPtr
  public :: kill

  !! Parameters
  character(nameLen), dimension(*), parameter :: AVAILABLE_GEOMETRIES = ['geometryStd ' ,&
                                                                         'geometryGrid']
  character(nameLen), dimension(*), parameter :: AVAILABLE_FIELDS = ['uniformScalarField',&
                                                                     'uniformVectorField']
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
  !!   name [in]   -> Name of the geometry
  !!   dict [in]   -> Dictionary with geometry definition
  !!   silent [in] -> Optional. Set to .true. to surpress console messeges. Default .false.
  !!     Note that errors will still be printed with silent=.true.
  !!
  !! Errors:
  !!   fatalError if geometry with the name was already defined
  !!
  subroutine addGeom(name, dict, silent)
    character(nameLen), intent(in)         :: name
    class(dictionary), intent(in)          :: dict
    logical(defBool), optional, intent(in) :: silent
    integer(shortInt)                      :: idx
    logical(defBool)                       :: silent_l
    integer(shortInt), parameter :: NOT_PRESENT = -7
    character(100), parameter    :: Here = 'addGeom (geometryReg_mod.f90)'

    ! Get silent value
    if (present(silent)) then
      silent_l = silent
    else
      silent_l = .false.
    end if


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
    call new_geometry(geometries(idx) % geom, dict, mm_namemap, silent_l)

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
  !! Add Field definition
  !!
  !! Args:
  !!   name [in]   -> Name of the field
  !!   dict [in]   -> Dictionary with field definition
  !!
  !! Errors:
  !!   fatalError if field with the name was already defined
  !!
  subroutine addField(name, dict)
    character(nameLen), intent(in)         :: name
    class(dictionary), intent(in)          :: dict
    integer(shortInt)                      :: idx
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
    call new_field(fields(idx) % kentta, dict)

  end subroutine addField

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
  !! Allocates and initialises a geometry from dictionary
  !!
  !! This is Factory procedure for geometries
  !!
  !! Args:
  !!   geom [out] -> Allocatable geometry to be allocated
  !!   dict [in]  -> Dictionary with geometry definition
  !!   mats [in]  -> Map of material names to matIdx
  !!   silent [in] -> Optional. Set to .true. to surpress console messeges. Default .false.
  !!     Note that errors will still be printed with silent=.true.
  !!
  !! Errors:
  !!   fatalError is type of geometry is unknown
  !!
  subroutine new_geometry(geom, dict, mats, silent)
    class(geometry), allocatable, intent(out) :: geom
    class(dictionary), intent(in)             :: dict
    type(charMap), intent(in)                 :: mats
    logical(defBool), optional, intent(in)    :: silent
    logical(defBool)                          :: silent_l
    character(nameLen)                        :: type
    character(100), parameter :: Here = 'new_geometry (geometryReg_mod.f90)'

    ! Get type
    call dict % get(type, 'type')

    ! Get silent flag
    if (present(silent)) then
      silent_l = silent
    else
      silent_l = .false.
    end if

    ! Allocate to right type
    select case (type)
      case ('geometryStd')
        allocate(geometryStd :: geom)

      case ('geometryGrid')
        allocate(geometryGrid :: geom)

      case default
        print '(A)', 'AVAILABLE GEOMETRIES'
        print '(A)', AVAILABLE_GEOMETRIES
        call fatalError(Here, trim(type)// ' is not valid geometry. See list above.')

    end select

    ! Initialise geometry
    call geom % init(dict, mats, silent_l)

  end subroutine new_geometry

  !!
  !! Allocates and initialises field from dictionary
  !!
  !! Args:
  !!   kentta [out] -> Field to be allocated (kentta ~= field (Fi.))
  !!   dict [in]   -> Dictionary with definition
  !!
  !! Errors:
  !!   fatalError is type of field is unknown
  !!
  subroutine new_field(kentta, dict)
    class(field), allocatable, intent(out) :: kentta
    class(dictionary), intent(in)          :: dict
    character(nameLen)                     :: type
    character(100), parameter :: Here = 'new_field (geometryReg_mod.f90) '

    ! Get type
    call dict % get(type, 'type')

    ! Build Field
    select case (type)
      case ('uniformScalarField')
        allocate(uniformScalarField :: kentta)

      case ('uniformVectorField')
        allocate(uniformVectorField :: kentta)

      case default
        print '(A)', "AVAILABLE FIELDS:"
        print '(A)', AVAILABLE_FIELDS
        call fatalError(Here, trim(type)//' is not valid field. See list above.')

    end select

    ! Initialise Field
    call kentta % init(dict)

  end subroutine new_field

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
