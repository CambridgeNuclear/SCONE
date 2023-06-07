!!
!! Module to build new fields. Called by geometryReg_mod
!!
module geometryFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use charMap_class,     only : charMap

  ! Geometry
  use geometry_inter,    only : geometry
  use geometryStd_class, only : geometryStd
  use geometryReg_mod,   only : gr_addGeom => addGeom

  ! Meterial interface
  use materialMenu_mod,  only : mm_nameMap => nameMap

  implicit none
  private


  !! Parameters
  character(nameLen), dimension(*), parameter :: AVAILABLE_GEOMETRIES = ['geometryStd']

  ! Public interface
  public :: new_geometry

contains

  !!
  !! Allocates and initialises a geometry from dictionary
  !!
  !! This is Factory procedure for geometries
  !!
  !! Args:
  !!   dict [in]  -> Dictionary with geometry definition
  !!   name [in]  -> Name of the geometry for the geometry registry
  !!   silent [in] -> Optional. Set to .true. to surpress console messeges. Default .false.
  !!     Note that errors will still be printed with silent=.true.
  !!
  !! Errors:
  !!   fatalError is type of geometry is unknown
  !!
  subroutine new_geometry(dict, name, silent)
    class(dictionary), intent(in)          :: dict
    character(nameLen), intent(in)         :: name
    logical(defBool), optional, intent(in) :: silent
    class(geometry), allocatable           :: geom
    logical(defBool)                       :: silent_l
    character(nameLen)                     :: type
    character(100), parameter :: Here = 'new_geometry (geometryFactory_func.f90)'

    ! Get silent flag
    if (present(silent)) then
      silent_l = silent
    else
      silent_l = .false.
    end if

    ! Get type
    call dict % get(type, 'type')

    ! Allocate to right type
    select case (type)
      case ('geometryStd')
        allocate(geometryStd :: geom)

      case default
        print '(A)', 'AVAILABLE GEOMETRIES'
        print '(A)', AVAILABLE_GEOMETRIES
        call fatalError(Here, trim(type)// ' is not valid geometry. See list above.')

    end select

    ! Initialise geometry
    call geom % init(dict, mm_nameMap, silent_l)

    ! Call geometry registry to add geometry
    call gr_addGeom(geom, name)

  end subroutine new_geometry


end module geometryFactory_func
