!!
!! Module to build geometries and store information (name + pointer) about them
!!
module geometryFactory_func

  use numPrecision
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary

  ! Nuclear data interface
  use nuclearData_inter,  only : nuclearData

  ! Abstract interfaces
  use geometry_inter,     only : geometry
  use cellGeometry_inter, only : cellGeometry

  ! Individual implementations
  use basicCellCSG_class, only : basicCellCSG

  implicit none
  private

  ! *** ADD NAME OF A NEW GEOMETRY HERE ***!
  ! List that contains all accaptable types of geometry
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_geometry = [ 'basicCellCSG']


  !!
  !! Public module interface
  !!
  public :: new_geometry_ptr
  public :: new_cellGeometry_ptr

contains

  !!
  !! Builds and allocates an instance of geometry from dictionary
  !!
  function new_geometry_ptr(dict, materials) result(new)
    class(dictionary), intent(in)         :: dict
    class(nuclearData), intent(in)        :: materials
    class(geometry), pointer              :: new
    character(nameLen)                    :: type
    character(100),parameter :: Here = 'new_geometry_ptr (geometryFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of geometry
    ! *** ADD CASE STATEMENT FOR A NEW GEOMETRY BELOW ***!
    ! **** AT THE MOMENT ALLOCATE + SELECT TYPE + INIT is very unelegant implementation
    ! **** Will have to be improved
    select case(type)
      case('basicCellCSG')
        ! Allocate and initialise
        allocate( basicCellCSG :: new)
        select type(new)
          type is (basicCellCSG)
            call new % init(dict,materials)
        end select

      case default
        print *, AVALIBLE_geometry
        call fatalError(Here, 'Unrecognised type of geometry: ' // trim(type))

    end select

  end function new_geometry_ptr

  !!
  !! Builds and allocates an instance of geometry from dictionary
  !!
  function new_cellGeometry_ptr(dict, materials) result(new)
    class(dictionary), intent(in)         :: dict
    class(nuclearData), intent(in)        :: materials
    class(cellGeometry), pointer          :: new
    class(geometry),pointer               :: temp
    character(100),parameter :: Here = 'new_cellGeometry_ptr (geometryFactory_func.f90)'

    temp => new_geometry_ptr(dict, materials)

    select type(temp)
      class is( cellGeometry)
        new => temp

      class default
        call fatalError(Here,'Requested geometry type is not cellGeometry')
    end select

  end function new_cellGeometry_ptr

end module geometryFactory_func
