module sourceFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! source interface
  use source_inter,      only : source

  ! source implementations
  use pointSource_class,     only : pointSource
  use fissionSource_class,   only : fissionSource
  use materialSource_class,  only : materialSource
  use bbSurfaceSource_class, only : bbSurfaceSource

  ! geometry
  use geometry_inter,    only : geometry

  implicit none
  private

  public :: new_source

  ! *** ADD NAME OF A NEW SOURCE HERE ***!
  ! List that contains all accaptable types of sources
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen),dimension(*),parameter :: AVAILABLE_sources = [ 'pointSource    ',&
                                                                     'fissionSource  ',&
                                                                     'materialSource ',&
                                                                     'bbsurfaceSource']

contains

  !!
  !! Allocate new allocatable source to a specific type
  !! If new is allocated it deallocates it
  !!
  subroutine new_source(new, dict, geom)
    class(source),allocatable, intent(inout) :: new
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    character(nameLen)                       :: type
    character(100),parameter :: Here = 'new_source (sourceFactory_func.f90)'

    ! Deallocate new if allocated
    if(allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of source
    ! *** ADD CASE STATEMENT FOR A NEW SOURCE BELOW ***!
    select case(type)
      case('pointSource')
        allocate(pointSource :: new)
        call new % init(dict, geom)

      case('fissionSource')
        allocate(fissionSource :: new)
        call new % init(dict, geom)

      case('materialSource')
        allocate(materialSource :: new)
        call new % init(dict, geom)

      case('bbSurfaceSource')
        allocate(bbSurfaceSource :: new)
        call new % init(dict, geom)

     !*** NEW SOURCE TEMPLATE ***!
     !case('<newSourceName>')
     !  allocate(<newSourceName> :: new)
     !  call new % init(dict, geom)
     !
     case default
       print *, AVAILABLE_sources
       call fatalError(Here, 'Unrecognised type of source: ' // trim(type))

    end select

  end subroutine new_source

end module sourceFactory_func
