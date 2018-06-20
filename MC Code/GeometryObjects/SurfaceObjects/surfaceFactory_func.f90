module surfaceFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! Diffrent surface classes
  use surface_inter, only : surface

  implicit none
  private

  ! *** ADD NAME OF A NEW SURFACE HERE ***!
  ! List that contains all accaptable types of surfaces
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_surfaces = [ 'keffActiveClerk  ',&
                                                                     'keffInactiveClerk' ]


contains

  !!
  !! Returns allocatable surface form dictionary
  !!
  function new_surface(dict) result(new)
    class(dictionary), intent(in) :: dict
    class(surface),allocatable    :: new
    character(nameLen)            :: type
    character(100),parameter      :: Here = 'new_surface (surfaceFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of surface
    ! *** ADD CASE STATEMENT FOR A NEW SURFACE BELOW ***!
    select case(type)
      case('keffActiveClerk')
      !  allocate(new, source = keffActiveClerk(dict) )

      case('keffInactiveClerk')
       ! allocate(new, source = keffInactiveClerk(dict) )

     !*** NEW SURFACE TEMPLATE ***!
     !case('<newSUrfaceName>')
     !  allocate(new, source = <newSurfaceName>(dict) )
     !
      case default
        print *, AVALIBLE_surfaces
        call fatalError(Here, 'Unrecognised type of surface: ' // trim(type))

    end select


  end function new_surface

  !!
  !! Returns pointer to allocated surface from dictionary
  !!
  function new_surface_ptr(dict) result(new)
    class(dictionary), intent(in) :: dict
    class(surface),pointer        :: new

    ! Allocate pointer and copy data from local allocatable
    allocate( new, source = new_surface(dict) )

  end function new_surface_ptr



end module surfaceFactory_func
