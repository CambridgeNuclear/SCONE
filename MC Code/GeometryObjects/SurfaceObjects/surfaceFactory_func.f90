module surfaceFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! Diffrent surface classes
  use surface_inter,            only : surface
  use infSurf_class,            only : infSurf
  ! Planes
  use plane_class,              only : plane
  use xPlane_class,             only : xPlane
  use yPlane_class,             only : yPlane
  use zPlane_class,             only : zPlane
  ! Cylinders and spheres
  use sphere_class,             only : sphere
  use xCylinder_class,          only : xCylinder
  use yCylinder_class,          only : yCylinder
  use zCylinder_class,          only : zCylinder
  ! Boxes & Square Cylinders
  use box_class,                only : box
  use xSquareCylinder_class,    only : xSquareCylinder
  use ySquareCylinder_class,    only : ySquareCylinder
  use zSquareCylinder_class,    only : zSquareCylinder
  ! TruncatedCylinders
  use xTruncatedCylinder_class, only : xTruncatedCylinder
  use yTruncatedCylinder_class, only : yTruncatedCylinder
  use zTruncatedCylinder_class, only : zTruncatedCylinder

  implicit none
  private

  public :: new_surface
  public :: new_surface_ptr

  ! *** ADD NAME OF A NEW SURFACE HERE ***!
  ! List that contains all accaptable types of surfaces
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_surfaces = [ 'infSurf        ',&
                                                                     'plane          ', &
                                                                     'xPlane         ', &
                                                                     'yPlane         ', &
                                                                     'zPlane         ', &
                                                                     'sphere         ', &
                                                                     'xCylinder      ', &
                                                                     'yCylinder      ', &
                                                                     'zCylinder      ', &
                                                                     'box            ' ,&
                                                                     'xSquareCylinder',&
                                                                     'ySquareCylinder',&
                                                                     'zSquareCylinder',&
                                                                     'xTruncCylinder ',&
                                                                     'yTruncCylinder ',&
                                                                     'zTruncCylinder ' ]


contains

  !!
  !! Returns allocatable surface form dictionary
  !!
  function new_surface(dict,name) result(new)
    class(dictionary), intent(in) :: dict
    character(nameLen),intent(in) :: name
    class(surface),allocatable    :: new
    character(nameLen)            :: type
    character(100),parameter      :: Here = 'new_surface (surfaceFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of surface
    ! *** ADD CASE STATEMENT FOR A NEW SURFACE BELOW ***!
    select case(type)
      case('infSurf')
        allocate(new, source = infSurf(dict,name) )

      case('plane')
        allocate(new, source = plane(dict,name) )

      case('xPlane')
        allocate(new, source = xPlane(dict,name) )

      case('yPlane')
        allocate(new, source = yPlane(dict,name) )

      case('zPlane')
        allocate(new, source = zPlane(dict,name) )

      case('sphere')
        allocate(new, source = sphere(dict,name) )

      case('xCylinder')
        allocate(new, source = xCylinder(dict,name) )

      case('yCylinder')
        allocate(new, source = yCylinder(dict,name) )

      case('zCylinder')
        allocate(new, source = zCylinder(dict,name) )

      case('box')
        allocate(new, source = box(dict,name) )

      case('xSquareCylinder')
        allocate(new, source = xSquareCylinder(dict,name) )

      case('ySquareCylinder')
        allocate(new, source = ySquareCylinder(dict,name) )

      case('zSquareCylinder')
        allocate(new, source = zSquareCylinder(dict,name) )

      case('xTruncCylinder')
        allocate(new, source = xTruncatedCylinder(dict,name) )

      case('yTruncCylinder')
        allocate(new, source = yTruncatedCylinder(dict,name) )

      case('zTruncCylinder')
        allocate(new, source = zTruncatedCylinder(dict,name) )

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
  function new_surface_ptr(dict,name) result(new)
    class(dictionary), intent(in) :: dict
    character(nameLen),intent(in) :: name
    class(surface),pointer        :: new

    ! Allocate pointer and copy data from local allocatable
    allocate( new, source = new_surface(dict,name) )

  end function new_surface_ptr



end module surfaceFactory_func
