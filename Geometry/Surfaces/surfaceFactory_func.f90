module surfaceFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! Surface interface
  use surface_inter,     only : surface

  ! Surfaces
  use aPlane_class,         only : aPlane
  use cylinder_class,       only : cylinder
  use plane_class,          only : plane
  use sphere_class,         only : sphere
  use box_class,            only : box
  use squareCylinder_class, only : squareCylinder
  use truncCylinder_class,  only : truncCylinder

  implicit none
  private

  ! List that contains all accaptable types of surfaces
  ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_SURFACE = ['xPlane         ',&
                                                                      'yPlane         ',&
                                                                      'zPlane         ',&
                                                                      'plane          ',&
                                                                      'xCylinder      ',&
                                                                      'yCylinder      ',&
                                                                      'zCylinder      ',&
                                                                      'sphere         ',&
                                                                      'box            ',&
                                                                      'xSquareCylinder',&
                                                                      'ySquareCylinder',&
                                                                      'zSquareCylinder',&
                                                                      'xTruncCylinder ',&
                                                                      'yTruncCylinder ',&
                                                                      'zTruncCylinder ' ]

  ! Public interface
  public :: new_surface_ptr
  public :: new_surface


contains

  !!
  !! Return a pointer to a new instance of allocated surface
  !!
  !! Args:
  !!   dict [in] -> Dictionary with surface definition
  !!
  !! Result:
  !!   class(surface) pointer to an allocated instance of the surface
  !!
  !! Errors:
  !!   fatalError if type of surface is unknown
  !!
  function new_surface_ptr(dict) result(new)
    class(dictionary), intent(in) :: dict
    class(surface), pointer       :: new
    character(nameLen)            :: type
    character(100), parameter :: Here = 'new_surface_ptr (surfaceFactory_func.f90)'

    ! Obtain type of the surface
    call dict % get(type, 'type')

    ! Allocate approperiate subclass
    select case (type)
      case ('xPlane', 'yPlane', 'zPlane')
        allocate (aPlane :: new)

      case ('plane')
        allocate (plane :: new)

      case ('sphere')
        allocate (sphere :: new)

      case ('xCylinder', 'yCylinder', 'zCylinder')
        allocate (cylinder :: new)

      case ('box')
        allocate (box :: new)

      case ('xSquareCylinder', 'ySquareCylinder', 'zSquareCylinder')
        allocate (squareCylinder :: new)

      case ('xTruncCylinder', 'yTruncCylinder', 'zTruncCylinder')
        allocate (truncCylinder :: new)

      case default
        print '(A)' , ' AVAILABLE SURFACES: '
        print '(A)' , AVAILABLE_SURFACE
        call fatalError(Here, 'Unrecognised type of a surface: '//trim(type))

      end select

      ! Initialise surface
      call new % init(dict)

  end function new_surface_ptr

  !!
  !! Allocate an allocatable instance of a surface
  !!
  !! Args:
  !!   new [out] -> Surface to be allocated
  !!   dict [in] -> Dictionary with the surface definition
  !!
  !! Errors:
  !!   fatalError if type of surface is unrecognised
  !!
  subroutine new_surface(new, dict)
    class(surface), allocatable, intent(out) :: new
    class(dictionary), intent(in)            :: dict
    class(surface), pointer                  :: temp

    temp => new_surface_ptr(dict)
    allocate (new, source = temp)
    call temp % kill()
    deallocate(temp)

  end subroutine new_surface

end module surfaceFactory_func
