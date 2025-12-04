module sourceFactory_func

  use numPrecision
  use errors_mod,        only : fatalError
  use dictionary_class,  only : dictionary

  ! source interface
  use source_inter,      only : source

  ! source implementations
  use pointSource_class,    only : pointSource
  use fissionSource_class,  only : fissionSource
  use materialSource_class, only : materialSource
  use mixSource_class,      only : mixSource, mixSource_TptrCast
  use externalSource_class, only : externalSource

  ! geometry
  use geometry_inter,    only : geometry

  implicit none
  private

  public :: new_source

  ! List that contains all accaptable types of sources
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all entries have the same length
  character(nameLen), dimension(*), parameter :: AVAILABLE_sources = [ 'pointSource   ',&
                                                                       'fissionSource ',&
                                                                       'materialSource',&
                                                                       'mixSource     ',&
                                                                       'externalSource' ]

contains

  !!
  !! Allocate new allocatable source to a specific type
  !! If new is allocated it deallocates it
  !!
  recursive subroutine new_source(new, dict, geom)
    class(source), allocatable, intent(inout) :: new
    class(dictionary), intent(in)             :: dict
    class(geometry), pointer, intent(in)      :: geom
    character(nameLen)                        :: type
    integer(shortInt)                         :: i
    type(mixSource), pointer                  :: mixSource_ptr
    character(nameLen), dimension(:), allocatable :: sourceNames
    character(100),parameter :: Here = 'new_source (sourceFactory_func.f90)'

    ! Deallocate new if allocated
    if (allocated(new)) deallocate(new)

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of source
    select case(type)
      case('pointSource')
        allocate(pointSource :: new)

      case('fissionSource')
        allocate(fissionSource :: new)

      case('materialSource')
        allocate(materialSource :: new)
      
      case('externalSource')
        allocate(externalSource :: new)

      case('mixSource')
        allocate(mixSource :: new)

        mixSource_ptr => mixSource_TptrCast(new)
        if (.not. associated(mixSource_ptr)) call fatalError(Here,'Failed to get mixSource pointer')

        ! Build mixture of sources here
        ! This is done to avoid circular dependences with mixSource_class
        call dict % get(sourceNames, 'sources')
        allocate(mixSource_ptr % sources(size(sourceNames)))

        do i = 1, size(sourceNames)
          call new_source(mixSource_ptr % sources(i) % sourceItem, dict % getDictPtr(sourceNames(i)), geom)
        end do

     case default
       print *, AVAILABLE_sources
       call fatalError(Here, 'Unrecognised type of source: ' // trim(type))

    end select

    ! Initialise new source
    call new % init(dict, geom)

  end subroutine new_source

end module sourceFactory_func
