module quadSurface_inter

  use surface_inter, only : surface

  implicit none
  private

  !!
  !! Abstract interface to group all quadratic surfaces.
  !!
  !! At the moment it adds nothing new, except grouping the quadratic
  !! surfaces into a familly of subclasses of surface.
  !!
  type, public, abstract, extends(surface) :: quadSurface
    ! Adds nothing yet
  end type

end module quadSurface_inter
