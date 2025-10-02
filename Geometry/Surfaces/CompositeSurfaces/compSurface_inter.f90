module compSurface_inter

  use surface_inter, only : surface

  implicit none
  private

  !!
  !! Abstract interface to group all composite surfaces.
  !!
  !! At the moment it adds nothing new, except grouping the composite
  !! surfaces into a familly of subclasses of surface.
  !!
  type, public, abstract, extends(surface) :: compSurface
    ! Adds nothing yet
  end type

end module compSurface_inter
