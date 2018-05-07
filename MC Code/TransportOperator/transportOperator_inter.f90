module transportOperator_inter
  use numPrecision
  use genericProcedures
  use universalVariables

  use particle_class
  use surface_class
  use geometry_class
  use cell_class
  use rng_class

  !use Nuclear_Data_MG_class           !Re-instate later!!!
  !use Nuclear_Data_CE_class

  implicit none
  private

  type, abstract, public :: transportOperator
    class(rng), pointer      :: random => null()         ! RNG - should this be associated to particle?
    !class(transport_Nuclear_data), pointer :: nuclearData => null()  ! nuclear data
    class(geometry), pointer :: geom => null()           ! references the geometry for cell searching
  contains
    procedure(transport), deferred :: transport
    procedure(applyBC), deferred   :: applyBC
  end type transportOperator

  abstract interface

    subroutine transport(self, p)
      import :: transportOperator, &
                particle
      class(transportOperator), intent(in) :: self
      class(particle), intent(inout)       :: p
    end subroutine transport

    subroutine applyBC(self, p, currentCell)
      import :: transportOperator, &
                particle, &
                cell_ptr
      class(transportOperator), intent(in) :: self
      class(particle), intent(inout)        :: p
      class(cell_ptr), intent(inout)        :: currentCell
    end subroutine applyBC

  end interface

end module transportOperator_inter

