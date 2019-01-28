module transportOperator_inter

  use numPrecision
  use universalVariables

  use particle_class,             only : particle
  use dictionary_class,           only : dictionary

  ! Geometry interfaces
  use cellGeometry_inter,         only : cellGeometry

  ! Tally interface
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearData_inter,          only : nuclearData
  use transportNuclearData_inter, only : transportNuclearData


  implicit none
  private

  !!
  !! Abstract interface for all collision operators
  !!
  type, abstract, public :: transportOperator
    !! Components below are not part of the interface!
    !! They are public ONLY to allow inheritance!
    class(transportNuclearData), pointer :: nuclearData => null()  ! nuclear data
    class(cellGeometry), pointer         :: geom        => null()  ! references the geometry for cell searching
    type(tallyAdmin), pointer            :: tally       => null()  ! Tally to recive reports
  contains
    procedure(transport), deferred :: transport
    procedure(init),deferred       :: init
  end type transportOperator

  abstract interface

    !!
    !! Move particle from collision to collision
    !!
    subroutine transport(self, p)
      import :: transportOperator, &
                particle
      class(transportOperator), intent(in) :: self
      class(particle), intent(inout)       :: p
    end subroutine transport

    !!
    !! Connect transport operator to other blocks.
    !! Each implementation checks if the geometry & nuclear data it gets it supports
    !!
    subroutine init(self,nucData,geom,settings)
      import :: transportOperator, &
                nuclearData, &
                cellGeometry, &
                dictionary
      class(transportOperator), intent(inout) :: self
      class(nuclearData),pointer, intent(in)  :: nucData
      class(cellGeometry),pointer, intent(in) :: geom
      class(dictionary),optional,intent(in)   :: settings
    end subroutine init

  end interface

end module transportOperator_inter

