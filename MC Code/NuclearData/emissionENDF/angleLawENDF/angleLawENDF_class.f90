module angleLawENDF_class

  use numPrecision
  use miuEndfPdf_class, only : miuEndfPdf_ptr

  implicit none
  private

  type,public :: angleLawENDF
    !! Class to contain emission angle propabilities for secondary particles
      private
      real(defReal),dimension(:),allocatable        :: energyGrid
      type(miuEndfPdf_ptr),dimension(:),allocatable :: miuEndfPdfs
    contains
      procedure :: init
  end type angleLawENDF

contains

  subroutine init(self,energyGrid,miuEndfPdfs)
    class(angleLawENDF),intent(inout)             :: self
    real(defReal),dimension(:), intent(in)        :: energyGrid
    type(miuEndfPdf_ptr),dimension(:), intent(in) :: miuEndfPdfs

    if(allocated(self % energyGrid))  deallocate(self % energyGrid)
    if(allocated(self % miuEndfPdfs)) deallocate(self % miuEndfPdfs)

    self % energyGrid  = energyGrid
    self % miuEndfPdfs = miuEndfPdfs
  end subroutine
    
end module angleLawENDF_class
