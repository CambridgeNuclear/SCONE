module matNucCDF_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  !!
  !! Object to store and transfer CDF of collision with a given nuclide in a material.
  !! It is cheating becouse it really stores just individual and total reaction cross sections.
  !! It is not fully safe becouse each cross-section needs to be loaded manually
  !!
  type, public :: matNucCDF
    integer(shortInt), dimension(:), allocatable, private :: nucIndices
    real(defReal), dimension(:), allocatable              :: nucTotalXS
    real(defReal)                                         :: matTotalXS
  contains
    procedure :: invert
    procedure :: init
  end type matNucCDF

contains

  !!
  !! Allocate space and load nuclide indices
  !!
  subroutine init(self,nucIndices)
    class(matNucCDF), intent(inout)            :: self
    integer(shortInt),dimension(:),intent(in)  :: nucIndices

    if(allocated(self % nucIndices)) deallocate (self % nucIndices)
    if(allocated(self % nucTotalXS)) deallocate (self % nucTotalXS)

    self % nucIndices = nucIndices
    allocate(self % nucTotalXS (size(nucIndices)))

  end subroutine init

  !!
  !! Invert probability distribution
  !!
  function invert(self,r) result(nucIdx)
    class(matNucCDF), intent(in) :: self
    real(defReal), intent(in)    :: r
    integer(shortInt)            :: nucIdx
    real(defReal)                :: r_scaled
    character(100),parameter     :: Here = 'invert (matNucCDF_class.f90)'
    integer(shortInt)            :: i

    ! Scale r to total reaction cross-section
    r_scaled = r * self % matTotalXS

    ! Search nucTotalXSs (scaled PDF) with linear search
    do i = 1,size(self % nucTotalXS)

      ! Decrement scaled random number
      r_scaled = r_scaled - self % nucTotalXS(i)

      if( r_scaled <= 0) then
        nucIdx = self % nucIndices(i)
        return

      end if
    end do

    call fatalError(Here,'Provided number to invert cdf must be < 1 ')

    ! Avoid compiler warning (Serpent territory?...)
    NucIdx = 0

  end function invert
    
end module matNucCDF_class
