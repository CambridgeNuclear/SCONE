module releaseMatrixMG_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError

  implicit none
  private


  !!
  !! Stores average emissions for MG material
  !! Memory inefficient implementation
  !!
  type, public :: releaseMatrixMG
    private
    real(defReal), dimension(:,:),allocatable :: scatterRelease
    real(defReal), dimension(:),allocatable   :: fissionRelease
  contains
    procedure :: init
    procedure :: releaseForScattering
    procedure :: releaseForFission
  end type releaseMatrixMG

contains

  !!
  !! Initialise
  !! For production matrix indexing convection is (G_out,G_in)
  !!
  subroutine init(self,nuVector, scatterProd)
    class(releaseMatrixMG), intent(inout)     :: self
    real(defReal), dimension(:), intent(in)   :: nuVector
    real(defReal), dimension(:,:), intent(in) :: scatterProd ! Scattering production matrix
    integer(shortInt)                         :: nG
    character(100), parameter                 :: Here = 'init (releaseMatrixMG_class.f90)'


    ! Read number of groups and deallocate arrays
    nG = size(nuVector)

    if(allocated(self % scatterRelease)) deallocate(self % scatterRelease)
    if(allocated(self % fissionRelease)) deallocate(self % fissionRelease)


    ! Check for invalid input
    if( nG /= size(scatterProd,1)) then
      call fatalError(Here, 'Provided nuVector and production matrix have diffrent group size')

    else if (size(scatterProd,1) /= size(scatterProd,2)) then
      call fatalError(Here, 'Scattering production amtrix is not square')

    else if (any(nuVector < 0.0 )) then
      call fatalError(Here, 'nuVector contains -ve entries')

    else if (any(scatterProd < 0.0)) then
      call fatalError(Here, 'Scattering production matrix contains -ve entries')

    end if

    ! Load data into storage
    self % fissionRelease = nuVector
    self % scatterRelease = scatterProd

  end subroutine init

  !!
  !! Get release for scattering from G_in to G_out
  !! Function interface is provided for compability with future sparse implementation.
  !! Currently most data is redundant and stores a lot of 1.0
  !!
  function releaseForScattering(self,G_in,G_out) result(nu)
    class(releaseMatrixMG), intent(in) :: self
    integer(shortInt), intent(in)      :: G_in
    integer(shortInt), intent(in)      :: G_out
    real(defReal)                      :: nu

    nu = self % scatterRelease(G_out,G_in)

  end function releaseForScattering

  !!
  !! Get release for fission in group G_in
  !!
  function releaseForFission(self,G_in) result(nu)
    class(releaseMatrixMG), intent(in) :: self
    integer(shortInt), intent(in)      :: G_in
    real(defReal)                      :: nu

    nu = self % fissionRelease(G_in)

  end function releaseForFission

    

end module releaseMatrixMG_class
