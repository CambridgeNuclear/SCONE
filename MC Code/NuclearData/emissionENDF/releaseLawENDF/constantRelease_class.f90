module constantRelease_class

  use numPrecision
  use releaseLawENDF_class, only : releaseLawENDF

  implicit none
  private

  interface constantRelease
   module procedure new_constantrelease
  end interface

  type, public, extends(releaseLawENDF) :: constantRelease
      private
      real(defReal) :: secondaryRelease = 1.0
    contains
      procedure :: init
      procedure :: releaseAt
  end type constantRelease

contains

  subroutine init(self,release)
    class(constantrelease), intent(inout) :: self
    real(defReal), intent(in)            :: release

    self % secondaryrelease = release

  end subroutine init

  function releaseAt(self,energy) result(release)
    class(constantrelease), intent(in) :: self
    real(defReal), intent(in)         :: energy
    real(defReal)                     :: release

    release = self % secondaryrelease

  end function releaseAt

  function new_constantRelease(release) result(newConstantRelease)
    real(defReal), intent(in)           :: release
    type(constantrelease),pointer        :: newConstantRelease

    allocate(newConstantRelease)
    call newConstantRelease % init(release)

  end function new_constantRelease

    
end module constantRelease_class
