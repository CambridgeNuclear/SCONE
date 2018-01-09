module constantRelease_class

  use numPrecision
  use releaseLawENDF_class, only : releaseLawENDF

  implicit none
  private

  interface constantRelease
   module procedure new_constantRelease
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
    class(constantRelease), intent(inout) :: self
    real(defReal), intent(in)            :: release

    self % secondaryRelease = release

  end subroutine init

  function releaseAt(self,energy) result(release)
    class(constantRelease), intent(in) :: self
    real(defReal), intent(in)         :: energy
    real(defReal)                     :: release

    release = self % secondaryRelease

  end function releaseAt

  function new_constantRelease(release) result(newConstantRelease)
    real(defReal), intent(in)           :: release
    type(constantRelease),pointer        :: newConstantRelease

    allocate(newConstantRelease)
    call newConstantRelease % init(release)

  end function new_constantRelease

    
end module constantRelease_class
