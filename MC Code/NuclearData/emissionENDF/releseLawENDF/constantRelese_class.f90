module constantRelese_class

  use numPrecision
  use releseLawENDF_class, only : releseLawENDF

  implicit none
  private

  interface constantRelese
   module procedure new_constantRelese
  end interface

  type, public, extends(releseLawENDF) :: constantRelese
      private
      real(defReal) :: secondaryRelese = 1.0
    contains
      procedure :: init
      procedure :: releseAt
  end type constantRelese

contains

  subroutine init(self,relese)
    class(constantRelese), intent(inout) :: self
    real(defReal), intent(in)            :: relese

    self % secondaryRelese = relese

  end subroutine init

  function releseAt(self,energy) result(relese)
    class(constantRelese), intent(in) :: self
    real(defReal), intent(in)         :: energy
    real(defReal)                     :: relese

    relese = self % secondaryRelese

  end function releseAt

  function new_constantRelese(relese) result(newConstantRelese)
    real(defReal), intent(in)           :: relese
    type(constantRelese),pointer        :: newConstantRelese

    allocate(newConstantRelese)
    call newConstantRelese % init(relese)

  end function new_constantRelese

    
end module constantRelese_class
