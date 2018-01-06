module releseLawENDF_class

  use numPrecision

  implicit none
  private

  type,abstract, public :: releseLawENDF
    !! Abstract interface class for all polymorfic classes to contain various ENDF laws for secondary
    !! neutron emissions
      private
    contains
      procedure(releseAt),deferred :: releseAt
  end type releseLawENDF

abstract interface

  function releseAt(self,energy) result(relese)
    import :: defReal,&
              releseLawENDF
    class(releseLawENDF), intent(in)  ::  self
    real(defReal), intent(in)         :: energy
    real(defReal)                     :: relese
  end function releseAt

end interface

contains
    
end module releseLawENDF_class
