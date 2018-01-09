module releaseLawENDF_class

  use numPrecision

  implicit none
  private

  type,abstract, public :: releaseLawENDF
    !! Abstract interface class for all polymorfic classes to contain various ENDF laws for secondary
    !! neutron emissions
      private
    contains
      procedure(releaseAt),deferred :: releaseAt
  end type releaseLawENDF

abstract interface

  function releaseAt(self,energy) result(release)
    import :: defReal,&
              releaseLawENDF
    class(releaseLawENDF), intent(in)  :: self
    real(defReal), intent(in)          :: energy
    real(defReal)                      :: release
  end function releaseAt

end interface

contains
    
end module releaseLawENDF_class
