module releaseLawENDF_inter

  use numPrecision

  implicit none
  private
  !!
  !! Abstract interface class for all polymorfic classes to contain various ENDF laws for secondary
  !! neutron emissions
  !!
  type,abstract, public :: releaseLawENDF
      private
    contains
      procedure(releaseAt),deferred :: releaseAt
  end type releaseLawENDF

  abstract interface

    !!
    !! Obtain average neutron emission for incedent energy E_in
    !!
    function releaseAt(self,E_in) result(release)
      import :: defReal,&
                releaseLawENDF
      class(releaseLawENDF), intent(in)  :: self
      real(defReal), intent(in)          :: E_in
      real(defReal)                      :: release
    end function releaseAt

  end interface

end module releaseLawENDF_inter
