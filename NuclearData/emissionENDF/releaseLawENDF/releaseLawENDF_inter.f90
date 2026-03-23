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
    procedure(hasEnergy),deferred :: hasEnergy
    procedure(kill),deferred      :: kill
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

    !!
    !! Return whether a given energy is available with the law
    !!
    function hasEnergy(self, E_in) result(has)
      import :: defReal, &
                defBool, &
                releaseLawENDF
      class(releaseLawENDF), intent(in)  :: self
      real(defReal), intent(in)          :: E_in
      logical(defBool)                   :: has
    end function hasEnergy

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: releaseLawENDF
      class(releaseLawENDF), intent(inout) :: self
    end subroutine kill

  end interface

end module releaseLawENDF_inter
