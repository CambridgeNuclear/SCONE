module constantRelease_class

  use numPrecision
  use genericProcedures,    only : fatalError
  use releaseLawENDF_inter, only : releaseLawENDF

  implicit none
  private

  interface constantRelease
   module procedure new_constantRelease
  end interface

  !!
  !! Constant release of neutrons independent of incedent energy
  !!
  type, public, extends(releaseLawENDF) :: constantRelease
      private
      real(defReal) :: secondaryRelease = ONE
    contains
      procedure :: init
      procedure :: releaseAt
      procedure :: kill
  end type constantRelease

contains

  !!
  !! Initialise
  !!
  subroutine init(self,release)
    class(constantRelease), intent(inout) :: self
    real(defReal), intent(in)             :: release
    character(100),parameter              :: Here='init (constantRelease_class.f90)'

    if( release < 0) call fatalError(Here,'-ve value of release provided!')
    self % secondaryRelease = release

  end subroutine init

  !!
  !! Release at energy E_in
  !!
  function releaseAt(self,E_in) result(release)
    class(constantRelease), intent(in) :: self
    real(defReal), intent(in)          :: E_in
    real(defReal)                      :: release

    release = self % secondaryRelease

  end function releaseAt

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(constantRelease), intent(inout) :: self

    self % secondaryRelease = ONE

  end subroutine kill

  !!
  !! Constructor
  !!
  function new_constantRelease(release) result(new)
    real(defReal), intent(in)            :: release
    type(constantRelease)                :: new

    call new % init(release)

  end function new_constantRelease


end module constantRelease_class
