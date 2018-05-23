module keffClerk_inter

  use numPrecision
  use tallyClerk_inter, only : tallyClerk

  implicit none
  private

  type, public,extends(tallyClerk),abstract :: keffClerk
    private
  contains
    procedure(keff),deferred :: keff
  end type keffClerk

  abstract interface
    !!
    !! Returns an estimate of keff
    !!
    function keff(self) result(k)
      import :: keffClerk,&
                defReal
      class(keffClerk), intent(in) :: self
      real(defReal)                :: k
    end function keff

  end interface

end module keffClerk_inter
