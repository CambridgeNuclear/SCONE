module noAngle_class

  use numPrecision
  use angleLawENDF_class, only : angleLawENDF
  use RNG_class, only : RNG

  implicit none
  private

  interface noAngle
    module procedure new_noAngle
  end interface

  type, public,extends(angleLawENDF) :: noAngle
    private
  contains
    procedure :: sample
    procedure :: probabilityOf
  end type noAngle

contains


  function sample(self,E,rand) result (miu)
    class(noAngle), intent(in)   :: self
    real(defReal), intent(in)    :: E
    class(RNG), intent(inout)    :: rand
    real(defReal)                :: miu

    miu = 1.0

  end function sample

  function probabilityOf(self,miu,E) result (prob)
    class(noAngle), intent(in)  :: self
    real(defReal), intent(in)   :: E, miu
    real(defReal)               :: prob

    if (miu == 1.0_defReal) then
      prob = 1.0
    else
      prob = 0.0
    end if

  end function probabilityOf

  function new_noAngle()
    type(noAngle), pointer  :: new_noAngle

    allocate(new_noAngle)

  end function

    
end module noAngle_class
