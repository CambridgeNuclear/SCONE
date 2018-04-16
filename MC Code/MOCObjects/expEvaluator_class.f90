!!
!! Evaluates the exponential used in MOC
!! Can be either a regular evaluation or based on a lookup table
!!
module expEvaluator_class

  use numPrecision
  use universalVariables
  use genericProcedures
  use lookup_class

  implicit none
  private

  type, public                 :: expEvaluator
    logical(defBool)           :: tabular
    class(lookup), private     :: table
    real(defReal)              :: maxOpticalLength
    real(defReal)              :: maxError
  contains
    procedure         :: init => initIntrinsic, initTabular
    procedure         :: eval
  end type expEvaluator

contains

  !!
  !! Constructor for default exponential evaluator
  !!
  subroutine initIntrinsic(self)
    class(expEvaluator), intent(inout) :: self
    self % tabular = .FALSE.

  end subroutine initIntrinsic

  !!
  !! Constructor for look-up exponential evaluator
  !!
  subroutine initTabular(self,maxOpticalLength,maxError)
    class(expEvaluator), intent(inout) :: self
    real(defReal), intent(in) :: maxError
    real(defReal), intent(in) :: maxOpticalLength
    self % tabular = .TRUE.

  end subroutine initTabular

  !!
  !! Evaluate the exponential 1 - exp(-tau)
  !!
  function eval(self,opticalLength)result(expVal)
    class(expEvaluator), intent(in) :: self
    real(defReal), intent(in) :: opticalLength
    real(defReal), intent(out) :: expVal

    if (self % tabular) then

    else
      expVal = 1 - exp(-opticalLength)
      return
    end if

  end function eval

end module expEvaluator_class
