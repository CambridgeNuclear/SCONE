module xsMacroSet_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError

  implicit none
  private

  !!
  !! Package of material macroscopic XSs. Will be used by tallies as well.
  !!
  type, public :: xsMacroSet
    real(defReal) :: totalXS   = 0.0
    real(defReal) :: scatterXS = 0.0
    real(defReal) :: captureXS = 0.0
    real(defReal) :: fissionXS = 0.0
  contains
    procedure :: invert

  end type xsMacroSet

contains

  !!
  !! Sample reaction using the xs in the set
  !!
  function invert(self,r) result(MT)
    class(xsMacroSet), intent(in) :: self
    real(defReal), intent(in)     :: r
    integer(shortInt)             :: MT
    real(defReal)                 :: r_scaled
    integer(shortInt)             :: i
    character(100),parameter      :: Here ='invert (xsMacroSet_class.f90)'

    ! Scale r to total reaction cross-section
    r_scaled = r * self % totalXS
    i=0

    ! Protect against -ve r
    if(r_scaled >= 0) i=i+1

    ! Decrement through all reactions
    r_scaled = r_scaled - self % scatterXS
    if(r_scaled > 0) i=i+1

    r_scaled = r_scaled - self % captureXS
    if(r_scaled > 0) i=i+1

    r_scaled = r_scaled - self % fissionXS
    if(r_scaled > 0) i=i+1

    ! Assign MT
    select case(i)
      case(1)
        ! Scattering
        MT = macroAllScatter
      case(2)
        ! Capture
        MT = macroCapture
      case(3)
        ! Fission
        MT = macroFission
      case(0)
        ! r < 0
        call fatalError(Here, 'Provided number to invert is -ve')

      case default
        ! r > 1 or wrong total xs
        if (r > 1.0_defReal) then
          call fatalError(Here,'Provided number to invert is greater then 1')

        else
          call fatalError(Here,'Total cross section must be too large')

        end if
    end select

  end function invert


    
end module xsMacroSet_class
