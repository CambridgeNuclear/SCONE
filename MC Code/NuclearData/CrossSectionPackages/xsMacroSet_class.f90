module xsMacroSet_class

  use numPrecision

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
  end type xsMacroSet

contains


    
end module xsMacroSet_class
