module xsMacroSet_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError, numToChar

  implicit none
  private

  !!
  !! Package of material macroscopic XSs. Will be used by tallies as well.
  !!
  type, public :: xsMacroSet
    real(defReal) :: totalXS     = 0.0
    real(defReal) :: scatterXS   = 0.0
    real(defReal) :: captureXS   = 0.0
    real(defReal) :: fissionXS   = 0.0
    real(defReal) :: nuFissionXS = 0.0
  contains
    procedure :: invert
    procedure :: xsOf
    procedure :: dummy

  end type xsMacroSet

  !!
  !! Pointer wrapper for a set to allow safe access through nuclearData interface
  !!
  type, public :: xsMacroSet_ptr
    private
    type(xsMacroSet), pointer :: ptr => null()
  contains
    ! Assigmnets
    generic   :: assignment(=) => shallowCopy, shallowCopy_pointer
    procedure :: shallowCopy
    procedure :: shallowCopy_pointer

    ! Procedures to Access Target
    procedure :: totalXS     => totalXS_ptr
    procedure :: scatterXS   => scatterXS_ptr
    procedure :: captureXS   => captureXS_ptr
    procedure :: fissionXS   => fissionXS_ptr
    procedure :: nuFissionXS => nuFissionXS_ptr
    procedure :: nu          => nu_ptr
    procedure :: invert      => invert_ptr
    procedure :: xsOf        => xsOf_ptr
  end type xsMacroSet_ptr

contains

  subroutine dummy(self)
    class(xsMacroSet), intent(inout) :: self
  end subroutine dummy

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
        MT = huge(MT)

      case default
        ! r > 1 or wrong total xs
        if (r > 1.0_defReal) then
          call fatalError(Here,'Provided number to invert is greater then 1')

        else
          call fatalError(Here,'Total cross section must be too large')

        end if
        MT = huge(MT)

    end select
  end function invert

  !!
  !! Return XS based on macroscopic MT number
  !!
  function xsOf(self,MT) result(xs)
    class(xsMacroSet), intent(in) :: self
    integer(shortInt), intent(in) :: MT
    real(defReal)                 :: xs
    character(100), parameter     :: Here ='xsOf (xsMacroSet_class.f90)'

    select case(MT)
      case(macroTotal)
        xs = self % totalXS

      case(macroCapture)
        xs = self % captureXS

      case(macroEscatter)
        call fatalError(Here,'Current design of data does not provide macroscopic Elastic scatter')
        xs = ZERO

      case(macroFission)
        xs = self % fissionXS

      case(macroNuFission)
        xs = self % nuFissionXS

      case(macroAbsorbtion)
        xs = self % fissionXS + self % captureXS

      case default
        call fatalError(Here,'Unknown macroscopic MT number: '//numToChar(MT))
        xs = ZERO

    end select
  end function xsOf

  !!************************************************************************************************
  !! Pointer Wrapper Procedures
  !!
  !!************************************************************************************************

  !!
  !! Assignment Pointer Wrapper to Pointer Wrapper
  !!
  subroutine shallowCopy(LHS,RHS)
    class(xsMacroSet_ptr), intent(inout) :: LHS
    type(xsMacroSet_ptr), intent(in)     :: RHS

    LHS % ptr => RHS % ptr

  end subroutine shallowCopy

  !!
  !! Assignmnet Pointer to Pointer Wrapper
  !!
  subroutine shallowCopy_pointer(LHS,RHS)
    class(xsMacroSet_ptr), intent(inout) :: LHS
    type(xsMacroSet),intent(in), pointer :: RHS

    LHS % ptr => RHS

  end subroutine shallowCopy_pointer

  !!
  !! Access Total Macroscopic XS
  !!
  function totalXS_ptr(self) result(xs)
    class(xsMacroSet_ptr), intent(in) :: self
    real(defReal)                     :: xs

    xs = self % ptr % totalXS

  end function totalXS_ptr

  !!
  !! Access Scattering Macroscopic XS. All scattering.
  !!
  function scatterXS_ptr(self) result(xs)
    class(xsMacroSet_ptr), intent(in) :: self
    real(defReal)                     :: xs

    xs = self % ptr % scatterXS

  end function scatterXS_ptr

  !!
  !! Access Capture Macroscopic XS
  !!
  function captureXS_ptr(self) result(xs)
    class(xsMacroSet_ptr), intent(in) :: self
    real(defReal)                     :: xs

    xs = self % ptr % captureXS

  end function captureXS_ptr


  !!
  !! Access Fission Macroscopic XS
  !!
  function fissionXS_ptr(self) result(xs)
    class(xsMacroSet_ptr), intent(in) :: self
    real(defReal)                     :: xs

    xs = self % ptr % fissionXS

  end function fissionXS_ptr

  !!
  !! Access nu*Fission Macroscopic XS
  !!
  function nuFissionXS_ptr(self) result(xs)
    class(xsMacroSet_ptr), intent(in) :: self
    real(defReal)                     :: xs

    xs = self % ptr % nuFissionXS

  end function nuFissionXS_ptr

  !!
  !! Obtain material average nu
  !!
  function nu_ptr(self) result(nu)
    class(xsMacroSet_ptr), intent(in) :: self
    real(defReal)                     :: nu

    if ( self % ptr % fissionXS == ZERO) then
      nu = ZERO

    else
      nu = self % ptr % nuFissionXS / self % ptr % fissionXS

    end if
  end function nu_ptr


  !!
  !! Access Invert procedure
  !!
  function invert_ptr(self,r) result(MT)
    class(xsMacroSet_ptr), intent(in) :: self
    real(defReal), intent(in)         :: r
    integer(shortInt)                 :: MT

    MT = self % ptr % invert(r)

  end function invert_ptr

  !!
  !! Acess xsOf procedure
  !!
  function xsOf_ptr(self,MT) result(xs)
    class(xsMacroSet_ptr), intent(in) :: self
    integer(shortInt), intent(in)     :: MT
    real(defReal)                     :: xs

    xs = self % ptr % xsOf(MT)

  end function xsOf_ptr
    
end module xsMacroSet_class
