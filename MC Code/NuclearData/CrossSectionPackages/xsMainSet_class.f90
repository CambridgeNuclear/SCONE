module xsMainSet_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError
  use xsSet_class,       only : xsSet

  implicit none
  private

  type, public,extends(xsSet) :: xsMainSet
    ! xss of main(lumped) reaction channels
    real(defReal)   :: total   = 0.0  ! Total xs
    real(defReal)   :: scatter = 0.0  ! anyScatter
    real(defReal)   :: capture = 0.0  ! anyCapture
    real(defReal)   :: fission = 0.0  ! anyFission

  contains
    procedure :: xsOf
    procedure :: invert
    procedure :: interpolate
    procedure :: interpolateTotal
    procedure :: interpolateTail
  end type xsMainSet

contains
  !!
  !! Obtain xs value gicven MT number. Give error if MT is invalid.
  !!
  function xsOf(self,MT) result(xs)
    class(xsMainSet), intent(in)  :: self
    integer(shortInt), intent(in) :: MT
    real(defReal)                 :: xs
    character(100),parameter      :: Here ='xsOf (xsMainSet_class.f90)'

    select case(MT)
      case(N_total)
        xs = self % total

      case(anyScatter)
        xs = self % scatter

      case(anyCapture)
        xs = self % capture

      case(anyFission)
        xs = self % fission

      case default
        call fatalError(Here,'Provided MT number does not match any reaction in the set')

    end select
    
  end function xsOf

  !!
  !! Sample reaction using the xs in the set
  !!
  function invert(self,r) result(MT)
    class(xsMainSet), intent(in) :: self
    real(defReal), intent(in)    :: r
    integer(shortInt)            :: MT
    real(defReal)                :: r_scaled
    integer(shortInt)            :: i
    character(100),parameter     :: Here ='invert (xsMainSet_class.f90)'

    ! Scale r to total reaction cross-section
    r_scaled = r * self % total
    i=0

    ! Protect against -ve r
    if(R_scaled >= 0) i=i+1

    ! Decrement through all reactions
    r_scaled = r_scaled - self % scatter
    if(r_scaled > 0) i=i+1

    r_scaled = r_scaled - self % capture
    if(r_scaled > 0) i=i+1

    r_scaled = r_scaled - self % fission
    if(r_scaled > 0) i=i+1

    ! Assign MT
    select case(i)
      case(1)
        ! Scattering
        MT = anyScatter
      case(2)
        ! Capture
        MT = anyCapture
      case(3)
        ! Fission
        MT = anyFission
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



  !!
  !! Perform linear interpolation between 2 xsMainSets (low and top) using the given
  !! interpolation factor (x-x_0)/(x_1-x_0) where 1- top value; 0- low value.
  !! Does not check sensibility of the output! May result in -ve xs!
  !!
  subroutine interpolate(self,low,top,f)
    class(xsMainSet), intent(inout)  :: self
    type(xsMainSet), intent(in)      :: low
    type(xsMainSet), intent(in)      :: top
    real(defReal), intent(in)        :: f
    real(defReal)                    :: f2

    ! Calculate (1-factor)
    f2 = 1.0 - f

    ! Interpolate xss
    self % total   = f2 * low % total   + f * top % total
    self % scatter = f2 * low % scatter + f * top % scatter
    self % capture = f2 * low % capture + f * top % capture
    self % fission = f2 * low % fission + f * top % fission

  end subroutine interpolate

  !!
  !! Perform linear interpolation of TOTAL XS ONLY between 2 xsMainSets (low and top) using
  !! the given interpolation factor (x-x_0)/(x_1-x_0) where 1- top value; 0- low value.
  !! Does not check sensibility of the output! May result in -ve xs!
  !!
  subroutine interpolateTotal(self,low,top,f)
    class(xsMainSet), intent(inout)  :: self
    type(xsMainSet), intent(in)      :: low
    type(xsMainSet), intent(in)      :: top
    real(defReal), intent(in)        :: f
    real(defReal)                    :: f2

    ! Calculate (1-factor)
    f2 = 1.0 - f

    ! Interpolate xss
    self % total   = f2 * low % total   + f * top % total

  end subroutine interpolateTotal

  !!
  !! Perform linear interpolation of ALL XS EXCEPT TOTAL between 2 xsMainSets (low and top) using
  !! the given interpolation factor (x-x_0)/(x_1-x_0) where 1- top value; 0- low value.
  !! Does not check sensibility of the output! May result in -ve xs!
  !!
  subroutine interpolateTail(self,low,top,f)
    class(xsMainSet), intent(inout)  :: self
    type(xsMainSet), intent(in)      :: low
    type(xsMainSet), intent(in)      :: top
    real(defReal), intent(in)        :: f
    real(defReal)                    :: f2

    ! Calculate (1-factor)
    f2 = 1.0 - f

    ! Interpolate xss
    self % scatter = f2 * low % scatter + f * top % scatter
    self % capture = f2 * low % capture + f * top % capture
    self % fission = f2 * low % fission + f * top % fission

  end subroutine interpolateTail


end module xsMainSet_class
