module xsMainSet_class

  use numPrecision
  use endfConstants
  use xsSet_class, only : xsSet

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

  subroutine interpolateTotal

  !!
  !! Perform linear interpolation of ALL XS EXCEPT TOTAL between 2 xsMainSets (low and top) using
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
    self % scatter = f2 * low % scatter + f * top % scatter
    self % capture = f2 * low % capture + f * top % capture
    self % fission = f2 * low % fission + f * top % fission

  subroutine interpolateTotal


end module xsMainSet_class
