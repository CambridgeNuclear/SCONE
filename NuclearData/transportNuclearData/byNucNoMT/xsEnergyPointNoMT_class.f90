module xsEnergyPointNoMT_class

  use numPrecision
  use xsMainCDF_class,   only : xsMainCDF
  use xsMainSet_class,   only : xsMainSet

  implicit none
  private

  !!
  !! Wrapper type to store all relevant information about nuclide reaction cross-sections at
  !! a single energy.
  !!
  type, public :: xsEnergyPointNoMT
    type(xsMainCDF) :: mainCDF     ! CDF to select main reaction channel
    type(xsMainSet) :: xs          ! reaction cross sections
  contains
    procedure :: interpolate
    procedure :: interpolateTotal
    procedure :: interpolateTail
    procedure :: load
  end type xsEnergyPointNoMT

contains

  !!
  !! Interpolate between two xsEnergyPoints (NoMT) using provided interpolation factor
  !! (x-x_0)/(x_1-x_0) where 1-top; 0-bottom. No error checks are performed so care must be taken.
  !!
  subroutine interpolate(self,low,top,f)
    class(xsEnergyPointNoMT), intent(inout) :: self
    type(xsEnergyPointNoMT), intent(in)     :: low
    type(xsEnergyPointNoMT), intent(in)     :: top
    real(defReal), intent(in)               :: f

    ! Interpolate CDF
    call self % mainCDF % interpolate ( low % mainCDF, &
                                        top % mainCDF, &
                                        f )
    ! Interpolate cross-sections
    call self % xs % interpolate( low % xs, &
                                  top % xs, &
                                  f )

  end subroutine interpolate

  !!
  !! Interpolate only total cross-section
  !!
  subroutine interpolateTotal(self,low,top,f)
    class(xsEnergyPointNoMT), intent(inout) :: self
    type(xsEnergyPointNoMT), intent(in)     :: low
    type(xsEnergyPointNoMT), intent(in)     :: top
    real(defReal), intent(in)               :: f

    call self % xs % interpolateTotal( low % xs, &
                                       top % xs, &
                                        f )
  end subroutine interpolateTotal

  !!
  !! Interpolate cdf and all cross-sections EXECPT TOTAL CROSS-SECTION
  !!
  subroutine interpolateTail(self,low,top,f)
    class(xsEnergyPointNoMT), intent(inout) :: self
    type(xsEnergyPointNoMT), intent(in)     :: low
    type(xsEnergyPointNoMT), intent(in)     :: top
    real(defReal), intent(in)               :: f

    ! Interpolate CDF
    call self % mainCDF % interpolate ( low % mainCDF, &
                                        top % mainCDF, &
                                        f )
    ! Interpolate cross-sections
    call self % xs % interpolateTail( low % xs, &
                                      top % xs, &
                                      f )

  end subroutine interpolateTail


  !!
  !! Load data into xsEnergyPointNoMT by provideing all cross-section individualy.
  !! It is not called init to avoid name conflict with its child class nuclideMemory, which
  !! needs to be initialised with diffrent data.
  !!
  elemental subroutine load(self,scatter,capture,fission)
    class(xsEnergyPointNoMT), intent(inout) :: self
    real(defReal), intent(in)               :: scatter
    real(defReal), intent(in)               :: capture
    real(defReal), intent(in),optional      :: fission
    real(defReal)                           :: total
    real(defReal)                           :: locFission

    ! Set local fission to 0.0 if fission xs is not given
    if (present(fission)) then
      locFission = fission

    else
      locFission = 0.0

    end if

    ! Calculate total xs
    total = scatter + capture + locFission

    ! Load xs and initialise cdf
    self % xs % total   = total
    self % xs % scatter = scatter
    self % xs % capture = capture
    self % xs % fission = locFission

    call self % mainCDF % init(scatter,capture,locFission)

  end subroutine load

    
end module xsEnergyPointNoMT_class
