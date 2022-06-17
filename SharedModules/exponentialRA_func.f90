module exponentialRA_func
  !! This module contains a function for efficiently computing an exponential
  !! for use in MoC implementations. RA stands for rational approximation.
  !! This is based on the implementation given in Minray:
  !! github.com/jtramm/minray/blob/master/cpu_srce/flux_attenuation_kernel.c
  !! I believe this originates from the M&C 2019 publication:
  !! "Adding a third level of parallelism to OpenMOC"
  !!

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private
  
  public :: exponential

  ! Numerator coefficients in rational approximation
  real(defReal), parameter :: c1n = -1.0000013559236386308, c2n = 0.23151368626911062025,&
          c3n = -0.061481916409314966140, c4n = 0.0098619906458127653020, c5n = -0.0012629460503540849940, &
          c6n = 0.00010360973791574984608, c7n = -0.000013276571933735820960

  ! Denominator coefficients in rational approximation
  real(defReal), parameter :: c0d = ONE, c1d = -0.73151337729389001396, c2d = 0.26058381273536471371, &
          c3d = -0.059892419041316836940, c4d = 0.0099070188241094279067, c5d = -0.0012623388962473160860, &
          c6d = 0.00010361277635498731388, c7d = -0.000013276569500666698498

contains

  !!
  !! Computes x = 1 - exp(-tau) for use in MoC calcs
  !! Tau is the optical distance
  !!
  elemental function exponential(tau) result(x)
    real(defReal), intent(in)    :: tau
    real(defReal)                :: x
    real(defReal)                :: den, num

    x = -tau
    den = c7d
    den = den * x + c6d
    den = den * x + c5d
    den = den * x + c4d
    den = den * x + c3d
    den = den * x + c2d
    den = den * x + c1d
    den = den * x + c0d

    num = c7n
    num = num * x + c6n
    num = num * x + c5n
    num = num * x + c4n
    num = num * x + c3n
    num = num * x + c2n
    num = num * x + c1n
    num = num * x

    x = num / den

  end function exponential

    
end module exponentialRA_func
