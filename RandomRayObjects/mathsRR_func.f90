module mathsRR_func
  
  !!
  !! This module contains maths used in random ray.
  !! First it has a function for efficiently computing an exponential
  !! for use in MoC implementations using a rational approximation.
  !! I believe this originates from the M&C 2019 publication:
  !! "Adding a third level of parallelism to OpenMOC"
  !!
  !! TODO: Add functions for efficient spherical harmonics evaluation.
  !!

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private
  
  public :: expF1, expF1Tau, expG, expG2

  ! Numerator coefficients in F1 rational approximation
  real(defFlt), parameter :: c1n = -1.0000013559236386308, c2n = 0.23151368626911062025,&
          c3n = -0.061481916409314966140, c4n = 0.0098619906458127653020, c5n = -0.0012629460503540849940, &
          c6n = 0.00010360973791574984608, c7n = -0.000013276571933735820960

  ! Denominator coefficients in F1 rational approximation
  real(defFlt), parameter :: c0d = 1.0_defFlt, c1d = -0.73151337729389001396, c2d = 0.26058381273536471371, &
          c3d = -0.059892419041316836940, c4d = 0.0099070188241094279067, c5d = -0.0012623388962473160860, &
          c6d = 0.00010361277635498731388, c7d = -0.000013276569500666698498

  ! Numerator coefficients in G rational approximation
  real(defFlt), parameter :: d0n = 0.5, d1n = 0.176558112351595, d2n = 0.04041584305811143, &
          d3n = 0.006178333902037397, d4n = 0.0006429894635552992 , d5n = 0.00006064409107557148

  ! Denominator coefficients in G rational approximation 
  real(defFlt), parameter :: d0d = 1.0, d1d = 0.6864462055546078, d2d = 0.2263358514260129, &
          d3d = 0.04721469893686252, d4d = 0.006883236664917246, d5d = 0.0007036272419147752 , d6d = 0.00006064409107557148 

  ! Coefficients for numerator in G2 rational approximation
  real(defFlt), parameter :: g1n = -0.08335775885589858, g2n = -0.003603942303847604, &
        g3n = 0.0037673183263550827, g4n = 0.00001124183494990467, g5n = 0.00016837426505799449

  ! Coefficients for denominator in G2 rational approximation
  real(defFlt), parameter :: g1d = 0.7454048371823628, g2d = 0.23794300531408347, &
        g3d = 0.05367250964303789, g4d = 0.006125197988351906, g5d = 0.0010102514456857377

contains

  !!
  !! Computes x = [1 - exp(-tau)]/tau for use in MoC calcs
  !! Tau is the optical distance.
  !! F1 is a common name in MoC literature
  !!
  elemental function expF1(tau) result(x)
    real(defFlt), intent(in)    :: tau
    real(defFlt)                :: x
    real(defFlt)                :: den, num

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
    ! Reintroduce this to give 1-exp(-tau)
    !num = num * x
    !x = num / den

    x = -num / den

  end function expF1
  
  !!
  !! Computes x = [1 - exp(-tau)] = F1*tau for use in MoC calcs
  !! Tau is the optical distance.
  !! F1 is a common name in MoC literature
  !!
  elemental function expF1Tau(tau) result(x)
    real(defFlt), intent(in)    :: tau
    real(defFlt)                :: x
    real(defFlt)                :: den, num

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

  end function expF1Tau

  !!
  !! Computes y = 1/x-(1-exp(-x))/x**2 using a 5/6th order rational approximation.
  !! From OpenMOC.
  !! Commonly referred to as G in MoC literature. Used to compute other exponentials.
  !!
  elemental function expG(tau) result(x)
    real(defFlt), intent(in)    :: tau
    real(defFlt)                :: x
    real(defFlt)                :: den, num

    x = tau

    den = d6d * x + d5d
    den = den * x + d4d
    den = den * x + d3d
    den = den * x + d2d
    den = den * x + d1d
    den = den * x + d0d

    num = d5n * x + d4n
    num = num * x + d3n
    num = num * x + d2n
    num = num * x + d1n
    num = num * x + d0n

    x = num / den

  end function expG

  !!
  !! Computes y = 2/3 - (1 + 2/x) * (1/x + 0.5 - (1 + 1/x) * (1-exp(-x)) / x)
  !! using a 5/5th order rational approximation,
  !! From OpenMoC.
  !!
  elemental function expG2(tau) result(x)
    real(defFlt), intent(in)    :: tau
    real(defFlt)                :: x
    real(defFlt)                :: den, num

    x = tau
    
    num = g5n*x + g4n
    num = num*x + g3n
    num = num*x + g2n
    num = num*x + g1n
    num = num*x

    ! Calculate denominator
    den = g5d*x + g4d
    den = den*x + g3d
    den = den*x + g2d
    den = den*x + g1d
    den = den*x + 1.0_defFlt

    x = num / den

  end function expG2
    
end module mathsRR_func
