module legendrePoly_func
  !! This module contains functions related to evaluating and sampling Legendre Polynomials
  !! All of the methods here are based on:
  !!   "Monte Carlo Particle Transport Methods: Neutron and Photon Calculations"
  !!    by I. Lux and L. Koblinger
  !!    ISBN 978-1-315-89573
  !!

  use numPrecision
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG

  implicit none
  private

  ! Public Interface
  public :: sampleLegendre

  ! Generic Interfaces

  !!
  !! Samples diffrent Legendre PDFs given number of coefficients and RNG
  !!   Diffrent Orders can use diffrent methods and have some limits on supported coefficients
  !!   Consult particular implementations for details
  !!
  interface sampleLegendre
    module procedure sampleLegendre_P1
  end interface sampleLegendre

contains

  !!
  !! Simple function to sample P1 Legendre PDF. Uses methods from Lux & Koblinger
  !!
  function sampleLegendre_P1(P1,rand) result(x)
    real(defReal), intent(in)    :: P1
    class(RNG), intent(inout)    :: rand
    real(defReal)                :: x
    real(defReal)                :: P1_loc
    real(defReal)                :: threshold
    integer(shortInt)            :: Low, Top, exec
    integer(shortInt), parameter :: UNIFORM = 1, LIN = 2, DELTA = 3
    character(100), parameter :: Here = 'sampleLegendre_P1 ( legendrePoly_func.f90)'

    ! Make local copy of P1 coeff. Take abs() to simplify code
    ! -ve P1 will be inverted at the end.
    P1_loc = abs(P1)

    ! Depending on whether P1 > 1 determine treshold and associated PDF for the mixing method
    ! For further details refer to Lux and Koblinger APPENDIX 3D
    ! If random number < threshold then Top is used.
    if ( P1_loc < ONE) then
      threshold = P1_loc
      Top = LIN
      Low = UNIFORM

    else if( P1_loc <= 3.0_defReal) then
      threshold = 0.5 * (P1_loc - ONE)
      Top = DELTA
      Low = LIN

    else
      call fatalError(Here,'P1 must have absolute value < 3.0')
      ! Avoid warnings
      threshold = ONE
      Top = 0
      Low = 0

    end if

    ! Use mixing method with the calculated Threshold
    if ( rand % get() < threshold ) then
      exec = Top
    else
      exec = Low
    end if

    ! Sample from UNIFORM ( PDF = 0.5); LIN ( PDF = 0.5 + 0.5 *mu) or DELATA ( PDF = DELTA(mu-1))
    select case(exec)
      case (UNIFORM)
        x = TWO * rand % get() - ONE

      case (LIN)
        ! Need to solve CDF(x) = 0.25 * x^2 + 0.5 * x + 0.25 = (0.5*x+0.5)^2)
        x = TWO * sqrt(rand % get()) - ONE

      case (DELTA)
        x = ONE

      case default
        call fatalError(Here,'This should never happen. WTF?')
        x = ZERO

    end select

    ! Invert result if P1 is -ve
    if ( P1 < ZERO ) x = -x

  end function sampleLegendre_P1

    
end module legendrePoly_func
