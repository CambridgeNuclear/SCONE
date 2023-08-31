!!
!! Module to store temperature and/or frequency-dependent equations for materials, especially
!! those too complicated to be easily read in from an input file. Also contains an energy grid to
!! allow materials to access particle energy group bounds for use in evaluating these equations.
!!
!! Also stores energy grid for multigroup problems for easy access by materials
!!
!! For a new set of material equations:
!!   -> Add name to AVAILABLE_equations
!!   -> Add case to evaluateCv and evaluateSigma
!!   -> Evaluate simple equations (e.g. 'marshak' or 'hohlraum') in these functions,
!!      or can link to new functions (e.g. 'olson1D')
!!
module materialEquations

  use numPrecision
  use universalVariables
  use genericProcedures,       only : fatalError
  use energyGrid_class,        only : energyGrid

  implicit none
  private

  character(nameLen),dimension(*),parameter :: AVAILABLE_equations = ['marshak ',&
                                                                      'hohlraum',&
                                                                      'olson1D ']

  public :: evaluateCv
  public :: evaluateSigma
  public :: normPlanckSpectrum

  ! Energy grid for multi-frequency problems for easy access by material classes
  type(energyGrid), pointer, public :: mgEnergyGrid => null()

  interface evaluateCv
    module procedure evaluateCv
  end interface

  interface evaluateSigma
    module procedure evaluateSigma
  end interface


  contains

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Non-specific interface equations - add new case to each
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Evaluate heat capacity equation at a given temperature
  !!
  !! Args:
  !!   equation -> name of the equation to use
  !!   T        -> temperature to evaluate at
  !!
  function evaluateCv(equation, T) result(cv)
    character(*), intent(in)        :: equation
    real(defReal), intent(in)       :: T
    real(defReal)                   :: cv
    character(100), parameter       :: Here = 'getCv (materialEquations.f90)'

    select case(equation)

      case('marshak')
        cv = 7.14

      case('hohlraum')
        cv = 0.3

      case('olson1D')
        cv = cvOlson1D(T)

      case default
        cv = ZERO
        print *, AVAILABLE_equations
        call fatalError(Here, 'Unrecognised equation: '//trim(equation))

    end select

  end function evaluateCv

  !!
  !! Evaluate opacity equation at a given temperature and frequency (energy E)
  !!
  !! Args:
  !!   equation -> name of the equation to use
  !!   T        -> temperature to evaluate at
  !!   nu       -> frequency of the particle
  !!
  function evaluateSigma(equation, T, E) result(sigma)
    character(*), intent(in)        :: equation
    real(defReal), intent(in)       :: T
    real(defReal), intent(in)       :: E
    real(defReal)                   :: sigma
    character(100), parameter       :: Here = 'getSigma (materialEquations.f90)'

    select case(equation)

      case('marshak')
        sigma = 10*T**(-3)

      case('hohlraum')
        sigma = 100*T**(-3)

      case('olson1D')
        sigma = sigmaOlson1D(T, E)

      case default
        sigma = ZERO
        print *, AVAILABLE_equations
        call fatalError(Here, 'Unrecognised equation: '//trim(equation))

    end select

  end function evaluateSigma

!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Custom material equations for various input files
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! 1-Dimensional Multi-frequency test by Olson (2020)
  !!
  !! Olson, G.L., 2020. Stretched and filtered multigroup pn transport for improved positivity
  !! and accuracy. Journal of Computational and Theoretical Transport 0, 1–18.
  !!
  function cvOlson1D(T) result(cv)
    real(defReal), intent(in) :: T
    real(defReal)             :: cv, root, alpha, dAlphadT

    root = sqrt(1+4*exp(0.1/T))
    alpha = 0.5*exp(-0.1/T)*(root-1)
    dAlphadT = 0.1*(alpha-1/root)/(T*T)

    cv = 0.1*(1+alpha+(T+0.1)*dAlphadT)

    ! Deal with numerical errors from poorly defined regions (e.g. T almost 0)
    if (cv /= cv .or. cv > INF) cv = ZERO

  end function cvOlson1D

  function sigmaOlson1D(T, E) result(sigma)
    real(defReal), intent(in)     :: T
    real(defReal), intent(in)     :: E
    real(defReal)                 :: sigma

    if (E < 0.008) then
      sigma = min(1e7_defReal, 1e9_defReal*T*T)
    else if (E < 0.3) then
      sigma = 192/(E*E*(1+200*T**1.5))
    else
      sigma = 192*sqrt(0.3)/(E**2.5*(1+200*T**1.5)) + 4e4_defReal*(0.3/E)**2.5/(1+8000*T**2)
    end if

    ! Multiply by density
    sigma = 0.001*sigma

  end function sigmaOlson1D


!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Commonly used equations
!!<><><><><><><>><><><><><><><><><><><>><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Evaluate frequency-normalised Planck spectrum
  !!
  !! Args:
  !!   nu -> frequency
  !!   T  -> temperature
  !!
  pure function normPlanckSpectrum(E, T) result(b)
    real(defReal), intent(in) :: E
    real(defReal), intent(in) :: T
    real(defReal)             :: b

    b = 15*E**3 / ((pi*T)**4 * (exp(E/T)-1))

  end function normPlanckSpectrum


end module materialEquations
