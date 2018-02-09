module maxwellEnergyPdf_class
  !! Bundle of procedures to obtain samples and probabilities from maxwellian distribution
  !! of energy.
  !!
  !! For Johnk's Theorem please refer to:
  !! Devroye, L., 1986. Non-Uniform Random Variate Generation. 1st ed. New York: Springer-Verlag.
  !!

  use numPrecision
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG

  implicit none
  private

  type, public :: maxwellEnergyPdf
    private
  contains
    procedure :: sample        => sample_Johnk
    generic   :: probabilityOf => probabilityOf_full ,&
                                  probabilityOf_withC

    procedure,private :: probabilityOf_full
    procedure,private :: probabilityOf_withC

  end type maxwellEnergyPdf
contains

  function sample_Johnk(self,kT,rand) result (E)
  !! Samples Maxwellian energy distribution ( gamma(3/2,1) distribution) using algorithms
  !! based on Johnk theorem. This is the MCNP, OpenMC and Serpent(?) approach.
    class(maxwellEnergyPdf), intent(in) :: self
    real(defReal),intent(in)            :: kT
    class(RNG), intent(inout)           :: rand
    real(defReal)                       :: E
    real(defReal)                       :: r1, r2, r3
    real(defReal)                       :: beta, gamma05, cosine

    ! Obtain all random numbers
    r1 = rand % get()
    r2 = rand % get()
    r3 = rand % get()

    ! Define B(a,b) as a RANDOM VARIABLE governed by beta(a,b) distribution
    ! Similarly define G(a,b) as a RANDOM VARIABLE governed by gamma(a,b) distribution
    ! Note that B(a,b)*G(a,b) is a product of RANDOM VARIABLES not of Probability Density Functions

    ! Obtain sample of B(0.5,0.5) distribution based on Johnk's Theorem. Sample of cosine
    ! instead of using rejection scheme

    cosine = cos(0.5*PI*r1)
    beta = cosine * cosine

    ! Obtain sample of G(0.5,1) using the fact that G(0.5,1) = B(0.5,0.5) * G(1,1)
    ! G(1,1) is just exponential distribution

    gamma05 = -log(r2) * beta

    ! Obtain sample of G(3/2,1) using the facte that G(3/2,1) = G(0.5,1) + G(1,1)

    E = (-log(r3) + gamma05) * kT

  end function sample_Johnk



  function probabilityOf_full(self,E_out,kT) result(prob)
  !! Returns probability of E_out.
    class(maxwellEnergyPdf), intent(in) :: self
    real(defReal), intent(in)           :: E_out
    real(defReal), intent(in)           :: kT
    real(defReal)                       :: prob
    real(defReal)                       :: Inv_Const
    real(defReal)                       :: funct

    funct = sqrt(E_out) * exp(-E_out/kT)

    ! Const^-1
    Inv_Const = kT * sqrt(kT) * 0.5 * sqrt(PI)

    ! funct * Const
    prob = funct / Inv_Const

  end function probabilityOf_full



  function probabilityOf_withC(self,E_out,kT,C) result(prob)
  !! Returns probability of E_out. Normalisation constant C is provided by the user.
  !! This function should be used when maximum energy is limited with restriction energy
  !! and PDF does not extends to infinity but reaches 0 at some other value
    class(maxwellEnergyPdf), intent(in) :: self
    real(defReal), intent(in)           :: E_out
    real(defReal), intent(in)           :: kT
    real(defReal), intent(in)           :: C
    real(defReal)                       :: prob

   prob = C * sqrt(E_out) * exp(-E_out/kT)

  end function probabilityOf_withC



end module maxwellEnergyPdf_class
