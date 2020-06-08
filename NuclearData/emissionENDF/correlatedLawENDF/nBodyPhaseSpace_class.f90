module nBodyPhaseSpace_class

  use numPrecision
  use genericProcedures,       only : fatalError, numToChar
  use RNG_class,               only : RNG
  use aceCard_class,           only : aceCard
  use maxwellEnergyPdf_class,  only : maxwellEnergyPdf
  use correlatedLawENDF_inter, only : correlatedLawENDF

  implicit none
  private

  !!
  !! Constructors
  !!
  interface nBodyPhaseSpace
    module procedure new_nBodyPhaseSpace_fromACE
    module procedure new_nBodyPhaseSpace_fromValues
  end interface

  !!
  !! N-body phase space distribution; ACE Law 66, ENDF Energy-Angle Distribution LAW=6
  !! Mu is isotropic
  !! Energy is given by Beta distribution with shape dependent on number of outgoing
  !! particles N in {3,4,5}. See ENDF Manual[1] for additional details.
  !!
  !! Outgoing energy is given by (see [1] for definition of Emax) :
  !! P(E') = C_n * sqrt(E') * (E_max - E') ^ (3*n/2 -4)
  !!
  !! Making substitution E' = T * E_max for T in [0,1] we can obtain
  !! P(T) = C_n' * sqrt(T) * (1 - T)^(3*n/2 -4)
  !!
  !! Which can be readily indentified as a beta distribution with following parameters depending
  !! on the value of n:
  !! n = 2 : P(T) = beta(1.5, 1.5)
  !! n = 3 : P(T) = beta(1.5, 3 )
  !! n = 4 : P(T) = beta(1.5, 4.5)
  !!
  !! Using the fact that beta distrubution can be represented by a ratio of gamma distributions
  !! as follows[2]:
  !!
  !! beta(A,B) = Ga /(Ga + Gb) with Ga = gamma(A,1); Gb = gamma(B,1)
  !!
  !! We can sample outgoing energy using standard algorithms for sampling gamma distributions.
  !! Note that gamma(1.5) is a Maxwellian energy distribution with temperature of 1.
  !! We note that gamma(N*1/2,1) distribution can be represented by a sum of gamma(1,1)
  !! distribution and a single beta(0.5, 0.5) distribution[2] e.g.:
  !!
  !! gamma(2.5,1) =  gamma(1,1) + gamma(1,1) + beta(0.5, 0,5) * gamma(1,1)
  !!
  !! Gamma(1,1) distribution is just exponential distribution which is trivial to sample by
  !! inversion. Beta(0.5,0.5) distribution can be sampled using a method based on Johnk's theorem[2]
  !! which gives that for U and V being uniform random numbers on [0,1]:
  !!
  !! Beta(0.5, 0.5) = U^2 / (U^2 + V^2) provided that U^2 + V^2 <= 1.0
  !!
  !! To avoid rejection we can we can sample a uniform random angle phi on [0,PI/2] (quarter circle)
  !! and random radius r on [0,1] and substitute:
  !! U = r * cos(phi); V = r * sin(phi)
  !! Becouse new U and V fill the requirement of Johnk's theorem we obtain:
  !!
  !! Beta(0.5, 0.5) = cos(phi)^2
  !!
  !! So we can obtain sample of Beta(0.5, 0.5) by sampling random angle phi and taking square of
  !! its cosine (or sine for that matter, it doesn't matter)
  !!
  !!
  !! [1] Trkov, A., M. Herman, and D. A. Brown. “ENDF-6 Formats Manual.”
  !! Data Formats and Procedures for the Evaluated Nuclear Data Files ENDF/B-VI and ENDF/B-VII,
  !! National Nuclear Data Center Brookhaven National Laboratory, Upton, NY, 2012, 11973–5000.
  !!
  !! [2] Devroye, Luc. Non-Uniform Random Variate Generation. New York: Springer, 1986.
  !!
  type, public, extends(correlatedLawENDF) :: nBodyPhaseSpace
    private
    integer(shortInt) :: N   = 0   ! Number of secondary particles
    real(defReal)     :: Q  = ZERO ! Q-value for the reaction
    real(defReal)     :: Ap = ZERO ! Total mass ration for N-particles
    real(defReal)     :: A  = ZERO ! Target Atomic weight ratio
  contains
    ! Interface implementation
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    ! Public interface
    procedure :: init_fromACE
    procedure :: invalid


  end type nBodyPhaseSpace

contains

  !!
  !! Samples mu and E_out givent incident energy E_in and random nummber generator
  !! See description of the class for details about the sampling algorithm
  !!
  subroutine sample(self, mu, E_out, E_in, rand)
    class(nBodyPhaseSpace), intent(in) :: self
    real(defReal), intent(out)         :: mu
    real(defReal), intent(out)         :: E_out
    real(defReal), intent(in)          :: E_in
    class(RNG), intent(inout)          :: rand
    real(defReal)                      :: G1, G2, r1, r2, r3, r4, r5, r6, Emax
    type(maxwellEnergyPdf)             :: maxwellPdf
    character(100),parameter :: Here ='smaple (nBodyPhaseSpace_class.f90)'

    ! Sample mu
    mu = TWO * rand % get() - ONE

    ! Sample energy
    ! Get sample of G1 = Gamma(1.5,1) from maxwell distribution
    G1 = maxwellPdf % sample(ONE, rand)

    select case(self % N)
      case(3)
        ! Get sample of G2 = Gamma(1.5,1) from maxwell distribution
        G2 = maxwellPdf % sample(ONE, rand)

      case(4)
        ! Get sample of G2 = Gamma(3,1) from 3 exponential distributions
        r1 = rand % get()
        r2 = rand % get()
        r3 = rand % get()
        G2 = -log(r1 * r2 * r3)

      case(5)
        ! Get sample of G2 = Gamma(4.5,1) from 4 exponential distribution and 1 beta(0.5, 0.5)
        r1 = rand % get()
        r2 = rand % get()
        r3 = rand % get()
        r4 = rand % get()
        r5 = rand % get()
        r6 = rand % get()
        G2 = -log(r1*r2*r3*r4) - cos(TWO*PI*r5) * cos(HALF*PI*r5) * log(r6)

      case default ! Should never happen
        call fatalError(Here,'Wrong number of 2nd-ary particles on run-time')
        G2 = ZERO
    end select

    ! Calculate maximum energy
    Emax = (self % Ap - ONE) / self % Ap * (E_in * self % A / (self % A + ONE) + self % Q)

    E_out = Emax * G1 / (G1 + G2)

  end subroutine sample

  !!
  !! Returns probability that neutron was emmited at mu & E_out given incident energy E_in
  !!
  function probabilityOf(self, mu, E_out, E_in) result(prob)
    class(nBodyPhaseSpace), intent(in) :: self
    real(defReal), intent(in)          :: mu
    real(defReal), intent(in)          :: E_out
    real(defReal), intent(in)          :: E_in
    real(defReal)                      :: prob
    real(defReal)                      :: Emax
    real(defReal)                      :: C
    character(100),parameter :: Here ='probabilityOf (nBodyPhaseSpace_class.f90)'

    ! Calculate maximum energy
    Emax = (self % Ap - ONE) / self % Ap * (E_in * self % A / (self % A + ONE) + self % Q)

    select case(self % N)
      case(3)
        C = 4.0/(PI * Emax * Emax)
        prob = C * sqrt(E_out) * sqrt(Emax - E_out)

      case(4)
        C = 105.0/(32.0 * sqrt(Emax**7))
        prob = C * sqrt(E_out) * (Emax - E_out) * (Emax - E_out)

      case(5)
        C = 256.0/(14.0 * PI * Emax**5)
        prob = C * sqrt(E_out) * sqrt((Emax - E_out)**7)

      case default ! Should never happen
        call fatalError(Here,'Wrong number of 2nd-ary particles on run-time')
        prob = ZERO
    end select

    ! Multiply by uniform probability of mu
    prob = HALF * prob

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(nBodyPhaseSpace), intent(inout) :: self

    ! Set constants to initial values
    self % N  = 0
    self % Q  = ZERO
    self % Ap = ZERO
    self % A  = ZERO

  end subroutine kill


  !!
  !! Checks if the objeck is valid
  !! Returns .true. if it isn't. Optionaly sets what is wrong with an object
  !!
  function invalid(self, msg) result(isWrong)
    class(nBodyPhaseSpace), intent(in)              :: self
    character(:),allocatable, intent(out), optional :: msg
    logical(defBool)                                :: isWrong

    isWrong = .false.

    ! Verify value of N
    if( 5 < self % N .or. self % N < 3 ) then
      isWrong = .true.
      if(present(msg))  then
        msg = 'Invalid number of outgoin particles: '// numToChar(self % N) //' must be 3, 4 or 5'
      end if
      return
    end if

    ! Verify sign of atomic rations
    if (self % A < ZERO .or. self % Ap < ZERO) then
      isWrong = .true.
      if(present(msg))  then
        msg = 'A: ' // numToChar(self % A) // 'or Ap: '// numToChar(self % A) // ' is -ve'
      end if
      return
    end if

  end function invalid

  !!
  !! Subroutine initialise from ACE, Q-value and atomic weight ratio of traget
  !! aceCard read head needs to be set to the beginning of the data
  !!
  subroutine init_fromACE(self, ACE, Q, A)
    class(nBodyPhaseSpace), intent(inout) :: self
    class(aceCard), intent(inout)         :: ACE
    real(defReal), intent(in)             :: Q
    real(defReal), intent(in)             :: A
    character(:),allocatable              :: msg
    character(100), parameter :: Here ='init_fromACE (nBodyPhaseSpace_class.f90)'

    ! Read and store data
    self % N  = ACE % readInt()
    self % Ap = ACE % readReal()
    self % Q  = Q
    self % A  = A

    ! Varify correctness
    if(self % invalid(msg)) call fatalError(Here, msg)

  end subroutine init_fromACE

  !!
  !! Constructor from ACE, Q and A
  !! aceCard read head needs to be set to the beginning of the data
  !!
  function new_nBodyPhaseSpace_fromACE(ACE, Q, A) result(new)
    class(aceCard), intent(inout) :: ACE
    real(defReal), intent(in)     :: Q
    real(defReal), intent(in)     :: A
    type(nBodyPhaseSpace)         :: new

    call new % init_fromACE(ACE, Q, A)

  end function new_nBodyPhaseSpace_fromACE

  !!
  !! Constructor from values
  !!
  function new_nBodyPhaseSpace_fromValues(N, Ap, Q, A) result(new)
    integer(shortInt), intent(in) :: N
    real(defReal), intent(in)     :: Ap
    real(defReal), intent(in)     :: Q
    real(defReal), intent(in)     :: A
    type(nBodyPhaseSpace)         :: new
    character(:),allocatable      :: msg
    character(100), parameter :: Here ='new_nBodyPhaseSpace_fromValues (nBodyPhaseSpace_class.f90)'

    new % N  = N
    new % Ap = Ap
    new % Q  = Q
    new % A  = A

    ! Varify correctness
    if(new % invalid(msg)) call fatalError(Here, msg)

  end function new_nBodyPhaseSpace_fromValues


end module nBodyPhaseSpace_class
