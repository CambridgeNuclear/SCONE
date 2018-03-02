module scatteringKernels_func

  use numPrecision
  use genericProcedures, only : rotateVector
  use RNG_class,         only : RNG

  implicit none
  private

  public  :: asymptoticScatter
  public  :: targetVelocity_constXS


  private :: sample_x2expx2
  private :: sample_x3expx2

contains


  !!
  !! Subroutine to perform scattering from a stationary target nucleus.
  !! Changes direction and energy in LAB system
  !! Takes mu in CM system, returns it in LAB system
  !!
  !! Dir needs to be normalised. Prodedure will produce wrong results without error message
  !! if it is not.
  !!
  !! Based on MCNP manual chapter 2.
  !!
  subroutine asymptoticScatter(dir,E,mu,A)
    real(defReal), dimension(3), intent(inout) :: dir  !! Normalised direction vector
    real(defReal), intent(inout)               :: E    !! Pre-collision energy in LAB
    real(defReal), intent(inout)               :: mu   !! Cosine of delection angle
    real(defReal), intent(in)                  :: A    !! Target mass [neutrons]
    real(defReal)                              :: E_in
    real(defReal)                              :: inv_Ap1

    ! Store initial energy and precalculate 1/(A+1)
    E_in = E
    inv_Ap1 = 1.0/ (A + 1.0)

    ! Find post-collision energy
    E  = (1.0 + A*A + 2 *mu) *E_in * inv_Ap1 * inv_Ap1

    ! Find deflection angle in LAB
    mu = mu + sqrt(E_in/E)* inv_Ap1

  end subroutine asymptoticScatter



  !!
  !! Function that returns a sample of target velocity using constant XS approximation
  !! V_t is a vector. The velocity is scaled by a factor sqrt(Mn/2) where Mn is mass of a neutron
  !! so that V_t*V_t=E_t with E_t beeing kinetic energy of a NEUTRON traveling with TARGET VELOCITY
  !! (note that it is not a kinetic energy of the target).
  !!
  !!
  function targetVelocity_constXS(E,dir,kT,A,rand) result (V_t)
    real(defReal), intent(in)               :: E
    real(defReal), dimension(3), intent(in) :: dir
    real(defReal), intent(in)               :: kT
    real(defReal), intent(in)               :: A
    class(RNG), intent(inout)               :: rand
    real(defReal),dimension(3)              :: V_t
    real(defReal)                           :: V_n
    real(defReal)                           :: alpha, mu, phi, P_acc
    real(defReal)                           :: X, Y
    real(defReal)                           :: r1, r2, r3, r4

    ! Calculate neutron Y = beta *V_n
    ! beta = sqrt(A*Mn/2kT). Note velocity scaling by sqrt(Mn/2).
    Y = sqrt(A*E/kT)

    ! Calculate treshhold factor alpha
    alpha = 2.0/(Y*sqrt(PI)+2.0)


    rejectionLoop: do

      ! Obtain random numbers
      r1 = rand % get()
      r2 = rand % get()
      r3 = rand % get()

      ! Sample X = beta * V_t
      if ( r1 > alpha ) then
        X = sample_x2expx2(rand)

      else
        X = sample_x3expx2(rand)

      end if

      ! Sample polar angle of target velocity wrt. neutron direction
      mu = 2.0 * r2 - 1.0;

      ! Calculate Acceptance Propability
      P_acc = sqrt(Y*Y + X*X - 2.0*X*Y*mu) / (Y+X)

      ! Accept or reject mu
      if (P_acc > r3) exit rejectionLoop

    end do rejectionLoop

    ! Calculate azimithal angle for traget and obtain target direction
    r4 = rand % get()
    phi = 2.0 *PI * r4

    V_t = dir
    call rotateVector(V_t, mu, phi)

    ! Scale target direction by magnitude of velocity
    V_t = V_t * (X * sqrt(kT/A))

  end function targetVelocity_constXS



  !!
  !! Helper function to sample x^2 * exp( - x^2) probability distribution
  !! Uses random numbers from provided RNG
  !! Changes variables to y = x^2 which transforms PDF to Gamma(3/2,1) distribution
  !! Then uses method based on Johnk's theorem and sum of Gamma distributed random variables
  !!
  function sample_x2expx2(rand) result(sample)
    class(RNG), intent(inout) :: rand
    real(defReal)             :: sample
    real(defReal)             :: r1, r2, r3
    real(defReal)             :: beta, gamma05, cosine

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
    sample = -log(r3) + gamma05

    ! Change variables back to x from y=x^2
    sample = sqrt(sample)

  end function sample_x2expx2


  !!
  !! Helper function to sample x^3 * exp( - x^2) probability distribution
  !! Uses random numbers from provided RNG
  !! Changes variables to y = x^2 which transforms PDF to Gamma(2,1) distribution
  !! Sampling Gamma(2,1) is trivial using sum of Gamma distributed random variables
  !!
  function sample_x3expx2(rand) result(sample)
    class(RNG), intent(inout) :: rand
    real(defReal)             :: sample
    real(defReal)             :: r1, r2

    ! Obtain random numbers
    r1 = rand % get()
    r2 = rand % get()

    ! Sample Gamma(2,1) by summing two samples of Gamma(1,1) [exponential distribution]
    sample = -log(r1) - log(r2)

    ! Change variables back to x
    sample = sqrt(sample)

  end function sample_x3expx2


end module scatteringKernels_func
