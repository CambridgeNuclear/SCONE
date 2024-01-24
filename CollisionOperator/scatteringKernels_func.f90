module scatteringKernels_func

  use numPrecision
  use genericProcedures, only : rotateVector
  use RNG_class,         only : RNG
  use particle_class,    only : particle

  ! Nuclear Data
  use ceNeutronNuclide_inter, only : ceNeutronNuclide

  implicit none
  private

  public  :: asymptoticScatter
  public  :: targetVelocity_constXS
  public  :: targetVelocity_DBRCXS
  public  :: asymptoticInelasticScatter

  private :: sample_x2expx2
  private :: sample_x3expx2
  private :: sample_targetVelocity

contains


  !!
  !! Subroutine to perform ELASTIC SCATTERING from a stationary target nucleus.
  !! Changes direction and energy in LAB system
  !! Takes mu in CM system, returns it in LAB system
  !!
  !! Dir needs to be normalised. Prodedure will produce wrong results without error message
  !! if it is not.
  !!
  !! Based on MCNP manual chapter 2.
  !!
  subroutine asymptoticScatter(E,mu,A)
    real(defReal), intent(inout)               :: E    !! Pre-collision energy in LAB
    real(defReal), intent(inout)               :: mu   !! Cosine of delection angle
    real(defReal), intent(in)                  :: A    !! Target mass [neutrons]
    real(defReal)                              :: E_in
    real(defReal)                              :: inv_Ap1

    ! Store initial energy and precalculate 1/(A+1)
    E_in = E
    inv_Ap1 = 1.0/ (A + 1.0)

    ! Find post-collision energy
    E  = (1.0 + A*A + 2 *A*mu) *E_in * inv_Ap1 * inv_Ap1

    ! Find deflection angle in LAB
    mu = (A*mu + 1)*sqrt(E_in/E)* inv_Ap1

    ! Correct possible nuclear data shortcomings
    if (mu > ONE) mu = ONE

  end subroutine asymptoticScatter

  !!
  !! Subroutine to perform INELASTIC SCATTERING from a stationary target nucleus.
  !! Changes direction and energy in LAB system
  !! Takes mu in CM system, and E_out in CM system.
  !! Returns mu in LAB system.
  !!
  !! Dir needs to be normalised. Prodedure will produce wrong results without error message
  !! if it is not.
  !!
  !! Based on MCNP manual chapter 2.
  !!
  subroutine asymptoticInelasticScatter(E,mu,E_out,A)
    real(defReal), intent(inout) :: E     !! Pre-collision energy in Lab
    real(defReal), intent(inout) :: mu    !! Cosine of delection angle
    real(defReal), intent(in)    :: E_out !! Post collision energy in CM
    real(defReal), intent(in)    :: A     !! Target mass [neutrons]
    real(defReal)                :: E_in
    real(defReal)                :: inv_Ap1

    ! Store initial energy and precalculate 1/(A+1)
    E_in = E
    inv_Ap1 = 1.0/ (A + 1.0)

    ! Find post-collision energy
    E = E_out + (E_in +TWO*mu*(A+ONE)*sqrt(E_in*E_out)) * inv_Ap1 * inv_Ap1

    ! Find deflection angle in LAB
    mu = mu * sqrt(E_out/E) + sqrt(E_in/E)* inv_Ap1

    ! Correct possible nuclear data shortcomings
    if (mu > ONE) mu = ONE

  end subroutine asymptoticInelasticScatter


  !!
  !! Function that returns a sample of target velocity using constant XS approximation
  !! V_t is a vector. The velocity is scaled by a factor sqrt(Mn/2) where Mn is mass of a neutron
  !! so that V_t*V_t=E_t with E_t beeing kinetic energy of a NEUTRON traveling with TARGET VELOCITY
  !! (note that it is not a kinetic energy of the target).
  !!
  function targetVelocity_constXS(E,dir,A,kT,rand) result (V_t)
    real(defReal), intent(in)               :: E
    real(defReal), dimension(3), intent(in) :: dir
    real(defReal), intent(in)               :: A
    real(defReal), intent(in)               :: kT
    class(RNG), intent(inout)               :: rand
    logical(defBool)                        :: accept
    real(defReal),dimension(3)              :: V_t
    real(defReal)                           :: alpha, mu, phi
    real(defReal)                           :: X, Y
    real(defReal)                           :: rel_v, r1

    ! Calculate neutron Y = beta *V_n
    ! beta = sqrt(A*Mn/2kT). Note velocity scaling by sqrt(Mn/2).
    Y = sqrt(A*E/kT)

    ! Calculate threshold factor alpha
    alpha = 2.0/(Y*sqrt(PI)+2.0)

    rejectionLoop: do

      ! Sample velocity and calculate angle and acceptance probability
      call sample_targetVelocity(X, accept, rel_v, mu, rand, Y, alpha)

      ! Accept or reject mu
      if (accept) exit rejectionLoop

    end do rejectionLoop

    ! Calculate azimithal angle for target and obtain target direction
    r1 = rand % get()
    phi = 2.0 *PI * r1

    V_t = rotateVector(dir, mu, phi)

    ! Scale target direction by magnitude of velocity
    V_t = V_t * (X * sqrt(kT/A))

  end function targetVelocity_constXS

  !!
  !! Function that returns a sample of target velocity with DBRC
  !! V_t is a vector. The velocity is scaled by a factor sqrt(Mn/2) where Mn is mass of a neutron
  !! so that V_t*V_t=E_t with E_t being kinetic energy of a NEUTRON traveling with TARGET VELOCITY
  !! (note that it is not a kinetic energy of the target).
  !!
  !!
  function targetVelocity_DBRCXS(nuc, E, dir, A, kT, rand, tempMaj) result (V_t)
    class(ceNeutronNuclide), intent(in)     :: nuc
    real(defReal), intent(in)               :: E
    real(defReal), dimension(3), intent(in) :: dir
    real(defReal), intent(in)               :: A
    real(defReal), intent(in)               :: kT
    real(defReal), intent(in)               :: tempMaj
    class(RNG), intent(inout)               :: rand
    logical(defBool)                        :: accept
    real(defReal),dimension(3)              :: V_t
    real(defReal)                           :: alpha, mu, phi, DBRC_acc
    real(defReal)                           :: X, Y
    real(defReal)                           :: r1, r2
    real(defReal)                           :: rel_v, rel_E, xs_rel_v

    ! Calculate neutron Y = beta *V_n
    ! beta = sqrt(A*Mn/2kT). Note velocity scaling by sqrt(Mn/2).
    Y = sqrt(A * E / kT)

    ! Calculate threshold factor alpha
    ! In MCNP, alpha is p1
    alpha = 2.0 / (Y * sqrt(PI) + 2.0)

    rejectionLoop: do

      ! Sample velocity and calculate angle and acceptance probability
      call sample_targetVelocity(X, accept, rel_v, mu, rand, Y, alpha)

      ! Accept or reject mu and target velocity
      if (.not. accept) cycle rejectionLoop

      ! Relative energy = relative velocity **2 due to sqrt(Mn/2) scaling factor
      rel_E = (rel_v**2 * kT / A)

      ! Find 0K scattering xs of target at relative energy
      xs_rel_v = nuc % elScatteringXS(rel_E)

      ! Introduce DBRC acceptance condition
      DBRC_acc = (xs_rel_v / tempMaj)
      r1 = rand % get()

      ! Accept or reject with DBRC
      if (DBRC_acc > r1) then
        exit rejectionLoop
      end if

    end do rejectionLoop

    ! Calculate azimithal angle for target and obtain target direction
    r2 = rand % get()
    phi = 2.0 * PI * r2

    V_t = rotateVector(dir, mu, phi)

    ! Scale target direction by magnitude of velocity
    V_t = V_t * (X * sqrt(kT / A))

  end function targetVelocity_DBRCXS


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

  !!
  !! Subroutine to sample the target velocity using the predefined helper functions,
  !! sample the scattering angle, and calculate the acceptance probability for the
  !! sampling loop.
  !!
  !! It is called used by different methods (contant xs, DBRC) when sampling the velocity
  !!
  subroutine sample_targetVelocity(X, accept, rel_v, mu, rand, Y, alpha)
    real(defReal), intent(out)    :: X
    logical(defBool), intent(out) :: accept
    real(defReal), intent(out)    :: rel_v
    real(defReal), intent(out)    :: mu
    class(RNG), intent(inout)     :: rand
    real(defReal), intent(in)     :: Y
    real(defReal), intent(in)     :: alpha
    real(defReal)                 :: r1, r2, r3, P_acc

    ! Obtain random numbers
    r1 = rand % get()
    r2 = rand % get()
    r3 = rand % get()

    ! Sample X = beta * V_t
    ! Uses helper functions below
    if ( r1 > alpha ) then
      X = sample_x2expx2(rand)
    else
      X = sample_x3expx2(rand)
    end if

    ! Sample polar angle of target velocity wrt. neutron direction
    mu = 2.0 * r2 - 1.0;

    ! Calculate relative velocity between neutron and target
    rel_v = sqrt(Y * Y + X * X - 2.0 * X * Y * mu)

    ! Calculate Acceptance Propability
    P_acc = rel_v / (Y + X)

    ! Verify acceptance condition
    accept = P_acc > r3

  end subroutine sample_targetVelocity


end module scatteringKernels_func
