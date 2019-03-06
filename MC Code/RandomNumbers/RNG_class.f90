module rng_class
  use numPrecision
  use iso_fortran_env, only : int64
  implicit none
  private

  !!
  !! Linear congruential 63 bit random number generator
  !!
  !! Follows recurrance formula: xi(i+1) = (g * xi(i) + c) mod M
  !!
  !! Global Parameters (values based on OpenMC):
  !!    g: multiplier = 2806196910506780709
  !!    c: increment = 1
  !!    M: modulus = 2**63
  !!
  !! NOTE: M chosen to be a power of 2 -> simplifies integer division (right shift instruction)
  !!       Definition of bitMask assumes that sign-bit is leftmost bit. Need test for that.
  !!
  type, public :: rng
    private
    integer(int64) :: rngSeed
    integer(int64) :: rngCount
    integer(int64) :: initialSeed
  contains
    procedure :: init
    procedure :: get
    procedure :: getInt
    procedure :: skip
    procedure :: getCount
    procedure :: getSeed
  end type rng

  !! Parameters
  !! NOTE: Bit enumeration in ibset begins at 0. (Section 13.3.1 of Fortran 2008 Standard)
  !!       63-th position is leftmost bit!
  integer(int64), parameter :: g = 2806196910506780709_int64 !2806191050678079_int64
  integer(int64), parameter :: c = 1_int64
  integer(int64), parameter :: M = ibset(0_int64,63)
  integer(int64), parameter :: bitMask = huge(0_int64) ! All bits but the 63th are 1
  real(defReal), parameter  :: norm = ONE / 2.0_defReal**(63)


contains

  !!
  !! Initialise RNG with a seed
  !!
  subroutine init(self, seed)
    class(rng), intent(inout)  :: self
    integer(int64), intent(in) :: seed

    self % rngSeed     = seed
    self % initialSeed = seed

  end subroutine init

  !!
  !! Returns value of random number on <0,1)
  !!
  function get(self) result(rand)
    class(rng), intent(inout) :: self
    real(defReal)             :: rand
    integer(int64)            :: seed

    ! Get current state of LCG
    seed = self % rngSeed

    ! Multiply by multiplier and keep rightmost 63 bits
    seed = iand(g * seed, bitMask)

    ! Add increment and keep rightmost 63 bits
    seed = iand(seed + c, bitMask)

    ! Convert integer LCG state to real number on <0,1)
    rand = seed * norm

    ! Update RNG state
    self % rngSeed  = seed
    self % rngCount = self % rngCount + 1

  end function get

  !!
  !! Return random integer instead of real
  !!
  function getInt(self) result(seed)
    class(rng), intent(inout) :: self
    integer(int64)            :: seed

    ! Get current state of LCG
    seed = self % rngSeed

    ! Multiply by multiplier and keep rightmost 63 bits
    seed = iand(g * seed, bitMask)

    ! Add increment and keep rightmost 63 bits
    seed = iand(seed + c, bitMask)

    ! Update RNG state
    self % rngSeed  = seed
    self % rngCount = self % rngCount + 1

  end function getInt

  !!
  !! Move state of the LCG by k forward
  !! -ve k moves state backwards
  !!
  !! Uses periodicity of LCG. Depends on the properties of modular arithmetic.
  !! [https://brilliant.org/wiki/modular-arithmetic/]
  !!
  !! The algorithim was developed by Forrest B. Brown:
  !! [https://laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/anl-rn-arb-stride.pdf]
  !! Justification for why the algorithim for evaluation of Ck works is made by me. [MAK]
  !!
  !! Assume that we are given a LCG in a state S0 and we are interested to find its state S_k after
  !! k steps. Then starting with recurrence relation we can expand the recurrance k times:
  !!
  !! S_k = g * S_(k-1) + c (mod M)
  !! S_k = g * ([ g * S_(k-2) + c (mod M)]) + c (mod M) = g**2 * S_(k-2) + c * ( g + 1) (mod M)
  !!
  !! Eventually using expression for a sum of geometric series:
  !!
  !! S_k = g**k * S0 + c * (g**k - 1) /(g - 1) (mod M) = Gk * S0 + Ck (mod M)
  !! S_k = Gk * S0 (mod M) + Ck (mod M)
  !!
  !! Thus the goal is now to estimate both of the terms in the sum by some smart algorithm.
  !! To calculate Gk consider binary representation of k and substitute it into exponent of g.
  !! Then it is possible to show that:
  !!
  !! Gk = g**k = (g**1)**k_0 * (g**2)**k_1 * (g**4)**k_2 * ... (mod M)
  !!
  !! Where k_i denotes the i-th bit of k binary representation and can be 0 or 1. Now the evaluation
  !! of Gk is trivial by noting that (g**n)**k = 1 if k=0 or (g**n)**k = g**n if k = 1.
  !!
  !! The evaluation of Ck is slightly more compilcated and it is based on two recurrance relations
  !! for the sum of geometric series. Denote L(k) to be a geometric series such that
  !! L(k) = 1 + g + g**2 + g**3 + ... + g**(k-1). And define L(0) = 0.
  !! Then following relations hold:
  !!
  !!   1) L(k)    = L(k-n) * g**n + L(n) for k >= n
  !!   2) L(k**2) = (g**k +1) * L(k)
  !!
  !! Relation 1) can be proven by noting that multiplication by g**n shifts geometric series
  !! n-times to the left and C(n) adds the missing elements on the right.
  !!
  !! Relation 2) can be proven by expanding the sum and using (x**2-1) = (x+1)(x-1)
  !!
  !! Using the binary expansion of k.  Denote highest bit index of k to be m.
  !! Than using the recurrance relation 1) it can be shown that:
  !!
  !! C(k) = C(k_m * 2**m + R_m)= C(R_m) * g**(2**m * k_m) + c * L(2**m * k_m)
  !!
  !! Where k_i is in {0;1}; k = k_m * 2**m + R_m and R_m = k_(m-1) * 2**(m-1) + R_(m-1) etc.
  !!
  !! With a relation 2) to evaluate L(2**m * k_m) this allows to evaluate C(k) in log(k) steps.
  !!   Note that if k_i == 0 then C(R_i) = C(R_(i-1))
  !!
  !! NOTE: variable names were changed from original Forrest Brown algorithim to be more
  !!       descreptive:
  !!       h -> gSq_to_i (g**(2**i))
  !!       f -> L (L as defined above)
  !!
  subroutine skip(self, k)
    class(rng)     :: self
    integer(int64) :: k         ! number of places to skip
    integer(int64) :: Gk        ! G**k (mod M)
    integer(int64) :: Ck        ! c*(g**k-1)/(g-1) (mod M)
    integer(int64) :: gSq_to_i  ! g_squared to power of i (g**(2**i))
    integer(int64) :: L         ! Sum of geometric series as defined above

    ! Set initial values
    Gk       = 1
    Ck       = 0
    gSq_to_i = g
    L        = c

    ! Can translate jump backwards to jump forwards due to periodicity of RNG
    ! For our settings period is M
    if(k < 0) k = k + M

    ! Unnecessary line. Sign bit of k is already 0
    k = iand(k + M, bitMask)

    do while( k > 0)
      if(iand(k, 1_int64)== 1) then ! Right-most bit is 1
        Gk = iand(Gk * gSq_to_i, bitMask)     ! Add to Gk
        Ck = iand(Ck * gSq_to_i + L, bitMask) ! Add to Ck

      end if

      L = iand(L * (gSq_to_i+1), bitMask)           ! Calculate next value of L
      gSq_to_i = iand(gSq_to_i * gSq_to_i, bitMask) ! Calculate next power of g**2
      k = ishft(k, -1)                              ! Right shift k by 1
    end do

    ! Jump forward
    self % rngSeed = iand(Gk * self % rngSeed + Ck, bitMask)

  end subroutine skip

  !!
  !! Return total number of psudo-random numbers generated
  !!
  function getCount(self) result (c)
    class(rng),intent(in) :: self
    integer(int64)        :: c

    c = self % rngCount

  end function getCount

  !!
  !! Returns value of seed used to initialise RNG
  !!
  function getSeed(self) result(seed)
    class(rng), intent(in) :: self
    integer(int64)         :: seed

    seed = self % initialSeed

  end function getSeed

end module rng_class
