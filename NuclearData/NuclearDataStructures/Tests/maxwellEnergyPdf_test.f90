module maxwellEnergyPdf_test

  use funit
  use iso_fortran_env,        only : int64
  use maxwellEnergyPdf_class, only : maxwellEnergyPdf
  use numPrecision
  use rng_class,              only : rng

  implicit none

contains

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! The generator's state space includes 0.0, so get() can return exactly 0.0 once
  !! per period. The seed below is chosen so that the SECOND draw of the sampler
  !! is exactly 0.0: with x(i+1) = (g*x(i) + 1) mod 2**63, g = 2806196910506780709. 
  !! An unguarded -log(r2) then yields +Inf energy; the sample must stay finite and 
  !! non-negative.
  !!
@Test
  subroutine testSampleAtZero()
    real(defReal)             :: E
    type(maxwellEnergyPdf)    :: mPdf
    type(rng)                 :: rand
    integer(int64), parameter :: SEED = 6426352591675910506_int64

    ! Initialise RNG.
    call rand % init(SEED)

    ! Sample energy and check that it is finite and non-negative.
    E = mPdf % sample(1.0_defReal, rand)
    @assertTrue(E == E)
    @assertTrue(E >= ZERO)
    @assertTrue(E < huge(E))

  end subroutine testSampleAtZero

end module maxwellEnergyPdf_test