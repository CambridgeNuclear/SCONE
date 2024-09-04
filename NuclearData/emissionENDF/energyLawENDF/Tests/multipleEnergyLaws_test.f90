module multipleEnergyLaws_test

  use numPrecision
  use RNG_class,                only : RNG
  use energyLawENDF_inter,      only : energyLawENDF
  use multipleEnergyLaws_class, only : multipleEnergyLaws
  use testEnergyLaw_class,      only : testEnergyLaw
  use funit

  implicit none

  ! Trivial fractional probability data from Pu-240 MT=37 JEFF 3.1.1
  real(defReal), dimension(*), parameter :: eGrid = [18.11606_defReal, 20.0_defReal]
  real(defReal), dimension(*), parameter :: pdf   = [0.25_defReal, 0.25_defReal]

@testCase
  type, extends(TestCase) :: test_multipleEnergyLaws
    private
    type(multipleEnergyLaws) :: law
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_multipleEnergyLaws


contains

  !!
  !! Sets up test_multipleEnergyLaws object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_multipleEnergyLaws), intent(inout) :: this
    type(testEnergyLaw)                           :: testEL
    class(energyLawENDF), allocatable             :: energyLaw
    integer(shortInt)                             :: i
    real(defReal),dimension(4), parameter :: E_Outs = [ 13.0_defReal, 18.0_defReal, 9.0_defReal, 13.0_defReal]

    ! Allocate space
    call this % law % init(6)

    ! Load energy laws
    do i=1,4
      testEL % E_out = E_outs(i)
      allocate(energyLaw, source = testEL)
      call this % law % addLaw(energyLaw, eGrid, pdf)
    end do

  end subroutine setUp

  !!
  !! Kills test_multipleEnergyLaws object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_multipleEnergyLaws), intent(inout) :: this

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test getting a probability
  !!
@Test
  subroutine testProbabilityOf(this)
    class(test_multipleEnergyLaws), intent(inout) :: this
    real(defReal), parameter :: TOL = 1.0E-7

    ! Immposible energy
    @assertEqual(ZERO, this % law % probabilityOf(1.0_defReal, 17.0_defReal), TOL)

    ! Valid energy with low capping
    @assertEqual(0.5_defReal, this % law % probabilityOf(13.0_defReal, 13.0_defReal), TOL)

    ! Valid energy with top capping
    @assertEqual(0.25_defReal, this % law % probabilityOf(18.0_defReal, 33.0_defReal), TOL)

    ! Valid energy without capping
    @assertEqual(0.25_defReal, this % law % probabilityOf(9.0_defReal, 19.0_defReal), TOL)

  end subroutine testProbabilityOf


  !!
  !! Sampling tests
  !!
@Test
  subroutine testSampling(this)
    class(test_multipleEnergyLaws), intent(inout) :: this
    real(defReal),dimension(1200)                 :: samples
    type(RNG)                                     :: rand
    real(defReal)                                 :: E_in
    integer(shortInt)                             :: i
    real(defReal)                                 :: C1, C2, C3, Chi

    ! Initialise RNG
    call rand % init(17_longInt)

    ! Draw samples with low capping
    E_in = 5.0_defReal
    do i=1,400
      samples(i) = this % law % sample(E_in, rand)
    end do

    ! Draw with top capping
    E_in = 30.1_defReal
    do i=401,800
      samples(i) = this % law % sample(E_in, rand)
    end do

    ! Draw with no capping
    E_in = 19.1_defReal
    do i=801,1200
      samples(i) = this % law % sample(E_in, rand)
    end do

    ! Preform some tests
    ! We can count number of each type of E_out.
    ! As the result we can apply Pearson's Chi-2 test
    !
    C1 = real(count(samples == 13.0_defReal), defReal)
    C2 = real(count(samples == 9.0_defReal), defReal)
    C3 = real(count(samples == 18.0_defReal), defReal)

    Chi = (C1-600.0_defReal)**2/600_defReal + &
          (C2-300.0_defReal)**2/300_defReal + &
          (C3-300.0_defReal)**2/300_defReal

    ! Value of 7.38 was chosen from 2-degree of freedom Chi-Sq distribution
    ! To obtain 97.5% confidence interval
    @assertGreaterThan(7.38_defReal, Chi)

  end subroutine testSampling

end module multipleEnergyLaws_test
