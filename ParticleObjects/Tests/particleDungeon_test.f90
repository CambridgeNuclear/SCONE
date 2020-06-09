module particleDungeon_test
  use numPrecision
  use RNG_class,             only : RNG
  use particle_class,        only : particle, particleState
  use particleDungeon_class, only : particleDungeon
  use pFUnit_mod

  implicit none


contains

  !!
  !! Test stack like access. Test is a dummy use case
  !!
@Test
  subroutine testStackInterface()
    type(particle)        :: p
    type(particleState)   :: phase
    type(particleDungeon) :: dungeon
    integer(shortInt)     :: i

    ! Is empty result for uninitialised dungeon
    @assertTrue(dungeon % isEmpty())

    ! Initialise
    call dungeon % init(10)

    ! Is empty initialised. No particles stored
    @assertTrue(dungeon % isEmpty())
    @assertEqual(0, dungeon % popSize())

    ! Store some particles with energy
    p % isMG = .false.
    p % E  = 8.6_defReal

    do i = 1,5
      call dungeon % detain(p)
    end do

    ! Store some phase coordinates with energy
    phase % isMG = .false.
    phase % E = 3.4_defReal

    do i = 1,4
      call dungeon % detain(phase)
    end do

    ! Verify size
    @assertFalse(dungeon % isEmpty())
    @assertEqual(9, dungeon % popSize())

    ! Remove particles
    do i = 1,4
      call dungeon % release(p)
      @assertEqual(3.4_defReal, p % E)
    end do

     do i = 1,5
      call dungeon % release(p)
      @assertEqual(8.6_defReal, p % E)
    end do

    ! Verify size
    @assertTrue(dungeon % isEmpty())
    @assertEqual(0, dungeon % popSize())

    ! Clean
    call dungeon % kill()
  end subroutine testStackInterface

  !!
  !! Test array like interface
  !!
@Test
  subroutine testArrayInterface
    type(particle)        :: p
    type(particleState)   :: phase
    type(particleDungeon) :: dungeon
    integer(shortInt)     :: i

    ! Is empty result for uninitialised dungeon
    @assertTrue(dungeon % isEmpty())

    ! Initialise to fixed size
    call dungeon % setSize(2)

    ! Is filled with random particles
    @assertFalse(dungeon % isEmpty())
    @assertEqual(2, dungeon % popSize())

    ! Extend size
    call dungeon % setSize(9)
    @assertFalse(dungeon % isEmpty())
    @assertEqual(9, dungeon % popSize())

    ! Store some particles with energy
    p % isMG = .false.
    p % E  = 8.6_defReal

    do i = 1,5
      call dungeon % replace(p, i)
    end do

    ! Store some phase coordinates with energy
    phase % isMG = .false.
    phase % E = 3.4_defReal

    do i = 1,4
      call dungeon % replace(phase, 5 + i)
    end do

    ! Verify size
    @assertFalse(dungeon % isEmpty())
    @assertEqual(9, dungeon % popSize())

    ! Raplace by particle and phaseCoords
    p % E = 13.0_defReal
    phase % E = 17.0_defReal
    call dungeon % replace(p, 9)
    call dungeon % replace(phase, 1)

    ! Verify particles by copies
    call dungeon % copy(p, 9)
    @assertEqual(13.0_defReal, p % E)

    do i=8,6,-1
      call dungeon % copy(p, i)
      @assertEqual(3.4_defReal, p % E)
    end do

    do i=5,2,-1
      call dungeon % copy(p, i)
      @assertEqual(8.6_defReal, p % E)
    end do
    call dungeon % copy(p, 1)
    @assertEqual(17.0_defReal, p % E)

    ! Verify that population has not changed
    @assertFalse(dungeon % isEmpty())
    @assertEqual(9, dungeon % popSize())

    ! Shrink size
    call dungeon % setSize(2)
    @assertFalse(dungeon % isEmpty())
    @assertEqual(2, dungeon % popSize())

    ! Clean population
    call dungeon % cleanPop()
    @assertTrue(dungeon % isEmpty())
    @assertEqual(0, dungeon % popSize())

    ! Clean memory
    call dungeon % kill()

  end subroutine testArrayInterface

  !!
  !! Test weight normalisation and inquiry
  !!
@Test
  subroutine testWeightNorm()
    type(particle)        :: p
    type(particleDungeon) :: dungeon
    integer(shortInt)     :: i
    real(defReal), parameter :: TOL = 1.0E-9

    ! Initialise
    call dungeon % init(10)
    @assertEqual(ZERO, dungeon % popWeight(), TOL)

    ! Store some particles with non-uniform weight
    do i = 1,5
      p % w = 0.5_defReal + i * 0.1_defReal
      call dungeon % detain(p)

    end do

    ! Verify total weight
    @assertEqual(4.0_defReal, dungeon % popWeight(), TOL)

    ! Normalise weight
    call dungeon % normWeight(2.0_defReal)
    @assertEqual(2.0_defReal, dungeon % popWeight(), TOL)

    ! Get particles and compare weight
    do i = 5,1,-1
      call dungeon % release(p)
      @assertEqual( 0.25_defReal + i * 0.05_defReal, p % w, TOL)
    end do

    ! Verify weight
    @assertEqual(0.0_defReal, dungeon % popWeight(), TOL)

    ! Clean
    call dungeon % kill()
  end subroutine testWeightNorm

  !!
  !! Test normalisation of population to smaller number
  !! Particles with non-uniform weight
  !!  NOTE: Weight preservation is disabled for now
  !!
@Test
  subroutine testNormPopDown()
    type(particleDungeon)    :: dungeon
    type(particle)           :: p
    type(RNG)                :: pRNG
    integer(shortInt)        :: i
    real(defReal), parameter :: TOL = 1.0E-9


    ! Initialise
    call dungeon % init(10)
    call pRNG % init(7865856_longInt)

    ! Store some particles with non-uniform weight
    do i = 1,10
      p % w = 0.5_defReal + i * 0.1_defReal
      call dungeon % detain(p)
    end do

    ! Normalise population
    call dungeon % normSize(5, pRNG)

    ! Verify size
    @assertEqual(5, dungeon % popSize())

    ! Verify weight *** DISABLED
    !@assertEqual(6.05_defReal, dungeon % popWeight(), TOL)

    ! Clean memory
    call dungeon % kill()

  end subroutine testNormPopDown

  !!
  !! Test normalisation of population to smaller number
  !! Particles with non-uniform weight
  !!  NOTE: Weight preservation is disabled for now
  !!
@Test
  subroutine testNormPopUp()
    type(particleDungeon)    :: dungeon
    type(particle)           :: p
    type(RNG)                :: pRNG
    integer(shortInt)        :: i
    real(defReal), parameter :: TOL = 1.0E-9

    ! Initialise
    call dungeon % init(20)
    call pRNG % init(435468_longInt)

    ! Store some particles with non-uniform weight
    do i = 1,10
      p % w = 0.5_defReal + i * 0.1_defReal
      call dungeon % detain(p)
    end do

    ! Normalise population
    call dungeon % normSize(15, pRNG)

    ! Verify size
    @assertEqual(15, dungeon % popSize())

    ! Verify weight *** DISABLED
    !@assertEqual(18.15_defReal, dungeon % popWeight(), TOL)

    ! Clean memory
    call dungeon % kill()

  end subroutine testNormPopUp


end module particleDungeon_test
