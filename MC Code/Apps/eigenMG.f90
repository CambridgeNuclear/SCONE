program eigenMG

  use numPrecision
  use RNG_class,                 only : RNG
  use byNucNoMT_class,           only : byNucNoMT
  use particle_class,            only : particle
  use particleDungeon_class,     only : particleDungeon

  use isotropicMG_class,         only : isotropicMG
  !use collisionOperatorMG_class, only : collisionOperatorMG
  use perMaterialCollisionOpMG_class, only : perMaterialCollisionOpMG

  use dictionary_class ,         only : dictionary
  use IOdictionary_class,        only : IOdictionary

  implicit none
  ! ### Declaration Begin Here ### !

  type(particle)            :: neutron
  !type(collisionOperatorMG) :: collisionPhysics
  type(perMaterialCollisionOpMG) :: collisionPhysics

  class(RNG), pointer       :: RNGptr
  type(isotropicMG),pointer :: mg

  integer(shortInt)          :: N, i
  !real(defReal)              :: Emax,Emin,Umax,Umin
  real(defReal)              :: k_old, k_new, ksum, ksum2, varK
  integer(shortInt)          :: nBins, idx
  integer(shortInt)          :: nInactive, nActive, startPop, endPop

  type(particleDungeon),pointer              :: cycle1, cycle2, cycleTemp
  integer(longInt), dimension(:),allocatable :: tally

  type(dictionary)      :: testDict
  type(IOdictionary)    :: IOdictTest

  ! ### Declarations end      ### !
  ! ### Main Programme Begins ### !

  call IOdictTest % initFrom('./RootMG')

  testDict = IOdictTest

  allocate(mg)
  call mg % init(testDict)

  allocate(RNGptr)
  call RNGptr % init(75785746574_longInt)

  !call collisionPhysics % attachXsData(mg)
  collisionPhysics % xsData => mg

  nBins = 3
  allocate(tally(nBins))
  tally = 0


  cycleTemp => null()
  nInactive = 300
  nActive   = 500
  N = 5000

  allocate(cycle1)
  allocate(cycle2)

  call cycle1 % init(int(3.0*N))
  call cycle2 % init(int(3.0*N))


 ! ##### Population initialisation
  neutron % pRNG => RNGptr
  neutron % xsData => mg
  do i=1,N
    neutron % G      = 1
    call neutron % teleport([0.0_8, 0.0_8, 0.0_8])
    call neutron % point([1.0_8, 0.0_8, 0.0_8])
    neutron % w      = 1.0
    neutron % isDead = .false.
    neutron % isMG   = .true.
    call cycle1 % detain(neutron)
  end do

! ##### Inactive cycles

  do i=1,nInactive
    startPop = cycle1 % popSize()
    generation: do

      call cycle1 % release(neutron)
      neutron % matIdx = 1

      History: do
        ! Tally energy
        !idx = 1 + int( nBins/(Umax-Umin) * (log(neutron % E) - Umin))
        !tally(idx) = tally(idx) + 1

        call collisionPhysics % collide(neutron,cycle1,cycle2)
        if(neutron % isDead) exit History

      end do History

     if(cycle1 % isEmpty() ) exit generation

    end do generation

   ! Calculate new k
    endPop = cycle2 % popSize()
    k_old  = cycle2 % k_eff
    k_new  = 1.0*endPop/startPop * k_old
   ! Normalise population
    call cycle2 % normSize(N, neutron % pRNG)

   ! Flip cycle dungeons
    cycleTemp => cycle2
    cycle2 => cycle1
    cycle1 => cycleTemp

   ! Load new k for normalisation
    cycle2 % k_eff = k_new

    print *, "Inactive cycle: ", i,"/",nInactive," k-eff (analog): ", k_new, "Pop: ", startPop, " -> ", endPop
  end do

! ************************************
! ****** Active cycles

  ksum  = 0.0
  ksum2 = 0.0
  varK = 0.0

  do i=1,nActive
    startPop = cycle1 % popSize()
    generationA: do

      call cycle1 % release(neutron)
      neutron % matIdx = 1

      HistoryA: do
        ! Tally energy
        idx = neutron % G
        tally(idx) = tally(idx) + 1
        call collisionPhysics % collide(neutron,cycle1,cycle2)
        if(neutron % isDead) exit HistoryA

      end do HistoryA

     if(cycle1 % isEmpty() ) exit generationA

    end do generationA

   ! Calculate new k
    endPop = cycle2 % popSize()
    k_old  = cycle2 % k_eff
    k_new  = 1.0*endPop/startPop * k_old

    ksum  = ksum  + k_new
    ksum2 = ksum2 + k_new * k_new

    k_new = ksum / i

   ! Normalise population
    call cycle2 % normSize(N, neutron % pRNG)

   ! Flip cycle dungeons
    cycleTemp => cycle2
    cycle2 => cycle1
    cycle1 => cycleTemp

   ! Load new k for normalisation
    cycle2 % k_eff = k_new

    if (i > 1 ) then
      varK = sqrt (1.0/(i*(i-1)) * (ksum2 - ksum*ksum/i))
    end if

    print *, "Active cycle: ", i,"/",nActive," k-eff (analog): ", k_new," +/- ", varK ," Pop: ", startPop, " -> ", endPop
  end do



end program eigenMG
