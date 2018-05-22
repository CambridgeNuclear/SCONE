program eigenCE

  use numPrecision
  use RNG_class,                         only : RNG
  use byNucNoMT_class,                   only : byNucNoMT
  use perNuclideNuclearDataCE_inter,     only : perNuclideNuclearDataCE
  !use collisionOperator_class,           only : collisionOperator
  use perNuclideCollisionOpCE_class,     only : perNuclideCollisionOpCE
  use perNuclideImplicitCaptureCE_class, only : perNuclideImplicitCaptureCE
  use particle_class,                    only : particle, phaseCoord
  use particleDungeon_class,             only : particleDungeon

  use dictionary_class ,       only : dictionary
  use IOdictionary_class,      only : IOdictionary

  use tallyAdminBase_class,    only : tallyAdminBase
  use keffActiveClerk_class,   only : keffActiveClerk
  use keffInactiveClerk_class, only : keffInactiveClerk
  use tallyInactiveAdmin_class, only : tallyInactiveAdmin
  use tallyActiveAdmin_class,   only : tallyActiveAdmin

  use tallyClerkFactory_func,  only : new_tallyClerk, new_tallyClerk_ptr

  implicit none

  type(particle)          :: neutron
  !type(collisionOperator) :: collisionPhysics_implement
  class(RNG), pointer     :: RNGptr
  class(byNucNoMT),pointer :: ce_implement

  !type(perNuclideImplicitCaptureCE) :: collisionPhysics
  type(perNuclideCollisionOpCE) :: collisionPhysics
  class(perNuclideNuclearDataCE),pointer :: ce


  integer(shortInt)          :: N, i
  real(defReal)              :: Emax,Emin,Umax,Umin
  real(defReal)              :: k_old, k_new, ksum, ksum2, varK
  real(defReal)              :: leakProb, r1
  integer(shortInt)          :: nBins, idx
  integer(shortInt)          :: nInactive, nActive, startPop, endPop

  type(particleDungeon),pointer              :: cycle1, cycle2, cycleTemp
  integer(longInt), dimension(:),allocatable :: tally
  type(phaseCoord) :: pre

  type(dictionary)      :: testDict
  type(dictionary)      :: clerkDict
  type(IOdictionary)    :: IOdictTest

  !type(tallyAdminBase)  :: tallyActive
  !type(tallyAdminBase)  :: tallyInactive
  type(tallyInactiveAdmin) :: tallyInactive
  type(tallyActiveAdmin) :: tallyActive
  type(keffActiveClerk)    :: k_imp
  type(keffInactiveClerk)  :: k_ana

  !### Declarations end
  !### Main Programme Begins



  call IOdictTest % initFrom('./materialInput')

  testDict = IOdictTest

  allocate(ce_implement)
  call ce_implement % init(testDict)

  allocate(RNGptr)
  call RNGptr % init(75785746574_longInt)



 ! allocate(collisionPhysics)
 ! call collisionPhysics_implement % attachXsData(ce_implement)

  ce => ce_implement

  !***** Create Tallies

  call clerkDict % init(1)
  call clerkDict % store('type','keffActiveClerk ')

  call tallyInactive % init()
  call tallyActive % init()

!  call tallyInactive % addTallyClerk(new_tallyClerk_ptr(clerkDict) )

 ! call tallyActive % addTallyClerk(new_tallyClerk(clerkDict))

  collisionPhysics % xsData => ce

   print *, 'Here'



  Emax = 20.0
  Emin = 1.0E-11
  Umax = log(Emax)
  Umin = log(Emin)

  nBins = 300
 !N = 1000000
  N = 5000
  allocate(tally(nBins))
  tally = 0
  leakProb = 0.0005_8
  leakProb = 0.0

  allocate(cycle1)
  allocate(cycle2)

  call cycle1 % init(int(100.0*N))
  call cycle2 % init(int(100.0*N))
  cycleTemp => null()
  nInactive = 300
  nActive   = 500

! ##### Population initialisation
  neutron % pRNG => RNGptr
  neutron % xsData => ce
  do i=1,N
    neutron % E      = 0.5
    call neutron % teleport([0.0_8, 0.0_8, 0.0_8])
    call neutron % point([1.0_8, 0.0_8, 0.0_8])
    neutron % w      = 1.0
    neutron % isDead = .false.
    call cycle1 % detain(neutron)
  end do


! ##### Inactive cycles

  do i=1,nInactive
    startPop = cycle1 % popSize()
    !*** Send report to tally
    call tallyInactive % reportCycleStart(cycle1)
    generation: do

      call cycle1 % release(neutron)
      neutron % matIdx = 4

      History: do
        ! Save beginning of history info
        pre = neutron

        ! Tally energy
        !idx = 1 + int( nBins/(Umax-Umin) * (log(neutron % E) - Umin))
        !tally(idx) = tally(idx) + 1

        ! Check if leaked
        r1 = RNGptr % get()
        if( r1 < leakProb ) then ! Neutron has leaked
          call tallyInactive % reportHist(pre,neutron,5001)
          exit History

        end if

        !** Send report to tally
        call tallyInactive  % reportInColl(neutron)
        !**
        call collisionPhysics % collide(neutron,cycle1,cycle2)
        if(neutron % isDead) then
          call tallyInactive % reportHist(pre,neutron,5000)
          exit History
        end if
      end do History

     if(cycle1 % isEmpty() ) exit generation

    end do generation

    !*** Send report to tally
    call tallyInactive % reportCycleEnd(cycle2)
    !***

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
    call tallyInactive % display()
  end do

! ************************************
! ****** Active cycles

  ksum  = 0.0
  ksum2 = 0.0
  varK = 0.0

  do i=1,nActive
    startPop = cycle1 % popSize()
    !*** Send report to tally
    call tallyActive % reportCycleStart(cycle1)
          !***
    generationA: do

      call cycle1 % release(neutron)
      neutron % matIdx = 4

      HistoryA: do
        ! Save beginning of history info
        pre = neutron

        ! Tally energy
        idx = 1 + int( nBins/(Umax-Umin) * (log(neutron % E) - Umin))
        tally(idx) = tally(idx) + 1

        ! Check if leaked
        r1 = RNGptr % get()
        if( r1 < leakProb ) then ! Neutron has leaked
          call tallyActive % reportHist(pre,neutron,5001)
          exit HistoryA

        end if

        !** Send report to tally
        call tallyActive % reportInColl(neutron)
        !**
        call collisionPhysics % collide(neutron,cycle1,cycle2)
        if(neutron % isDead) then
          call tallyActive % reportHist(pre,neutron,5000)
          exit HistoryA
        end if

      end do HistoryA

     if(cycle1 % isEmpty() ) exit generationA

    end do generationA

    !*** Send report to tally
    call tallyActive % reportCycleEnd(cycle2)
    !***


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
    call tallyActive % display()
  end do





print *, 'S = ['
do i =1,size(tally)
  print *, tally(i)
end do

print *, '];'

print *,"SAMPLES COUNT: ", neutron % pRNG % getCount()
print *,"Number of Collisions: ", collisionPhysics % collCount

end program eigenCE
