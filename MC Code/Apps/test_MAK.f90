program test

  use numPrecision
  use RNG_class
  use genericProcedures
  use byNucNoMT_Data_class, only : byNucNoMT_Data
  use byNucNoMT_class,      only : byNucNoMT
  use aceNoMT_class

  use endfTable_class, only : endfTable

  use releaseLawENDF_class
  use constantRelease_class, only : constantRelease
  use polynomialRelease_class, only : polynomialRelease
  use tabularRelease_class, only : tabularRelease

  use muEndfPdf_class,   only : muEndfPdf, muEndfPdf_ptr
  use isotropicmu_class, only : isotropicmu
  use equiBin32mu_class, only : equiBin32mu
  use tabularmu_class,   only : tabularmu
  use tabularPdf_class,    only : tabularPdf
  use angleLawENDF_class,   only : angleLawENDF
  use tabularAngle_class, only : tabularAngle
  use noAngle_class,       only: noAngle

  use tabularEnergy_class, only: tabularEnergy
  use contTabularEnergy_class, only : contTabularEnergy
  use energyLawENDF_class,     only : energyLawENDF
  use noEnergy_class,          only : noEnergy

  use dictionary_class , only: dictionary, dictContent
  use IOdictionary_class, only : IOdictionary

  use maxwellEnergyPdf_class, only: maxwellEnergyPdf
  use maxwellSpectrum_class,  only: maxwellSpectrum

  use xsCDF_class,     only : xsCDF
  use xsMainCDF_class, only : xsMainCDF
  use xsMainSet_class, only : xsMainSet

  use xsEnergyPointNoMT_class, only : xsEnergyPointNoMT

  use xsNucMacroSet_class,     only : xsNucMacroSet
  use xsMacroSet_class, only : xsMacroSet

!  use collisionOperator_class, only : collisionOperator
  use particle_class,          only : particle

  use scatteringKernels_func,  only : targetVelocity_constXS,asymptoticScatter
  use particleDungeon_class,   only : particleDungeon
  use outscatterCDF_class,     only : outscatterCDF

  use perMaterialMgXs_inter, only : perMaterialMgXs
  use isotropicMG_class,       only : isotropicMG
 ! use collisionOperatorMG_class, only : collisionOperatorMG

  use datalessMaterials_class,     only : datalessMaterials

  implicit none

  type myType
    real(defReal) :: c
    real(defReal) :: a
    real(defReal) :: b
  end type myType


  integer(kind=shortInt) :: i
  integer(longInt)       :: longI, longI2, rate
  integer(kind=shortInt),dimension(99),target :: A
  integer(kind=shortInt),dimension(:),pointer :: B
  INTEGER(SHORTiNT),DIMENSION(:), allocatable :: C
  type(ByNucNoMT_Data)  :: CEdata
  character(len=pathLen)      :: matInput="./testInput"
  character(len=pathLen)      :: isoInput="/home/mak60/myACE/JEF311.aceXS"
  character(len=99)      :: format
  character(len=99),dimension(2) :: Ach
  character(len=pathLen)    :: acePath = " /home/mak60/myACE/acedata/92235JEF311.ace "
  character(pathLen)       :: testDictFile = "./testDictInput"
  integer(shortInt)         :: firstLine = 1170286
  type(aceNoMT)             :: isotope
  real(defReal) :: kl, eps, km, kn
  real(defReal),dimension(:),allocatable :: x, pdf, x2, pdf2,R

  type(RNG) :: random


  class(releaseLawENDF),pointer  :: release
  class(endfTable),pointer       :: eTable

  real(defReal), dimension(20) ::  second
  integer(shortInt), dimension(4) :: bounds, ENDF

  real, pointer :: p1,p2,p3

  type(muEndfPdf_ptr) :: myPtr, myPtr2
  type(tabularEnergy), pointer :: tabPtr
  class(angleLawENDF),pointer :: angle
  class(energyLawENDF), pointer :: energyT
  type(tabularEnergy),dimension(:),allocatable   :: tables
  type(contTabularEnergy) :: yjfj

  class(*),pointer   :: GP(:), GP2(:)
  real,pointer       :: a1(:),a2(:)
  real,pointer       :: rp1(:),rp2(:)
  character(:),pointer :: char_ptr

  type(dictionary) :: testDict
  type(dictionary) :: testDict2
  type(dictionary) :: testDict3
  type(dictionary) :: testDict4

  character(nameLen) :: abc = "aaa", bc
  character(20),dimension(:),allocatable  :: charT
  character(20),dimension(:),pointer      :: cA_ptr
  class(*),dimension(:),pointer           :: Uptr
  class(*),dimension(:),pointer           :: Uptr2
  character(20),dimension(:),pointer      :: cA_ptr2
  character(20),dimension(:),pointer      :: cA_ptr3
  character(20),dimension(:),pointer      :: cA_ptr4
  character(20),dimension(:),pointer      :: localPointer
  character(20),dimension(:), pointer     :: newMemory

  character(20),dimension(:),allocatable :: cA_alloc1
  character(20),dimension(:),allocatable :: cA_alloc2
  type(dictContent)     :: dictNode
  type(IOdictionary)    :: IOdictTest
  integer(shortInt) :: err
  character(100) :: msg
  type(maxwellEnergyPdf) :: maxwell
  type(maxwellSpectrum), pointer :: maxSpec

  type(myType),target   :: testType
  real(defReal),pointer :: tt

  type(xsMainCDF) :: myCDF
  type(xsMainCDF) :: uCDF
  type(xsMainCDF) :: bCDF


  class(xsMainSet), pointer :: mySet
  type(xsMacroSet) :: bSet
  type(xsMainSet) :: uSet
  type(xsEnergyPointNoMT) :: Epoint
  type(xsEnergyPointNoMT) :: Bpoint
  type(xsEnergyPointNoMT) :: Upoint

  type(byNucNoMT),pointer :: ce
  real(defReal),dimension(:),allocatable :: energy
  integer(shortInt)                       :: N

  type(xsNucMacroSet),pointer  :: nuclideInvert
  type(xsMacroSet),pointer :: MacroXS

  type(particle)          :: neutron
!  type(collisionOperator) :: collisionPhysics
!  type(collisionOperatorMG) :: MGcoll
  class(RNG), pointer     :: RNGptr
  real(defReal),dimension(3) :: V_ter
  real(defReal) :: mu, E_as, E_fg
  real(defReal) :: Emax,Emin,Umax,Umin
  integer(shortInt) :: nBins, idx
  integer(longInt), dimension(:),allocatable :: tally
  type(particleDungeon),pointer :: cycle1, cycle2, cycleTemp
  integer(shortInt)     :: nInactive, nActive, startPop, endPop
  real(defReal)         :: k_old, k_new, ksum, ksum2, varK

  type(outscatterCDF) :: outCDF
  integer(shortInt), dimension(:,:), allocatable :: ReSh
  class(perMaterialMgXS),allocatable :: MGData
  type(datalessMaterials)  :: matNoDat


!  call IOdictTest % initFrom('./materialInput')
!
!  testDict = IOdictTest
!
!
!  call matNoDat % init(testDict)
!  print *, matNoDat % materials
!
!  print *, matNoDat % getName(4)
!
!
!
!  stop
!**********************************************************************!

!  bSet % total    = 4.0
!  bSet % scatter  = 3.0
!  bSet % capture  = 1.5
!  bSet % fission  = 0.5
!
!  uSet % total   = 6.0
!  uSet % scatter = 3.5
!  uSet % capture = 1.5
!  uSet % fission = 1.0
!
!  call uCDF % init(uSet % scatter, uSet % capture, uSet % fission)
!  call bCDF % init(bSet % scatter, bSet % capture, bSet % fission)
!
! ! print *, uCDF % cdf
! ! call myCDF % interpolate( bCDF,uCDF, 0.5_8)
! ! print *, myCDF % cdf
! ! print *, bCDF % cdf
!
!  call mySet % interpolate(bSet,uSet,1.0_8)
!
!  print *, bSet % total, bSet % scatter, bSet % capture, bSet % fission
!  print *, mySet % total, mySet % scatter, mySet % capture, mySet % fission
!  print *, uSet % total, uSet % scatter, uSet % capture, uSet % fission
!  call myCDF % init (1.0_8,2.0_8,0.1_8)
!! print *, myCDF % cdf
!! print *, myCDF % invert(0.0_8)
!  call bPoint % init(2.0_8,1.0_8)
!  call uPoint % init(3.0_8,1.0_8,9.9_8)
!
!  call ePoint % interpolate (bPoint,uPoint,-1.0_8)
!
!  print*,  bPoint % xs % total,bPoint % xs % scatter, bPoint % xs % capture, bPoint % xs % fission
!  print*,  ePoint % xs % total,ePoint % xs % scatter, ePoint % xs % capture, ePoint % xs % fission
!  print*,  uPoint % xs % total,uPoint % xs % scatter, uPoint % xs % capture, uPoint % xs % fission
!  call isotope % init(acePath,243050)
!
!  do i = 1,size(isotope % energyGrid)
!    print *,isotope % energyGrid(i) ,isotope % xsData(i) % xs % total,isotope % xsData(i) % xs % scatter, &
!            isotope % xsData(i) % xs % capture, isotope % xsData(i) % xs % fission
!  end do
!!
!call IOdictTest % initFrom('./RootMG')
!testDict = IOdictTest
!
!allocate( isotropicMG :: MGData)
!
!call MGData % init(testDict)
!print *, MGData % matData
!stop
!Set % scatterXS = 2.0
!Set % captureXS = 1.0
!Set % fissionXS = 0.5
!Set % totalXS = bSet % scatterXS + bSet % captureXS + bSet % fissionXS
!print *, bSet % invert(0.9_8)
!print *, binarySearch([1.0_8],1.1_8)
!all outCDF % init([2.0_8,1.0_8,0.5_8])
!rint *, outCDF % invert(0.85_8)
!eSh = reshape([1, 2, 3, 4, 5, 6,7,8,9],[3,3])
!rint *, ReSh(:,2)
!rint *, sum(ReSh,1)
!top
!all IOdictTest % initFrom('./materialInput')
!print *, IOdictTest % keysDict()
!estDict = IOdictTest !% getDict('myFourthMat')
!print *, testDict % getChar('type')
!print *, testDict % keysReal()
!print *, testDict % getReal('temp')
!charT = testDict % keysReal()
!print *, charT
!print *, IOdictTest % keysDict_type('material')
!llocate(ce)
!llocate(ce % dataBlock)
!all ce % dataBlock % init(testDict)
!all ce % dataBlock % printMe()
!top
!!
!allocate(ce)
!call ce % readFrom(matInput,isoInput)
!allocate(RNGptr)
!call RNGptr % init(5875757_8)
!neutron % pRNG => RNGptr
!neutron % E = 7.0
!!neutron % dir = [1.0, 0.0 , 0.0]
!neutron % matIdx = 4
!neutron % isDead = .false.
!call collisionPhysics % attachXsData(ce)
!collisionPhysics % locRNG => neutron % pRNG
!call RNGptr % init(75785746574_longInt)
!Emax = 20.0
!Emin = 1.0E-11
!Umax = log(Emax)
!Umin = log(Emin)
!nBins = 300
!!N = 1000000
! N = 5000
!allocate(tally(nBins))
!tally = 0
!allocate(cycle1)
!allocate(cycle2)
!call cycle1 % init(int(2.0*N))
!call cycle2 % init(int(2.0*N))
!cycleTemp => null()
!nInactive = 300
!nActive   = 2000
! ##### Population initialisation
!do i=1,N
!  neutron % E      = 0.5
!  call neutron % teleport([0.0_8, 0.0_8, 0.0_8])
!  call neutron % point([1.0_8, 0.0_8, 0.0_8])
!  neutron % w      = 1.0
!  neutron % isDead = .false.
!  call cycle1 % throw(neutron)
!end do
!##### Fixed Source Calculation
!###########################################################
!do
!!   neutron % E      = 10.0
!!   call neutron % point([1.0_8, 0.0_8, 0.0_8])
!!   neutron % matIdx = 4
!!   neutron % isDead = .false.
!   call cycle1 % release(neutron)
!   neutron % matIdx = 4
!
!
!   History: do
!     ! Tally energy
!     idx = 1 + int( nBins/(Umax-Umin) * (log(neutron % E) - Umin))
!     tally(idx) = tally(idx) + 1
!
!     call collisionPhysics % collide(neutron,cycle1,cycle1)
!     if(neutron % isDead) exit History
!   end do History
!
!   if(cycle1 % isEmpty() ) exit
!end do
!##### Eigenvalue calculation
!####################################################
!  *** Inactive cycles
! do i=1,nInactive
!   startPop = cycle1 % popSize()
!   generation: do
!     call cycle1 % release(neutron)
!     neutron % matIdx = 4
!     History: do
!       ! Tally energy
!       !idx = 1 + int( nBins/(Umax-Umin) * (log(neutron % E) - Umin))
!       !tally(idx) = tally(idx) + 1
!       call collisionPhysics % collide(neutron,cycle1,cycle2)
!       if(neutron % isDead) exit History
!     end do History
!    if(cycle1 % isEmpty() ) exit generation
!   end do generation
!  ! Calculate new k
!   endPop = cycle2 % popSize()
!   k_old  = cycle2 % k_eff
!   k_new  = 1.0*endPop/startPop * k_old
!  ! Normalise population
!   call cycle2 % normSize(N, neutron % pRNG)
!  ! Flip cycle dungeons
!   cycleTemp => cycle2
!   cycle2 => cycle1
!   cycle1 => cycleTemp
!  ! Load new k for normalisation
!   cycle2 % k_eff = k_new
!   print *, "Inactive cycle: ", i,"/",nInactive," k-eff (analog): ", k_new, "Pop: ", startPop, " -> ", endPop
! end do
! ************************************
! ****** Active cycles
! ksum  = 0.0
! ksum2 = 0.0
! varK = 0.0
! do i=1,nActive
!   startPop = cycle1 % popSize()
!   generationA: do
!     call cycle1 % release(neutron)
!     neutron % matIdx = 4
!     HistoryA: do
!       ! Tally energy
!       idx = 1 + int( nBins/(Umax-Umin) * (log(neutron % E) - Umin))
!       tally(idx) = tally(idx) + 1
!       call collisionPhysics % collide(neutron,cycle1,cycle2)
!       if(neutron % isDead) exit HistoryA
!     end do HistoryA
!    if(cycle1 % isEmpty() ) exit generationA
!   end do generationA
!  ! Calculate new k
!   endPop = cycle2 % popSize()
!   k_old  = cycle2 % k_eff
!   k_new  = 1.0*endPop/startPop * k_old
!  ksum  = ksum  + k_new
!  ksum2 = ksum2 + k_new * k_new
!  k_new = ksum / i
!  ! Normalise population
!   call cycle2 % normSize(N, neutron % pRNG)
!  ! Flip cycle dungeons
!   cycleTemp => cycle2
!   cycle2 => cycle1
!   cycle1 => cycleTemp
!  ! Load new k for normalisation
!   cycle2 % k_eff = k_new
!   if (i > 1 ) then
!     varK = sqrt (1.0/(i*(i-1)) * (ksum2 - ksum*ksum/i))
!   end if
!   print *, "Active cycle: ", i,"/",nActive," k-eff (analog): ", k_new," +/- ", varK ," Pop: ", startPop, " -> ", endPop
! end do
!rint *, 'S = ['
!o i =1,size(tally)
! print *, tally(i)
!nd do
!rint *, '];'
!stop
!print *,"S=["
!do i=1,5000000
!  neutron % E = 0.06e-6_8
!
!  call collisionPhysics % scatterFromFixed(mu,neutron, neutron % E ,0.999170_8,2,1)
! !mu = 2.0* RNGptr % get() - 1.0
!  E_as = neutron % E
!  !call asymptoticScatter(E_as,mu,236.005800_8)
!  neutron % E = 0.06e-6_8
!
!  call collisionPhysics % scatterFromMoving(mu,neutron, neutron % E ,0.999170_8,7.755597306E-8_8,2,1)
!  E_fg = neutron % E
!  print *, E_as, E_fg
!
! ! V_ter = targetVelocity_constXS(0.001E-6_8, [1.0_8,0.0_8,0.0_8] ,2.585199E-8_8, 300.0_8, RNGptr)
! !print *, sqrt(sum(V_ter *V_ter))
!! print '(E100.90)', sqrt(V_ter(1)**2+V_ter(2)**2+V_ter(3)**2)
!
!end do
!print *,"];"
! do i=1,100
!
!   call collisionPhysics % collide(neutron)
!
! end do
!!call ce % dataBlock % print()
!top
!****************************************************************************
! ***** Test play code to interpolate XSs
!
!      ! Set energy points for interpolation
!      ! min exponent -7
!      ! max exponent 1.30
!      N = 100000
!      allocate (energy(N))
!      energy = [((1.30+7.0)/N * i - 7.0, i=1,N)]
!      energy = 10**energy
!      !energy = [(20.0/N * i, i=1,N)]
!
!
!      print *, ce % isInCMframe(18,2)
!
!      do i=1,N
!         !kl = ce % matShelf(1) % getTotal(energy(i))
!         !nuclideInvert => ce % matShelf(1) % nucCdf
!
!         print *, energy(i), ce % getMajorantXS(energy(i)),ce % getTotalMatXS(energy(i),1), &
!                  ce % getTotalMatXS(energy(i),2),ce % getTotalMatXS(energy(i),3),ce % getTotalMatXS(energy(i),4)
!
!
!         !print *, energy(i)
!         !call ce % matShelf(1) % setEnergy(energy(i))
!         !MacroXS => ce % matShelf(1) % XS
!         !print *, energy(i), MacroXS % totalXS, MacroXS % scatterXS, MacroXS % captureXS, MacroXS % fissionXS
!         !print *, energy(i), ce % matShelf(1) % getTotal(energy(i))
!         !call ce % getMainXS(mySet,energy(i),4)
!         !rint *, energy(i), mySet % total, mySet % scatter, mySet % capture, mySet % fission
!      end do
!      stop
!*********************************
! !call CEdata % readFrom(matInput,isoInput)
! !call CEdata % print()
!  stop
!
!  testType % a = 9.9
!  testType % b = 1.1
!
!  tt => testType % a
!
!  print *, tt
!  !deallocate(tt)
!
!  print *, testType % b
!! stop
!  call IOdictTest % initFrom(testDictFile)

!  print *, "TOP DICTIONARY"
!  print *, IOdictTest % getRealArray('list1')
!  print *, IOdictTest % getCharArray('list3')
!  print *, IOdictTest % getReal('keyword1')
!  print *, IOdictTest % getReal('keyword2')
!  print *, IOdictTest % getChar('keyword3')
!  print *, IOdictTest % getInt('keyword4')
!  print *, IOdictTest % getInt('keyword7')
!
!  testDict = IOdictTest % getDict('dictionary1')
!
!  print *, "NESTED DICTIONARY "
!  print *, testDict % getInt('keyword')
!  print *, testDict % getReal('keyword2')
!  print *, testDict % getChar('keyword3')
! print *, testDict % getIntArray('list2')
!
!  testDict2 = testDict % getDict('dictionary3')
!
!  print *, "NESTED DICTIONARY2 "
!  print *, testDict2 % getInt('key')
!  print *, testDict2 % getInt('key2')
!
!  testDict3 = testDict2 % getDict('dictionary4')
!
!  print *, "NESTED DICTIONARY3"
!  print *, testDict3 % getInt('key')
!  print *, testDict3 % getReal('key2')
!
!  testDict4 = testDict2 % getDict('dictionary5')
!
!  print *, "NESTED DICTONARY 3.1"
!  print *, testDict4 % getReal('key')
!  print *, testDict4 % getInt('key2')
!print *, IOdictTest % keysReal()
!print *, IOdictTest % keysInt()
!print *, IOdictTest % keysChar()
!print *, IOdictTest % keysRealArray()
!print *, IOdictTest % keysIntArray()
!print *, IOdictTest % keysCharArray()
!print *, IOdictTest % keysDict()
! print *, IOdictTest % getCharArray('list3')
!do i=1,100000
!  !print *, i
!  call testDict % kill()
!  call IOdictTest % get(testDict,'dictionary1')
!end do
!stop

!  maxSpec => maxwellSpectrum( [1.0E-11_8, 20.0_8], [1.33974_8, 1.33974_8], -20.0_8)
!
!
!  do i =1,80000
!   ! kl = 20.0/80000 * (i-1)
!    print *, maxSpec % sample(1.0_8,random)
!  end do
!
!stop

!  call isotope % init(acePath,1)

  call CEdata % readFrom(matInput,isoInput)
  call CEdata % printMe()

  stop
  !call isotope % init(acePath,firstLine)


  !R= [1.1_8, 2.3_8, 3.6_8, 9.6_8, 11.9_8]

  !print *, binaryFloorIdxC(R,1.1_8)

  !C = [(2*i,i=1,10)]



  energy = [(real(i),i=1,20)]
  second = [(i*5.0,i=1,20)]
  bounds = [ 5, 10, 15, 20]
  ENDF = [5,1,3,4]
   !bounds = [20]
   !ENDF = [1]


!  second = [1.0_8,1.5_8,2.0_8,1.5_8,1.7_8,2.5_8,2.3_8,2.0_8,2.1_8,1.5_8]
!  energy = [(real(i),i=1,10)]
!  bounds = [3,7,10]
!  ENDF = [1,2,1]
!  !energy(10)= energy(9)
!  energy(4) = energy(3)
!  energy(6) = energy(5)
!

!  release => tabularRelease(energy,second,bounds,ENDF)!
!
!  eTable  => endfTable(energy,second,bounds,ENDF)
!
!  do i=1,100
!    kl = (19.99-1.0)/100.0 * i + 1.0
!    print *, kl, release % releaseAt(kl), eTable % at(kl), release % releaseAt(kl)-eTable % at(kl)
!  end do


!  x = [(-1.0+2.0/10*i,i=0,10)]
!  pdf = abs(x)
!
!  x2   = [ -1.0_8, 1.0_8]
!  pdf2 = [ 0.5_8, 0.5_8]
!
!  x = [ -1.0_8, 0.0_8, 1.0_8]
!  pdf = [ 1.0_8/3 ,2.0_8/3, 76876.0_8]

  !call table % init(x,pdf,0)


!
!  R(1) = -1.0
!  myPtr = equiBin32mu(R)
!  myPtr  = tabularmu(x2,pdf2,1)
!  myPtr2 = tabularmu(x,pdf,1)
!
!  allocate(tables(2))
!  call tables(1) % init (x2,pdf2,1)
!  call tables(2) % init (x,pdf,1)
!
!
!  angle  => tabularAngle([0.0_8, 1.0_8],[myPtr, myPtr2])
!  energyT => contTabularEnergy([0.0_8, 1.0_8],tables)
!
!  deallocate(energyT)
!
!  energyT => noEnergy()
!  angle => noAngle()
! ! myPtr = tabularmu(x,pdf,0)
!
!  do i=0,1000
!    kl = 2.0/1000 * i - 1.0
!     !eps = random % get()
!
!  !   print *, kl, angle % probabilityOf(kl,1.0_8), energyT % probabilityOf(kl,0.5_8)
!
!  end do



  !print *, binarySearch(energy,1.0_8)

  !print *, linearFloorIdx(C,21)
!  call system_clock(count_rate=rate)
!  call system_clock(count=longI)
!  do i=1,1
!    kl= interpolate(0.00001_8, 1.0_8, 0.01_8, 10.0_8, 0.5_8)
!  end do
!  call system_clock(count=longI2)
!  print *, 'end', (longI2-longI)/real(rate)
end program test

