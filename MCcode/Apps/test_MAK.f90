program test

  use numPrecision
  use dictionary_class,         only : dictionary
  use RNG_class,                only : RNG
  use IOdictionary_class,       only : IOdictionary
  use aceNeutronDatabase_class, only : aceNeutronDatabase
  use nuclearDatabase_inter,    only : nuclearDatabase
  use materialMenu_mod,         only : mm_init => init, mm_kill => kill
  use particle_class,           only : particle
  use ceNeutronCache_mod
  use ceNeutronDatabase_inter
  
  implicit none
!
!  ! Material definitions
!  character(*),parameter :: MAT_INPUT_STR = "  & 
!  & water { temp 273;                          &
!  &       composition {                        &
!  &       1001.03 5.028E-02;                   &
!  &       8016.03 2.505E-02;                   &
!  &       }                                    &
!  &     }                                      &
!  & uo2  { temp 1;                             &
!  &        composition {                       &
!  &        92233.03 2.286E-02;                 &
!  &        8016.03  4.572E-02;                 &
!  &        }                                   &
!  &      }"
!
!  ! CE Neutron Database specification
!  character(*),parameter :: ACE_INPUT_STR = " & 
!  & aceLibrary ./IntegrationTestFiles/testLib; "
!
!  class(ceNeutronDatabase), allocatable, target :: data
!  class(nuclearDatabase), pointer  :: ptr
!  type(IOdictionary)               :: matDict
!  type(IOdictionary)               :: dataDict
!  type(particle)                   :: p
!  type(RNG),target                 :: rand
!  real(defReal)                    :: t1, t2
!  real(defReal)                    :: Emin, Emax, E, step, L_low, L
!  integer(shortInt)                :: i, N
!
!  allocate(aceNeutronDatabase :: data)
!
!  ! Prepare dictionaries
!  call matDict % initFromChar(MAT_INPUT_STR)
!  call dataDict % initFromChar(ACE_INPUT_STR)
!
!  ! Build material menu
!  call mm_init(matDict)
!
!  ! Initialise data
!  ptr => data
!  select type(data)
!    type is (aceNeutronDatabase)
!      call data % init(dataDict, ptr, silent = .true.)
!      call data % activate([1,2])
!  end select
!
!  !call data % energyBounds (t1,t2)
!  !print *, t1, t2
!
!  call p % build([ZERO,ZERO,ZERO],[ONE,ZERO,ZERO], 1.0E-6_defReal, ONE)
!  p % pRNG => rand
!
!  !print *, data % materials(2) % nuclides
  !call data % updateTotalMatXS(1.0E-6_8, 1, rand)
  !print *, materialCache(1) % xss % total

  !call data % updateMajorantXS(1.0E-6_8, rand)
  !print *, majorantCache(1) % xs

  !print *, size(materialCache), size(nuclideCache), size(majorantCache)
!  Emin = 1.1E-11_8
!  Emax = 20.0_8
!  N = 20000
!
!  L_low = log(Emax/Emin)
!  step = L_low / N
!
!  print *, "A = ["
!  do i=1,N
!    L = L_low - step * i
!    E = Emax * exp(-L)
!    p % E = E
!    print *, E, ptr % getTotalMatXS(p,2)
!
!  end do
!  print *, "];"
!
!




end program test




