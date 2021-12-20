module IMCMGstd_class

  use numPrecision
  use endfConstants
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_NEUTRON
  use particleDungeon_class,         only : particleDungeon

  ! Abstract interface
  use collisionProcessor_inter,      only : collisionProcessor, collisionData ,init_super => init

  ! Nuclear Data Interface
  use nuclearDataReg_mod,            only : ndReg_getIMCMG => getIMCMG
  use nuclearDatabase_inter,         only : nuclearDatabase
  use mgIMCDatabase_inter,           only : mgIMCDatabase
  use mgIMCMaterial_inter,           only : mgIMCMaterial, mgIMCMaterial_CptrCast
  use reactionHandle_inter,          only : reactionHandle
  use multiScatterMG_class,          only : multiScatterMG, multiScatterMG_CptrCast
  !use fissionMG_class,               only : fissionMG, fissionMG_TptrCast

  ! Cross section packages
  use IMCXsPackages_class,       only : IMCMacroXSs


  ! Nuclear Data
  !use nuclearData_inter,              only : nuclearData
  !use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Cross-section packages to interface with nuclear data
  !use xsMacroSet_class,               only : xsMacroSet, xsMacroSet_ptr

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for MG IMCs
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary IMCs)
  !!   -> Does not create secondary non-IMC projectiles
  !!
  !! Settings:
  !!  NONE
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type            IMCMGstd;
  !!   }
  !!
  type, public, extends(collisionProcessor) :: IMCMGstd
    private
    class(mgIMCDatabase), pointer, public :: xsData => null()
    class(mgIMCMaterial), pointer, public :: mat    => null()
  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: elastic
    procedure :: inelastic
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs
  end type IMCMGstd

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(IMCMGstd), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(100), parameter :: Here = 'init (IMCMGstd_class.f90)'

    ! Call superclass
    call init_super(self, dict)

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(IMCMacroXSs)                :: macroXSs
    real(defReal)                        :: r
    character(100),parameter :: Here =' sampleCollision (IMCMGstd_class.f90)'

    ! Verify that particle is MG NEUTRON
    if( .not. p % isMG .or. p % type /= P_NEUTRON) then             ! (not yet integrated new particle type)
      call fatalError(Here, 'Supports only MG NEUTRON. Was given CE '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getIMCMG()
    if(.not.associated(self % xsData)) call fatalError(Here, "Failed to get active database for MG IMC")

    ! Get and verify material pointer
    self % mat => mgIMCMaterial_CptrCast( self % xsData % getMaterial( p % matIdx()))
    if(.not.associated(self % mat)) call fatalError(Here, "Failed to get MG IMC Material")

    ! Select Main reaction channel
    call self % mat % getMacroXSs(macroXSs, p % G, p % pRNG)
    r = p % pRNG % get()

    collDat % MT = macroXSs % invert(r)

  end subroutine sampleCollision

  !!
  !! Preform implicit treatment
  !!
  subroutine implicit(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(IMCMacroXSs)                :: macroXSs
    !type(fissionMG),pointer              :: fission
    type(particleState)                  :: pTemp
    real(defReal),dimension(3)           :: r, dir
    integer(shortInt)                    :: G_out, n, i
    real(defReal)                        :: wgt, w0, rand1, mu, phi
    real(defReal)                        :: sig_tot, k_eff, sig_nufiss
    character(100),parameter :: Here = 'implicit (IMCMGstd_class.f90)'

    !if ( self % mat % isFissile()) then
      ! Obtain required data
    !  wgt   = p % w                ! Current weight
    !  w0    = p % preHistory % wgt ! Starting weight
    !  k_eff = p % k_eff            ! k_eff for normalisation
    !  rand1 = p % pRNG % get()     ! Random number to sample sites

    !  call self % mat % getMacroXSs(macroXSs, p % G, p % pRNG)

    !  sig_tot    = macroXSs % total
    !  sig_nuFiss = macroXSs % nuFission

    !  ! Sample number of fission sites generated
    !  !n = int(wgt * sig_nuFiss/(sig_tot*k_eff) + r1, shortInt)
    !  n = int(abs( (wgt * sig_nuFiss) / (w0 * sig_tot * k_eff)) + rand1, shortInt)

    !  ! Shortcut if no particles were samples
    !  if (n < 1) return

    !  ! Get Fission reaction object
    !  fission => fissionMG_TptrCast( self % xsData % getReaction(macroFission, collDat % matIdx))
    !  if (.not.associated(fission)) call fatalError(Here, 'Failed to getrive fissionMG reaction object')

    !  ! Store new sites in the next cycle dungeon
    !  wgt =  sign(w0, wgt)
    !  r   = p % rGlobal()

    !  do i=1,n
    !    call fission % sampleOut(mu, phi, G_out, p % G, p % pRNG)
    !    dir = rotateVector(p % dirGlobal(), mu, phi)

    !    ! Copy extra detail from parent particle (i.e. time, flags ect.)
    !    pTemp       = p

    !    ! Overwrite position, direction, energy group and weight
    !    pTemp % r   = r
    !    pTemp % dir = dir
    !    pTemp % G   = G_out
    !    pTemp % wgt = wgt

    !    call nextCycle % detain(pTemp)
    !  end do
    !end if

  end subroutine implicit

  !!
  !! Elastic Scattering
  !!
  subroutine elastic(self, p , collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    ! Do nothing. Should not be called

  end subroutine elastic

  !!
  !! Preform scattering
  !!
  subroutine inelastic(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    class(multiScatterMG),pointer        :: scatter
    integer(shortInt)                    :: G_out   ! Post-collision energy group
    real(defReal)                        :: phi     ! Azimuthal scatter angle
    real(defReal)                        :: w_mul   ! Weight multiplier
    character(100),parameter :: Here = "inelastic (IMCMGstd_class.f90)"

    ! Assign MT number
    collDat % MT = macroIEscatter

    ! Get Scatter object
    scatter => multiScatterMG_CptrCast( self % xsData % getReaction(macroIEscatter, collDat % matIdx))
    if(.not.associated(scatter)) call fatalError(Here, "Failed to get scattering reaction object for MG IMC")

    ! Sample Mu and G_out
    call scatter % sampleOut(collDat % muL, phi, G_out, p % G, p % pRNG)

    ! Read scattering multiplicity
    w_mul = scatter % production(p % G, G_out)

    ! Update IMC state
    p % G = G_out
    p % w = p % w * w_mul
    call p % rotate(collDat % muL, phi)

  end subroutine inelastic

  !!
  !! Preform capture
  !!
  subroutine capture(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead = .true.

  end subroutine capture

  !!
  !! Preform fission
  !!
  subroutine fission(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

  !  p % isDead = .true.

  end subroutine fission

  !!
  !! Applay cutoffs or post-collision implicit treatment
  !!
  subroutine cutoffs(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    ! Do nothing

  end subroutine cutoffs

end module IMCMGstd_class
