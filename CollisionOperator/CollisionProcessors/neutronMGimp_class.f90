module neutronMGimp_class

  use numPrecision
  use endfConstants
  use universalVariables,            only : nameWW
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_NEUTRON
  use particleDungeon_class,         only : particleDungeon

  ! Abstract interface
  use collisionProcessor_inter,      only : collisionProcessor, collisionData ,init_super => init

  ! Nuclear Data Interface
  use nuclearDataReg_mod,            only : ndReg_getNeutronMG => getNeutronMG
  use nuclearDatabase_inter,         only : nuclearDatabase
  use mgNeutronDatabase_inter,       only : mgNeutronDatabase
  use mgNeutronMaterial_inter,       only : mgNeutronMaterial, mgNeutronMaterial_CptrCast
  use reactionHandle_inter,          only : reactionHandle
  use multiScatterMG_class,          only : multiScatterMG, multiScatterMG_CptrCast
  use fissionMG_class,               only : fissionMG, fissionMG_TptrCast

  ! Cross section packages
  use neutronXsPackages_class,       only : neutronMacroXSs

  ! Geometry and fields
  use geometryReg_mod,                only : gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use weightWindowsField_class,       only : weightWindowsField, weightWindowsField_TptrCast


  ! Nuclear Data
  !use nuclearData_inter,              only : nuclearData
  !use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Cross-section packages to interface with nuclear data
  !use xsMacroSet_class,               only : xsMacroSet, xsMacroSet_ptr

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for MG neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  NONE
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type            neutronMGimp;
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronMGimp
    private
    class(mgNeutronDatabase), pointer, public :: xsData => null()
    class(mgNeutronMaterial), pointer, public :: mat    => null()

    ! Variance reduction options
    integer(shortInt) :: maxSplit
    logical(defBool)  :: weightWindows
    type(weightWindowsField), pointer :: weightWindowsMap

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

    ! Variance reduction procedures
    procedure, private :: split
    procedure, private :: russianRoulette

  end type neutronMGimp

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronMGimp), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    integer(shortInt)                  :: idx
    character(100), parameter :: Here = 'init (neutronMGimp_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    ! Obtain settings for variance reduction
    call dict % getOrDefault(self % maxSplit,'maxSplit', 1000)
    call dict % getOrDefault(self % weightWindows,'weightWindows', .false.)

    ! Sets up the weight windows field
    if (self % weightWindows) then
      idx = gr_fieldIdx(nameWW)
      self % weightWindowsMap => weightWindowsField_TptrCast(gr_fieldPtr(idx))
    end if

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMacroXSs)                :: macroXSs
    real(defReal)                        :: r
    character(100),parameter :: Here =' sampleCollision (neutronMGimp_class.f90)'

    ! Verify that particle is MG neutron
    if( .not. p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Supports only MG Neutron. Was given CE '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronMG()
    if(.not.associated(self % xsData)) call fatalError(Here, "Failed to get active database for MG Neutron")

    ! Get and verify material pointer
    self % mat => mgNeutronMaterial_CptrCast( self % xsData % getMaterial( p % matIdx()))
    if(.not.associated(self % mat)) call fatalError(Here, "Failed to get MG Neutron Material")

    ! Select Main reaction channel
    call self % mat % getMacroXSs(macroXSs, p % G, p % pRNG)
    r = p % pRNG % get()

    collDat % MT = macroXSs % invert(r)

  end subroutine sampleCollision

  !!
  !! Preform implicit treatment
  !!
  subroutine implicit(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMacroXSs)                :: macroXSs
    type(fissionMG),pointer              :: fission
    type(particleState)                  :: pTemp
    real(defReal),dimension(3)           :: r, dir
    integer(shortInt)                    :: G_out, n, i
    real(defReal)                        :: wgt, w0, rand1, mu, phi
    real(defReal)                        :: sig_tot, k_eff, sig_nufiss
    character(100),parameter :: Here = 'implicit (neutronMGimp_class.f90)'

    if ( self % mat % isFissile()) then
      ! Obtain required data
      wgt   = p % w                ! Current weight
      w0    = p % preHistory % wgt ! Starting weight
      k_eff = p % k_eff            ! k_eff for normalisation
      rand1 = p % pRNG % get()     ! Random number to sample sites

      call self % mat % getMacroXSs(macroXSs, p % G, p % pRNG)

      sig_tot    = macroXSs % total
      sig_nuFiss = macroXSs % nuFission

      ! Sample number of fission sites generated
      !n = int(wgt * sig_nuFiss/(sig_tot*k_eff) + r1, shortInt)
      n = int(abs( (wgt * sig_nuFiss) / (w0 * sig_tot * k_eff)) + rand1, shortInt)

      ! Shortcut if no particles were samples
      if (n < 1) return

      ! Get Fission reaction object
      fission => fissionMG_TptrCast( self % xsData % getReaction(macroFission, collDat % matIdx))
      if (.not.associated(fission)) call fatalError(Here, 'Failed to getrive fissionMG reaction object')

      ! Store new sites in the next cycle dungeon
      wgt =  sign(w0, wgt)
      r   = p % rGlobal()

      do i=1,n
        call fission % sampleOut(mu, phi, G_out, p % G, p % pRNG)
        dir = rotateVector(p % dirGlobal(), mu, phi)

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        pTemp       = p

        ! Overwrite position, direction, energy group and weight
        pTemp % r   = r
        pTemp % dir = dir
        pTemp % G   = G_out
        pTemp % wgt = wgt

        call nextCycle % detain(pTemp)
      end do
    end if

  end subroutine implicit

  !!
  !! Elastic Scattering
  !!
  subroutine elastic(self, p , collDat, thisCycle, nextCycle)
    class(neutronMGimp), intent(inout)   :: self
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
    class(neutronMGimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    class(multiScatterMG),pointer        :: scatter
    integer(shortInt)                    :: G_out   ! Post-collision energy group
    real(defReal)                        :: phi     ! Azimuthal scatter angle
    real(defReal)                        :: w_mul   ! Weight multiplier
    character(100),parameter :: Here = "inelastic (neutronMGimp_class.f90)"

    ! Assign MT number
    collDat % MT = macroIEscatter

    ! Get Scatter object
    scatter => multiScatterMG_CptrCast( self % xsData % getReaction(macroIEscatter, collDat % matIdx))
    if(.not.associated(scatter)) call fatalError(Here, "Failed to get scattering reaction object for MG neutron")

    ! Sample Mu and G_out
    call scatter % sampleOut(collDat % muL, phi, G_out, p % G, p % pRNG)

    ! Read scattering multiplicity
    w_mul = scatter % production(p % G, G_out)

    ! Update neutron state
    p % G = G_out
    p % w = p % w * w_mul
    call p % rotate(collDat % muL, phi)

  end subroutine inelastic

  !!
  !! Preform capture
  !!
  subroutine capture(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGimp), intent(inout)   :: self
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
    class(neutronMGimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead = .true.

  end subroutine fission

  !!
  !! Applay cutoffs or post-collision implicit treatment
  !!
  subroutine cutoffs(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    real(defReal), dimension(3)          :: val
    real(defReal)                        :: minWgt, maxWgt, avWgt

    if (p % isDead) then
      ! Do nothing !

    ! Weight Windows treatment
    elseif (self % weightWindows) then
      val = self % weightWindowsMap % at(p)
      minWgt = val(1)
      maxWgt = val(2)
      avWgt  = val(3)

      ! If a particle is outside the WW map and all the weight limits
      ! are zero nothing happens. NOTE: this holds for positive weights only
      if ((p % w > maxWgt) .and. (maxWgt /= ZERO) .and. (p % splitCount < self % maxSplit)) then
        call self % split(p, thisCycle, maxWgt)
      elseif (p % w < minWgt) then
        call self % russianRoulette(p, avWgt)
      end if

    end if

  end subroutine cutoffs

  !!
  !! Perform Russian roulette on a particle
  !!
  subroutine russianRoulette(self, p, avWgt)
    class(neutronMGimp), intent(inout) :: self
    class(particle), intent(inout)     :: p
    real(defReal), intent(in)          :: avWgt

    if (p % pRNG % get() < (ONE - p % w/avWgt)) then
      p % isDead = .true.
    else
      p % w = avWgt
    end if

  end subroutine russianRoulette

  !!
  !! Split particle which has too large a weight
  !!
  subroutine split(self, p, thisCycle, maxWgt)
    class(neutronMGimp), intent(inout)    :: self
    class(particle), intent(inout)        :: p
    class(particleDungeon), intent(inout) :: thisCycle
    real(defReal), intent(in)             :: maxWgt
    integer(shortInt)                     :: mult, i, splitCount

    ! This value must be at least 2
    mult = ceiling(p % w/maxWgt)

    ! Limit maximum split
    if (mult > self % maxSplit - p % splitCount + 1) then
      mult = self % maxSplit  - p % splitCount + 1
    end if

    ! Decrease weight
    p % w = p % w/mult

    ! Save current particle splitCount
    splitCount = p % splitCount

    ! Add split particle's to the dungeon
    do i = 1,mult-1
      p % splitCount = 0
      call thisCycle % detain(p)
    end do

    p % splitCount = splitCount + mult

  end subroutine split

end module neutronMGimp_class
