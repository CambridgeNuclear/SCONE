module IMCMGstd_class

  use numPrecision
  use endfConstants
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_PHOTON
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

  ! Cross section packages
  use IMCXsPackages_class,           only : IMCMacroXSs


  ! Nuclear Data
  !use nuclearData_inter,              only : nuclearData
  !use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Cross-section packages to interface with nuclear data
  !use xsMacroSet_class,               only : xsMacroSet, xsMacroSet_ptr

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for MG IMC
  !! Determines type of collision as either absorption or effective scattering
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
    class(IMCMGstd), intent(inout)     :: self
    class(dictionary), intent(in)      :: dict
    character(100), parameter :: Here = 'init (IMCMGstd_class.f90)'

    ! Call superclass
    call init_super(self, dict)

  end subroutine init

  !!
  !! Samples collision
  !!
  !! Absorption with probability equal to fleck factor, otherwise
  !!  effective scattering
  !!
  !! Physical scattering is omitted as in reference paper "Four Decades of Implicit Monte Carlo"
  !!  (Allan B Wollaber) but may be included later if desired 
  !!
  subroutine sampleCollision(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)       :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(IMCMacroXSs)                    :: macroXSs
    real(defReal)                        :: r, fleck
    character(100),parameter :: Here =' sampleCollision (IMCMGstd_class.f90)'

    ! Verify that particle is MG PHOTON
    if( .not. p % isMG .or. p % type /= P_PHOTON) then
      call fatalError(Here, 'Supports only MG PHOTON. Was given NEUTRON and/or CE '//printType(p % type))
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

    fleck = self % mat % getFleck()

    if( r < fleck ) then
      ! Effective absoprtion
      collDat % MT = macroCapture
    else
      ! Effective scattering
      collDat % MT = macroIEScatter
    end if

  end subroutine sampleCollision

  !!
  !! Preform implicit treatment
  !!
  subroutine implicit(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)       :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    ! Do nothing.

  end subroutine implicit

  !!
  !! Elastic Scattering
  !!
  subroutine elastic(self, p , collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)       :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    character(100), parameter :: Here = 'elastic (IMCMGstd_class.f90)'

    ! Do nothing. Should not be called

    call fatalError(Here, "elastic subroutine should not be called")

  end subroutine elastic

  !!
  !! Preform scattering - Currently this is for effective scattering, and energy weights
  !!                       are unchanged (so is actually elastic)
  !!
  subroutine inelastic(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)       :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    class(multiScatterMG),pointer        :: scatter
    integer(shortInt)                    :: G_out   ! Post-collision energy group
    real(defReal)                        :: phi, mu     ! Azimuthal scatter angle
    real(defReal)                        :: w_mul   ! Weight multiplier
    real(defReal), dimension(3)          :: dir
    character(100),parameter :: Here = "inelastic (IMCMGstd_class.f90)"

    ! Assign MT number
    collDat % MT = macroIEScatter 

    ! Sample Direction - chosen uniformly inside unit sphere
    mu = 2 * p % pRNG % get() - 1
    phi = p % pRNG % get() * 2*pi
    dir(1) = mu
    dir(2) = sqrt(1-mu**2) * cos(phi)
    dir(3) = sqrt(1-mu**2) * sin(phi)

    !p % coords % dir = dir
    call p % rotate(mu, phi)

  end subroutine inelastic

  !!
  !! Preform capture
  !!
  subroutine capture(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)       :: self
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
    class(IMCMGstd), intent(inout)       :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    character(100), parameter :: Here = 'fission (IMCMGstd_class.f90)'

    ! Do nothing. Should not be called

    call fatalError(Here, "elastic subroutine should not be called")

  end subroutine fission

  !!
  !! Applay cutoffs or post-collision implicit treatment
  !!
  subroutine cutoffs(self, p, collDat, thisCycle, nextCycle)
    class(IMCMGstd), intent(inout)       :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    ! Do nothing

  end subroutine cutoffs

end module IMCMGstd_class
