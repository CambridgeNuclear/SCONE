module neutronMGstd_class

  use numPrecision
  use endfConstants
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, phaseCoord, printType, P_NEUTRON
  use particleDungeon_class,         only : particleDungeon

  ! Abstract interface
  use collisionProcessor_inter,      only : collisionProcessor, collisionData ,init_super => init

  ! Nuclear Data
  use nuclearData_inter,              only : nuclearData
  use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Cross-section packages to interface with nuclear data
  use xsMacroSet_class,               only : xsMacroSet, xsMacroSet_ptr

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
  !!   type            neutronMGstd;
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronMGstd
    private
    class(perMaterialNuclearDataMG), pointer, public :: xsData => null()
  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: scatter
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs
  end type neutronMGstd

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronMGstd), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(100), parameter :: Here = 'init (neutronMGstd_class.f90)'

    ! Call superclass
    call init_super(self, dict)

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    real(defReal)                        :: r
    type(xsMacroSet_ptr)                 :: materialXS
    character(100),parameter :: Here =' sampleCollision (neutronMGstd_class.f90)'

    ! Verify that particle is MG neutron
    if( .not. p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Supports only MG Neutron. Was given CE '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    associate ( xs => p % xsData)
      select type(xs)
        class is(perMaterialNuclearDataMG)
          self % xsData => xs

        class default
          call fatalError(Here, 'Unsupported type of Nuclear Data interface. &
                               & Only perMaterialNuclearDataMG is accepted')
      end select
    end associate

    ! Select Main reaction channel
    call self % xsData % getMatMacroXS(materialXS, p % G, collDat % matIdx)

    r = p % pRNG % get()

    collDat % MT = materialXS % invert(r)

  end subroutine sampleCollision

  !!
  !! Preform implicit treatment
  !!
  subroutine implicit(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(xsMacroSet_ptr)                 :: materialXS
    type(phaseCoord)                     :: pTemp
    real(defReal),dimension(3)           :: r, dir
    integer(shortInt)                    :: G_out, n, i
    real(defReal)                        :: wgt, w0, rand1, mu, phi
    real(defReal)                        :: sig_tot, k_eff, sig_nufiss

    if ( self % xsData % isFissileMat(collDat % matIdx)) then
      ! Obtain required data
      wgt   = p % w                ! Current weight
      w0    = p % preHistory % wgt ! Starting weight
      k_eff = p % k_eff            ! k_eff for normalisation
      rand1 = p % pRNG % get()     ! Random number to sample sites

      call self % xsData % getMatMacroXS(materialXS, p % G, collDat % matIdx)

      sig_tot    = materialXS % totalXS()
      sig_nuFiss = materialXS % nuFissionXS()

      r   = p % rGlobal()
      dir = p % dirGlobal()

      ! Sample number of fission sites generated
      !n = int(wgt * sig_nuFiss/(sig_tot*k_eff) + r1, shortInt)
      n = int(abs( (wgt * sig_nuFiss) / (w0 * sig_tot * k_eff) + rand1), shortInt)

      ! Store new sites in the next cycle dungeon
      wgt =  sign(w0, wgt)

      do i=1,n
        call self % xsData % sampleMuGout(mu, G_out, p % G, p % pRNG, macroFission, collDat % matIdx)
        phi = TWO*PI * p % pRNG % get()
        dir = rotateVector(dir, mu, phi)

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
  !! Preform scattering
  !!
  subroutine scatter(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    integer(shortInt)                    :: G_out   ! Post-collision energy group
    real(defReal)                        :: phi     ! Azimuthal scatter angle
    real(defReal)                        :: w_mul   ! Weight multiplier

    ! Assign MT number
    collDat % MT = macroAllScatter

    ! Sample Mu and G_out
    call self % xsData % sampleMuGout(collDat % muL, G_out, p % G, p % pRNG, collDat % MT, collDat % matIdx)
    phi = TWO*PI * p % pRNG % get()

    ! Read scattering multiplicity
    w_mul = self % xsData % releaseAt(p % G, G_out, collDat % MT, collDat % matIdx)

    ! Update neutron state
    p % G = G_out
    p % w = p % w * w_mul
    call p % rotate(collDat % muL, phi)

  end subroutine scatter

  !!
  !! Preform capture
  !!
  subroutine capture(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGstd), intent(inout)   :: self
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
    class(neutronMGstd), intent(inout)   :: self
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
    class(neutronMGstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    ! Do nothing

  end subroutine cutoffs
    
end module neutronMGstd_class
