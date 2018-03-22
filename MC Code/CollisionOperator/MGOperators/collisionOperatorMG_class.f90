module collisionOperatorMG_class

  use numPrecision
  use endfConstants

  use genericProcedures,      only : fatalError, rotateVector
  use RNG_class,              only : RNG
  use particle_class,         only : particle
  use particleDungeon_class,  only : particleDungeon
  use perMaterialMgXS_inter,  only : perMaterialMgXS

  ! Cross-section packages to interface with nuclear data
  use xsMacroSet_class,       only : xsMacroSet

  implicit none
  private

  type, public :: collisionOperatorMG
   !* private ** DEBUG
    class(perMaterialMgXS), pointer :: xsData => null()
    class(RNG), pointer             :: locRNG => null()

  contains
    procedure :: attachXSData
    procedure :: collide

    ! Private procedures

    procedure :: performScattering
    procedure :: performCapture
    procedure :: performFission
    procedure :: createFissionSites

  end type collisionOperatorMG

contains

  !!
  !! Initialise XS operator by providing pointer to XS data block
  !!
  subroutine attachXsData(self,xsData)
    class(collisionOperatorMG), intent(inout)   :: self
    class(perMaterialMgXS),pointer, intent(in)  :: xsData
    character(100), parameter               :: Here =' attachXsData (collisionOperatorMG_class.f90)'

    if(.not.associated(xsData)) call fatalError(Here,'Allocated xs data must be provided')

    self % xsData => xsData

  end subroutine attachXSData


  !!
  !! Subroutine to collide a neutron.
  !!
  !! Takes as an input two particleDungeons to be able to add new particles for this or next Cycle
  !!
  subroutine collide(self,p,thisCycle,nextCycle)
    class(collisionOperatorMG), intent(inout) :: self
    class(particle), intent(inout)            :: p
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    integer(shortInt)                         :: G
    integer(shortInt)                         :: matIdx
    integer(shortInt)                         :: MT
    real(defReal)                             :: r
    class(xsMacroSet), pointer                :: materialXS => null()

    ! Retrive Pointer to the random number generator
    self % locRNG => p % pRNG

    ! Load neutron energy group and material
    G = p % G
    matIdx = p % matIdx

    ! Select Main reaction channel
    call self % xsData % getMatMacroXS(materialXS,G,matIdx)

    r = self % locRNG % get()
    MT = materialXS % invert(r)

    ! Generate fission sites if nuclide is fissile
    if ( self % xsData % isFissileMat(matIdx)) then
      call self % createFissionSites(nextCycle,p,matIdx,materialXS)

    end if

    ! Call procedure to do reaction processing
    select case (MT)
      case(macroAllScatter)
       ! call self % performScattering(p,matIdx)

      case(macroCapture)
        call self % performCapture(p,matIdx)

      case(macroFission)
        call self % performFission(p,matIdx)

    end select


  end subroutine collide

  !!
  !! Change particle state in scattering reaction (mu is assumed to be in LAB frame)
  !!
  subroutine performScattering(self,p,matIdx)
    class(collisionOperatorMG), intent(inout) :: self
    type(particle), intent(inout)             :: p
    integer(shortInt), intent(in)             :: matIdx
    integer(shortInt)                         :: MT
    integer(shortInt)                         :: G         ! Pre-collision energy group
    integer(shortInt)                         :: G_out     ! Post-collision energy group
    real(defReal)                             :: muL       ! Cosine of scattering in LAB frame
    real(defReal)                             :: phi       ! Azimuthal scatter angle
    real(defReal)                             :: w_mul   ! Weight multiplier

    ! Assign MT number
    MT = macroAllScatter
    G = p % G

    ! Sample Mu and G_out
    call self % xsData % sampleMuGout(muL, G_out, G, self % locRNG, MT, matIdx)
    phi = TWO*PI * self % locRNG % get()

    ! Read scattering multiplicity
    w_mul = self % xsData % releaseAt(G, G_out, MT, matIdx)

    ! Update neutron state
    p % G = G_out
    p % w = p % w * w_mul
    call p % rotate(muL,phi)


  end subroutine performScattering




  !!
  !! Change particle state in capture reaction
  !!
  subroutine performCapture(self,p,matIdx)
    class(collisionOperatorMG), intent(inout) :: self
    type(particle), intent(inout)             :: p
    integer(shortInt),intent(in)              :: matIdx

    p % isDead =.true.

  end subroutine performCapture

  !!
  !! Change particle state in fission reaction
  !!
  subroutine performFission(self,p,matIdx)
    class(collisionOperatorMG), intent(inout) :: self
    type(particle), intent(inout)             :: p
    integer(shortInt),intent(in)              :: matIdx

    p % isDead =.true.

  end subroutine performFission

  !!
  !!
  !!
  subroutine createFissionSites(self,nextCycle,p,matIdx,Xss)
    class(collisionOperatorMG), intent(inout)  :: self
    type(particleDungeon), intent(inout)       :: nextCycle
    type(particle), intent(in)                 :: p
    integer(shortInt), intent(in)              :: matIdx
    type(xsMacroSet),pointer, intent(in)       :: Xss
    type(particle)                             :: pTemp
    real(defReal),dimension(3)                 :: r, dir
    integer(shortInt)                          :: G, G_out, n, i
    real(defReal)                              :: nu, wgt, r1, mu, phi
    real(defReal)                              :: sig_fiss, sig_tot, k_eff

    ! Obtain required data
    G     = p % G
    wgt   = p % w
    nu    = self % xsData % releaseAt(G,G, macroFission, matIdx)
    k_eff = nextCycle % k_eff
    r1    = self % locRNG % get()

    sig_fiss = Xss % fissionXS
    sig_tot  = Xss % totalXS

    r   = p % rGlobal()
    dir = p % globalDir()

    ! Sample number of fission sites generated
    n = int(wgt * nu * sig_fiss/(sig_tot*k_eff) + r1, shortInt)

    ! Throw new sites to the next cycle dungeon
    wgt = 1.0
    do i=1,n
      call self % xsData % sampleMuGout(mu, G_out, G, self % locRNG, macroFission, matIdx)
      phi = 2*PI * self % locRNG % get()
      dir = rotateVector(dir,mu,phi)

      call pTemp % build(r,dir,G_out,wgt)
      call nextCycle % throw(pTemp)

    end do

  end subroutine createFissionSites


end module collisionOperatorMG_class
