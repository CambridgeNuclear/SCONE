!!
!! Transport operator for implicit Monte Carlo scheme using delta tracking
!!
module transportOperatorIMC_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle, P_PHOTON
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Superclass
  use transportOperator_inter,    only : transportOperator

  ! Geometry interfaces
  use geometry_inter,             only : geometry

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase
  use mgIMCMaterial_inter,        only : mgIMCMaterial, mgIMCMaterial_CptrCast

  implicit none
  private

  !!
  !! Transport operator that moves a particle with IMC tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorIMC
    class(mgIMCMaterial), pointer, public :: mat    => null()
  contains
    procedure          :: transit => imcTracking
    procedure, private :: materialTransform
  end type transportOperatorIMC

contains

  subroutine imcTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorIMC), intent(inout) :: self
    class(particle), intent(inout)             :: p
    type(tallyAdmin), intent(inout)            :: tally
    class(particleDungeon), intent(inout)      :: thisCycle
    class(particleDungeon), intent(inout)      :: nextCycle
    real(defReal)                              :: majorant_inv, sigmaT, dTime, dColl
    character(100), parameter :: Here = 'IMCTracking (transportOperatorIMC_class.f90)' 

    ! Get majornat XS inverse: 1/Sigma_majorant
    !majorant_inv = ONE / self % xsData % getMajorantXS(p)

    IMCLoop:do

      ! Deal with material particles, only relevant for ISMC
      if(p % getType() == P_MATERIAL_MG) then
        call self % materialTransform(p, tally)
        if(p % fate == TIME_FATE) exit IMCLoop
      end if

      if(p % getType() .ne. P_PHOTON_MG) call fatalError(Here, 'Particle is not MG Photon')

      majorant_inv = ONE / self % xsData % getMajorantXS(p)

      ! Find distance to time boundary
      dTime = lightSpeed * (p % timeMax - p % time)

      ! Sample distance to move particle before potential collision
      dColl = -log( p % pRNG % get() ) * majorant_inv

      ! Determine which distance to move particle
      if (dColl < dTime) then
        ! Move partice to potential collision location
        call self % geom % teleport(p % coords, dColl)
        p % time = p % time + dColl / lightSpeed
      else
        ! Move particle to end of time step location
        call self % geom % teleport(p % coords, dTime)
        p % fate = TIME_FATE
        p % time = p % timeMax
      end if

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
        p % isDead = .true.
        return
      end if

      if (p % fate == TIME_FATE) exit IMCLoop

      ! Check for void
      if( p % matIdx() == VOID_MAT) cycle IMCLoop

      ! Obtain the local cross-section
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

      ! Protect Against Sillines
      !if( sigmaT*majorant_inv < ZERO .or. ONE < sigmaT*majorant_inv) then
      !  call fatalError(Here, "TotalXS/MajorantXS is silly: "//numToChar(sigmaT*majorant_inv))
      !end if

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real
      if (p % pRNG % get() < sigmaT*majorant_inv) exit IMCLoop

    end do IMCLoop

    call tally % reportTrans(p)
  end subroutine imcTracking


  !!
  !! Transform material particles into radiation photons with
  !! probability per unit time of c*sigma_a*fleck*eta
  !!
  !! Used only for ISMC, not for standard IMC
  !!
  subroutine materialTransform(self, p, tally)
    class(transportOperatorIMC), intent(inout) :: self
    class(particle), intent(inout)             :: p
    type(tallyAdmin), intent(inout)            :: tally
    real(defReal)                              :: sigmaT, fleck, eta, mu, phi
    real(defReal), dimension(3)                :: dir
    character(100), parameter                  :: Here = 'materialTransform (transportOperatorIMC_class.f90)'
 
    ! Confirm that time = 0
    !if (p % time .ne. 0) call fatalError(Here, 'Material particle should have time = 0')

    ! Get and verify material pointer
    self % mat => mgIMCMaterial_CptrCast( self % xsData % getMaterial( p % matIdx()))
    if(.not.associated(self % mat)) call fatalError(Here, "Failed to get MG IMC Material")

    sigmaT = self % xsData % getTransMatXS(p, p % matIdx())     !! Should be sigma_a, may need changing when sorting out cross-sections
    fleck  = self % mat % getFleck()
    eta    = self % mat % getEta()

    ! Sample time to transform into radiation photon
    p % time = p % time - log(p % pRNG % get()) / (sigmaT*fleck*eta*lightSpeed)

    ! Exit loop if particle remains material until end of time step
    if (p % time >= p % timeMax) then
      p % fate = TIME_FATE
      p % time = p % timeMax
      ! Tally energy for next temperature calculation
      call tally % reportHist(p)
    else
      p % type = P_PHOTON
      ! Resample direction
      mu = 2 * p % pRNG % get() - 1
      phi = p % pRNG % get() * 2*pi
      dir(1) = mu
      dir(2) = sqrt(1-mu**2) * cos(phi)
      dir(3) = sqrt(1-mu**2) * sin(phi)
      call p % point(dir)
    end if

  end subroutine materialTransform


end module transportOperatorIMC_class
