!!
!! Transport operator for implicit Monte Carlo scheme using delta tracking
!!
module transportOperatorIMC_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle
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

  implicit none
  private

  !!
  !! Transport operator that moves a particle with IMC tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorIMC
  contains
    procedure :: transit => imcTracking
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
    majorant_inv = ONE / self % xsData % getMajorantXS(p)

    IMCLoop:do

      ! Find distance to time boundary
      dTime = lightSpeed * (p % timeMax - p % time)

      ! Sample distance to move particle before potential collision
      dColl = -log( p% pRNG % get() ) * majorant_inv

      ! Determine which distance to move particle
      if (dColl < dTime) then
        ! Move partice to potential collision location
        call self % geom % teleport(p % coords, dColl)
        p % time = p % time + dColl / lightSpeed
      else
        ! Move particle to end of time step location
        call self % geom % teleport(p % coords, dTime)
        p % fate = AGED_FATE
        p % time = p % timeMax
      end if

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
        p % isDead = .true.
        return
      end if

      if (p % fate == AGED_FATE) exit IMCLoop

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


end module transportOperatorIMC_class
