!!
!! Transport operator for implicit Monte Carlo tracking
!!
module transportOperatorIMC_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng
  use coord_class,                only : coordList

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
    real(defReal)                              :: majorant_inv, sigmaT, distance
    real(defReal)                              :: dTime, dGeom, dColl
    integer(shortInt)                          :: event
    type(coordList)                            :: p_coords
    character(100), parameter :: Here = 'IMCTracking (transportOperatorIMC_class.f90)' 

    ! Get majornat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getMajorantXS(p)

    IMCLoop:do

      ! Obtain the local cross-section
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

      ! Find distance to time boundary
      dTime = lightSpeed * (thisCycle % endOfStepTime - p % time)

      ! Find distance to cell boundary
      dGeom = 1000000
      p_coords = p % coords
      call self % geom % move_noCache(p % coords, dGeom, event)   ! Better way to do this?
      p % coords = p_coords 

      ! Sample distance to collision
      dColl = -log( p % pRNG % get() ) / sigmaT


      ! Find lowest value
      if ( dTime < dGeom .and. dTime < dColl) then
        print *, 'Time'
      else if ( dGeom < dColl ) then
        print *, 'Geom'
      else
        print *, 'Coll'
      end if

      !print *, 'dTime =', dTime, 'dGeom =', dGeom, 'dColl =', dColl


      distance = -log( p% pRNG % get() ) * majorant_inv

      ! Move partice in the geometry
      call self % geom % teleport(p % coords, distance)

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
        p % isDead = .true.
        return
      end if

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
