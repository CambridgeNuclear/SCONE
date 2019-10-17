!!
!! Transport operator for delta tracking
!!
module transportOperatorST_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use RNG_class,                  only : RNG

  ! Superclass
  use transportOperator_inter,    only : transportOperator

  ! Geometry interfaces
  use cellGeometry_inter,         only : cellGeometry

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !!
  !! Transport operator that moves a particle with delta tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorST
  contains
    procedure :: transit => surfaceTracking
  end type transportOperatorST

contains

  !!
  !! Performs surface tracking
  !!
  subroutine surfaceTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorST), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    logical(defBool)                          :: isColl
    real(defReal)                             :: sigmaT, dist

    STLoop: do

      ! Obtain the local cross-section
      if( p % matIdx() == VOID_MAT) then
        dist = INFINITY

      else
        sigmaT = self % xsData % getTransMatXS(p, p % matIdx())
        dist = -log( p % pRNG % get()) / sigmaT
      end if

      ! Save state before movement
      call p % savePrePath()

      ! Move to the next stop. NOTE: "move" resets dist to distanced moved!
      call self % geom % move(p % coords, dist, isColl)

      ! Send tally report for a path moved
      call tally % reportPath(p, dist)

      ! Kill particle if it has leaked
      if( p % matIdx() == OUTSIDE_FILL) then
        p % isDead = .true.
        p % fate = LEAK_FATE
      end if

      ! Return if particle stoped at collision (not cell boundary)
      if( isColl .or. p % isDead) exit STLoop

    end do STLoop

    call tally % reportTrans(p)

  end subroutine surfaceTracking

end module transportOperatorST_class
