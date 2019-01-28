!!
!! Transport operator for delta tracking
!!
module transportOperatorST_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError
  use particle_class,             only : particle, phaseCoord
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Superclass
  use transportOperator_inter,    only : transportOperator

  ! Geometry interfaces
  use cellGeometry_inter,         only : cellGeometry

  ! Tally interface
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearData_inter,          only : nuclearData
  use transportNuclearData_inter, only : transportNuclearData

  implicit none
  private

  type, public, extends(transportOperator) :: transportOperatorST
  contains
    procedure :: init
    procedure :: transport => surfaceTracking
  end type transportOperatorST

contains

  !!
  !! Initialise transportOperatorST
  !!
  subroutine init(self, nucData, geom, settings)
    class(transportOperatorST), intent(inout) :: self
    class(nuclearData), pointer, intent(in)   :: nucData
    class(cellGeometry), pointer, intent(in)  :: geom
    class(dictionary), optional, intent(in)   :: settings
    character(100),parameter :: Here ='init (transportOperatorDT_class.f90)'

    ! Check that nuclear data type is supported
    select type(nucData)
      class is (transportNuclearData)
        self % nuclearData => nucData

      class default
        call fatalError(Here,'Class of provided nuclear data is not supported by transportOperator')

    end select

    ! Attach geometry
    self % geom => geom

    ! TO DO: include settings, e.g., variance reduction, majorant adjustment

  end subroutine init

  !!
  !! Performs surface tracking
  !!
  subroutine surfaceTracking(self,p)
    class(transportOperatorST), intent(in) :: self
    class(particle), intent(inout)         :: p
    logical(defBool)                       :: isColl
    real(defReal)                          :: sigmaT, dist

    ! Save pre-Transition state
    call p % savePreTransition()

    STLoop: do

      ! Obtain the local cross-section
      sigmaT = self % nuclearData % getTransXS(p, p % matIdx())

      ! Sample particle flight distance
      dist = -log( p % pRNG % get()) / sigmaT

      ! Save state before movement
      call p % savePrePath()

      ! Move to the next stop
      call self % geom % move(p % coords, dist, isColl)

      ! Send tally report
      call self % tally % reportPath(p, dist)

      ! Return if particle stoped at collision (not cell boundary)
      if( isColl ) exit STLoop

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % isDead = .true.
        ! TODO: REPORT HISTORY END
        exit STLoop
      end if
    end do STLoop

    ! Send transition report
    call self % tally % reportTrans(p)

  end subroutine surfaceTracking

end module transportOperatorST_class
