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

  !** Debug
  use coord_class,  only : coordList

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

    STLoop: do

      ! Obtain the local cross-section
      sigmaT = self % nuclearData % getTransXS(p, p % matIdx())

      ! Sample particle flight distance
      dist = -log( p % pRNG % get()) / sigmaT

      ! Move to the next stop
      call self % geom % move(p % coords, dist, isColl)

      ! * TODO SEND TALLY REPORT

      ! Return if particle stoped at collision (not cell boundary)
      if( isColl ) return

    end do STLoop

  end subroutine surfaceTracking

end module transportOperatorST_class
