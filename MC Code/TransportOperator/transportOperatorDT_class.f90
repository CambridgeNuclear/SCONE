!!
!! Transport operator for delta tracking
!!
module transportOperatorDT_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError
  use particle_class,             only : particle
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

  type, public, extends(transportOperator) :: transportOperatorDT
  contains
    procedure :: init
    procedure :: transport => deltaTracking
  end type transportOperatorDT

contains

  !!
  !! Initialise transportOperatorDT
  !!
  subroutine init(self, nucData, geom, settings) !return nuclearData at some point!
    class(transportOperatorDT), intent(inout) :: self
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
  !!
  !!
  subroutine deltaTracking(self,p)
    class(transportOperatorDT), intent(in) :: self
    class(particle), intent(inout)         :: p
    real(defReal)                          :: majorant, sigmaT, distance

    majorant = self % nuclearData % getMajorantXS(p)

    ! Save pre-transition state
    call p % savePreTransition()

    DTLoop:do
      distance = -log( p% pRNG % get() )/majorant

      ! Move partice in the geometry
      call self % geom % teleport(p % coords, distance)

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % isDead = .true.
        ! TODO: REPORT HISTORY END
        return
      end if

      ! Obtain the local cross-section
      sigmaT = self % nuclearData % getTransXS(p, p % matIdx())

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real
      if (p % pRNG % get() < sigmaT/majorant) exit DTLoop

    end do DTLoop

    ! Tally transition
    call self % tally % reportTrans(p)

  end subroutine deltaTracking


end module transportOperatorDT_class
