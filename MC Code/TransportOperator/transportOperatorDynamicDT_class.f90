!!
!! Transport operator for delta tracking
!!
module transportOperatorDynamicDT_class
  use numPrecision
  use universalVariables
  use tallyCodes,                 only : leak_FATE, aged_FATE

  use genericProcedures,          only : fatalError
  use particle_class,             only : particle
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Superclass
  use transportOperator_inter,    only : transportOperator

  ! Geometry interfaces
  use cellGeometry_inter,         only : cellGeometry

  ! Nuclear data interfaces
  use nuclearData_inter,          only : nuclearData
  use transportNuclearData_inter, only : transportNuclearData

  implicit none
  private

  type, public, extends(transportOperator) :: transportOperatorDynamicDT
  contains
    procedure :: init
    procedure :: transport => deltaTracking
  end type transportOperatorDynamicDT

contains

  !!
  !! Initialise transportOperatorDT
  !!
  subroutine init(self, nucData, geom, settings) !return nuclearData at some point!
    class(transportOperatorDynamicDT), intent(inout) :: self
    class(nuclearData), pointer, intent(in)          :: nucData
    class(cellGeometry), pointer, intent(in)         :: geom
    class(dictionary), optional, intent(in)          :: settings
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

  subroutine deltaTracking(self,p)
    class(transportOperatorDynamicDT), intent(in) :: self
    class(particle), intent(inout)         :: p
    real(defReal)                          :: majorant, sigmaT, distance, &
                                              timeMax, flightTime, speed

    majorant = self % nuclearData % getMajorantXS(p)
    timeMax = p % timeMax
    speed = lightSpeed * sqrt(TWO * p % E / neutronMass)

    DTLoop:do
      distance = -log(p%pRNG%get())/majorant

      ! Calculate particle flight time
      flightTime = distance / speed

      ! Does particle cross time boundary?
      if (flightTime + p % time >= timeMax) then
        flightTime = timeMax - p % time
        distance = speed * flightTime
        p % time = timeMax
        p % fate = aged_FATE
      else
        p % time = p % time + flightTime
      end if

      ! Move partice in the geometry
      call self % geom % teleport(p % coords, distance)

      ! If particle has leaked exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % isDead = .true.
        p % fate = leak_FATE
        ! TODO: REPORT HISTORY END
        return
      end if

      ! Check whether particle has reached the time boundary
      if ( p % fate == aged_FATE) then
        p % isDead = .true.
        exit DTLoop
      end if

      ! Obtain the local cross-section
      sigmaT = self % nuclearData % getTransXS(p, p % matIdx())

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real
      if (p%pRNG%get() < sigmaT/majorant) exit DTLoop

    end do DTLoop

  end subroutine deltaTracking

end module transportOperatorDynamicDT_class
