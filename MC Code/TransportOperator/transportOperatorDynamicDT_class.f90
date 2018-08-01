!!
!! Transport operator for delta tracking
!!
module transportOperatorDynamicDT_class
  use numPrecision
  use universalVariables
  use tallyCodes                  only : leak_FATE, aged_FATE

  use genericProcedures,          only : fatalError
  use particle_class,             only : particle
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Superclass
  use transportOperator_inter,    only : transportOperator

  ! Geometry interfaces
  use geometry_class,             only : geometry
  use surface_class,              only : surface_ptr
  use cell_class,                 only : cell_ptr

  ! Nuclear data interfaces
  use nuclearData_inter,          only : nuclearData
  use transportNuclearData_inter, only : transportNuclearData

  implicit none
  private

  type, public, extends(transportOperator) :: transportOperatorDynamicDT
  contains
    procedure :: init
    procedure :: transport => deltaTracking
    procedure :: applyBC
  end type transportOperatorDynamicDT

contains

  !!
  !! Initialise transportOperatorDT
  !!
  subroutine init(self, nucData, geom, settings) !return nuclearData at some point!
    class(transportOperatorDynamicDT), intent(inout) :: self
    class(nuclearData), pointer, intent(in)          :: nucData
    class(geometry), pointer, intent(in)             :: geom
    class(dictionary), optional, intent(in)          :: settings
    character(100),parameter :: Here ='init (transportOperatorDynamicDT_class.f90)'

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
    class(transportOperatorDT), intent(in) :: self
    class(particle), intent(inout)         :: p
    real(defReal)                          :: majorant, sigmaT, distance, &
                                              tmax, time, flightTime, speed
    type(cell_ptr)                         :: currentCell

    majorant = self % nuclearData % getMajorantXS(p)
    tmax = p % timeMax

    DTLoop:do
      distance = -log(p%pRNG%get())/majorant

      ! Calculate particle flight time
      speed = lightSpeed * sqrt(TWO * p % E / neutronMass)
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

      ! Move particle to new location in the global co-ordinate systems
      call p % moveGlobal(distance)

      ! Find the new cell which the particle occupies
      currentCell = self % geom % whichCell(p%coords)

      ! If the particle is outside the geometry, apply boundary conditions
      if (.not. currentCell % insideGeom()) then
        call self % applyBC(p, currentCell)
        ! End the transport step if the particle is killed
        if (p % isDead) then
          p % fate = leak_FATE
          exit DTLoop
        end if
      end if

      ! Check whether particle has reached the time boundary
      if ( p % fate == aged_FATE) then
        p % isDead
        exit DTLoop
      end if

      ! Obtain the local cross-section
      call p % updateLocation()
      sigmaT = self % nuclearData % getTransXS(p, p% matIdx)

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real
      if (p%pRNG%get() < sigmaT/majorant) then
        exit DTLoop
      end if
    end do DTLoop

    call currentCell % kill()

  end subroutine deltaTracking

  !!
  !! Apply boundary conditions when using delta tracking
  !! Given location in a particular halfspace, surface should
  !! communicate the necessary transformations required
  !!
  subroutine applyBC(self, p, currentCell)
    class(transportOperatorDT), intent(in) :: self
    class(particle), intent(inout)         :: p
    class(cell_ptr), intent(inout)         :: currentCell
    type(surface_ptr)                      :: currentSurface
    logical(defBool)                       :: killParticle

    killParticle = .FALSE.

    ! Iterate until particle is either inside the geometry or dead
    do while (.NOT. (currentCell % insideGeom() .OR. p % isDead))

      ! Find the bounding surface
      currentSurface = self % geom % boundarySurface

      ! Reset nesting and apply appropriate transformations to the particle
      call p % resetNesting()
      call currentSurface % boundaryTransform(p%coords%lvl(1)%r, p%coords%lvl(1)%dir, killParticle)

      ! Kill particle if it crossed a vacuum boundary
      if(killParticle) then
        p % isDead = .TRUE.
      else
        ! Identify which cell the particle now occupies
        currentCell = self % geom % whichCell(p%coords)
      end if

    end do

  end subroutine applyBC

end module transportOperatorDT_class
