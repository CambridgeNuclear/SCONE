!!
!! Transport operator for delta tracking
!!
module transportOperatorDT_class
  use numPrecision
  use genericProcedures
  use universalVariables
  use dictionary_class

  use transportOperator_inter

  use particle_class
  use surface_class
  use geometry_class
  use cell_class
  use rng_class

  !use Nuclear_Data_MG_class           !Re-instate later!!!
  !use Nuclear_Data_CE_class

  implicit none
  private

  type, public, extends(transportOperator) :: transportOperatorDT
  contains
    procedure :: init
    procedure :: transport => deltaTracking
    procedure :: applyBC
  end type transportOperatorDT

contains

  !!
  !! Initialise for delta tracking
  !!
  subroutine init(self, geom, settings) !return nuclearData at some point!
    class(transportOperatorDT), intent(inout) :: self
    !class(Nuclear_data_MG), target :: nuclearData
    class(geometry), target                   :: geom
    class(dictionary), optional               :: settings

    !self%MGData => nuclearData
    self%geom => geom

    ! TO DO: include settings, e.g., variance reduction, majorant adjustment

  end subroutine init

  subroutine deltaTracking(self,p)
    class(transportOperatorDT), intent(in) :: self
    class(particle), intent(inout)         :: p
    real(defReal)                          :: majorant, sigmaT, distance
    type(cell_ptr)                         :: currentCell

    !majorant = self % getMajorantXS(p)
    majorant = 0.8
    DTLoop:do
      distance = -log(p%pRNG%get())/majorant

      ! Move particle to new location in the global co-ordinate systems
      call p % moveGlobal(distance)

      ! Find the new cell which the particle occupies
      currentCell = self % geom % whichCell(p%coords)

      ! If the particle is outside the geometry, apply boundary conditions
      if (.not. currentCell % insideGeom()) then
        call self % applyBC(p, currentCell)
        ! End the transport step if the particle is killed
        if (p % isDead) then
          exit DTLoop
        end if
      end if

      ! Obtain the local cross-section
      call p % updateLocation()
      !sigmaT = self % getTotalXS(p)
      sigmaT = 0.6

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
