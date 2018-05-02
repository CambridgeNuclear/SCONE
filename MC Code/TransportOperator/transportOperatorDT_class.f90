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
  subroutine init(self, random, geom, settings) !return nuclearData at some point!
    class(transportOperatorDT), intent(inout) :: self
    !class(Nuclear_data_MG), target :: nuclearData
    class(rng), target                        :: random
    class(geometry), target                   :: geom
    class(dictionary), optional               :: settings

    self%random => random
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
      distance = -log(self%random%get())/majorant

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
      !sigmaT = self % getTotalXS(p)
      sigmaT = 0.6

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real
      if (self%random%get() < sigmaT/majorant) then
        exit DTLoop
      end if
    end do DTLoop

    call currentCell % kill()

  end subroutine deltaTracking

  !!
  !! Apply boundary conditions when using delta tracking
  !!
  subroutine applyBC(self, p, currentCell)
    class(transportOperatorDT), intent(in) :: self
    class(particle), intent(inout)         :: p
    class(cell_ptr), intent(inout)         :: currentCell
    type(surface_ptr)                      :: currentSurface

    ! Iterate until particle is either inside the geometry or dead
    do while (.NOT. currentCell % insideGeom())
      ! Identify which surface the particle crossed at the highest geometry level
      currentSurface = currentCell % whichSurface(p%rGlobal(), -p%dirGlobal())

      ! Check the boundary conditions on the surface
      ! If vacuum, kill the particle
      if (currentSurface % isVacuum()) then
        p % isDead = .TRUE.

      ! If reflective or periodic, must ensure that the surface is a plane!
      else if (currentSurface % isReflective()) then

        ! Return particle to global coordinates and apply the reflective transform
        call p % resetNesting()
        call currentSurface % reflectiveTransform(p%coords%lvl(1)%r, p%coords%lvl(1)%dir)

        ! Identify which cell the particle now occupies
        currentCell = self % geom % whichCell(p%coords)

      else if (currentSurface % isPeriodic()) then

        ! Apply the periodic translation associated with the surface to global co-ordinates
        call p % teleport(p % rGlobal() + currentSurface % periodicTranslation())

        ! Identify which cell the particle now occupies
        currentCell = self % geom % whichCell(p%coords)

      else
        call fatalError('applyBC, transportOperatorDT',&
        'Could not identify correct boundary conditions')
      end if

    end do

    ! The particle should no longer be sat on a surface
    call currentSurface % kill()

  end subroutine applyBC

  !!
  !! Apply boundary conditions when using delta tracking
  !! Given location in a particular halfspace, surface should
  !! communicate the necessary transformations required
  !!
  subroutine applyBCTransformation(self, p, currentCell)
    class(transportOperatorDT), intent(in) :: self
    class(particle), intent(inout)         :: p
    class(cell_ptr), intent(inout)         :: currentCell
    type(surface_ptr)                      :: currentSurface
    logical(defBool)                       :: killParticle

    killParticle = .FALSE.

    ! Iterate until particle is either inside the geometry or dead
    do while (.NOT. (currentCell % insideGeom() .AND. p % isDead))
      ! Identify which surface the particle crossed at the highest geometry level
      currentSurface = currentCell % whichSurface(p%rGlobal(), -p%dirGlobal())

      ! Check the boundary conditions on the surface
      ! If vacuum, kill the particle
      if (currentSurface % isVacuum()) then
        p % isDead = .TRUE.
        exit

      ! For compound surfaces, apply the transformation corresponding
      ! to the particle's location in halfspace
      else if (currentSurface % isCompound) then

        call p % resetNesting()
        call currentSurface % applyTransformation()

        if(killParticle) then
          p % isDead = .TRUE.
        else
          ! Identify which cell the particle now occupies
          currentCell = self % geom % whichCell(p%coords)
        end if

      else
        call fatalError('applyBC, transportOperatorDT',&
        'Could not identify boundary conditions to apply')

      end if

    end do

  end subroutine applyBCTransformation

end module transportOperatorDT_class
