!!
!! Transport operator for delta tracking
!!
module transportOperatorST_class
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
  use lattice_class

  !use Nuclear_Data_MG_class           !Re-instate later!!!
  !use Nuclear_Data_CE_class

  implicit none
  private

  type, public, extends(transportOperator) :: transportOperatorST
  contains
    procedure :: init
    procedure :: transport => surfaceTracking
    procedure :: applyBC
  end type transportOperatorST

contains

  !!
  !! Initialise for delta tracking
  !!
  subroutine init(self, geom, settings)
    class(transportOperatorST), intent(inout) :: self
    !class(Nuclear_data_MG), target :: nuclearData
    class(geometry), target                 :: geom
    class(dictionary), optional             :: settings

    !self%MGData => nuclearData
    self%geom => geom

    ! TO DO: include settings, e.g., variance reduction

  end subroutine init

  !!
  !! Performs surface tracking
  !!
  subroutine surfaceTracking(self,p)
    class(transportOperatorST), intent(in) :: self
    class(particle), intent(inout)         :: p
    real(defReal)                          :: sigmaT, distance, boundaryDistance, &
                                              testDistance, latDistance, lDist
    integer(shortInt)                      :: i, n, nMin, latIdx, iLat
    type(cell_ptr)                         :: c, currentCell
    type(lattice_ptr)                      :: lat
    logical(defBool)                       :: moveUp

    STLoop: do

      ! Calculate boundary distance: descend the different co-ordinate levels starting from the highest
      ! Ensures bounds of parent cells are not exceeded
      n = p % nesting()
      nMin = 1
      boundaryDistance = INFINITY
      latDistance = INFINITY
      do i = 1,n
        c = self % geom % cells(p % getCellIdx(i))
        testDistance = c % getDistance(p%rLocal(i), p%dirLocal(i))

        ! Check if the particle is in a lattice cell
        latIdx = p % getLatIdx(i)
        if (latIdx > 0) then
          lat = self % geom % lattices(latIdx)
          lDist = lat % getDistance(p%rLocal(i), p%dirLocal(i))
          if (latDistance >= lDist) then
            latDistance = lDist
            iLat = i - 1
          end if
        end if

        if (boundaryDistance > testDistance) then
          boundaryDistance = testDistance
          nMin = i
          ! Must move up a level if crossing a lattice boundary
          if (boundaryDistance > latDistance) then
            nMin = iLat
            boundaryDistance = latDistance
          end if
        end if
      end do
      call c % kill()
      call lat % kill()
      n = nMin

      ! Obtain the local cross-section
      call p % updateLocation()
      !sigmaT = self % getTotalXS(p)
      sigmaT = 0.6

      ! Sample particle flight distance
      distance = -log(p%pRNG%get())/sigmaT

      ! The particle escapes the cell and moves to the next
      if (boundaryDistance <= distance) then

        ! Move particle to the surface with a small nudge to move across the boundary
        call p % moveLocal(boundaryDistance + NUDGE, n)

        ! Find the new base cell which the particle occupies
        currentCell = self % geom % whichCell(p%coords, n)

        ! If the particle is outside the geometry, apply boundary conditions
        do while (.not. currentCell % insideGeom())
          print *,'boundary'
          call self % applyBC(p, currentCell)
          ! End transport if the particle is killed
          if (p % isDead) then
            exit STLoop
          end if
        end do

        ! Continue performing surface tracking
        cycle

      else
        ! Move particle to new location
        call p % moveLocal(distance, n)
        exit STLoop
      end if
    end do STLoop

    call currentCell % kill()

  end subroutine surfaceTracking

  !!
  !! Apply boundary conditions when using surface tracking
  !! This routine may need checking for cases in which the boundary
  !! is crossed repeatedly
  !!
  subroutine applyBC(self, p, currentCell)
    class(transportOperatorST), intent(in) :: self
    class(particle), intent(inout)       :: p
    class(cell_ptr), intent(inout)       :: currentCell
    type(surface_ptr)                    :: currentSurface

    ! Return to global coordinates - this may be a superfluous call!!
    call p % resetNesting()

    ! Identify which surface the particle crossed at the highest geometry level
    currentSurface = currentCell % whichSurface(p%rGlobal(), -p%dirGlobal())

    ! Check the boundary conditions on the surface
    ! If vacuum, kill the particle
    if (currentSurface % isVacuum()) then
      p % isDead = .TRUE.
      call currentSurface % kill()

    ! If reflective or periodic, must ensure that the surface is a plane!
    else if (currentSurface % isReflective()) then

      ! Move particle back to surface which it crossed at the highest geometry level
      call p % moveGlobal(-NUDGE)

      ! Reflect the particle's direction of travel at the highest geometry level
      call currentSurface % reflect(p%coords%lvl(1)%r, p%coords%lvl(1)%dir)

      ! Nudge the particle to avoid problems at corners
      call p % moveGlobal(NUDGE)

      ! Identify which cell the particle now occupies
      currentCell = self % geom % whichCell(p%coords)

      ! The particle should no longer be sat on a surface
      call currentSurface % kill()

    else if (currentSurface % isPeriodic()) then

      ! Apply the periodic translation associated with the surface
      call p % teleport(p%rGlobal() + currentSurface % periodicTranslation() &
                        + NUDGE * p%dirGlobal())

      ! Identify which cell the particle now occupies
      currentCell = self % geom % whichCell(p%coords)

      ! The particle should no longer be sat on a surface
      call currentSurface % kill()

    else
      call fatalError('applyBC, transportOperatorST',&
      'Could not identify correct boundary conditions')
    end if

  end subroutine applyBC

end module transportOperatorST_class
