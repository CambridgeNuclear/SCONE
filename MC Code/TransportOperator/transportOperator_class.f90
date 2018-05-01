!
! Transport operator - determines how a neutron moves from point r' to r
! given initial energy and direction
!
!

module transportOperator_class
  use numPrecision
  use genericProcedures
  use universalVariables

  use coord_class
  use particle_class
  use surface_class
  use geometry_class
  use cell_class
  use universe_class
  use rng_class
  use lattice_class

  !use Nuclear_Data_MG_class           !Re-instate later!!!
  !use Nuclear_Data_CE_class
  use geometry_class

  implicit none
  private


  type, public :: transportOperator
    private
    class(rng), pointer      :: random => null()         ! RNG - should this be associated to particle?
    !class(Nuclear_data_MG), pointer :: MGData => null()  ! multi-group data
    !class(Nuclear_data_CE), pointer :: CEData => null()  ! continuous energy data
    class(geometry), pointer :: geom => null()           ! references the geometry for cell searching
    logical                  :: isDT = .true.            ! perform delta tracking?
  contains
    procedure :: initMG
    procedure :: performTransport
    procedure :: getSigmaT
    procedure :: getMajorant
    procedure :: deltaTracking
    procedure :: surfaceTracking
    procedure :: applyBCsDT
    procedure :: applyBCsST
    procedure :: walk
  end type transportOperator

contains

  !!
  !! Initialise for multi-group data
  !!
  subroutine initMG(self, random, geom) !return nuclearData at some point!
    class(transportOperator), intent(inout) :: self
    !class(Nuclear_data_MG), target :: nuclearData
    class(rng), target                      :: random
    class(geometry), target                 :: geom

    self%random => random
    !self%MGData => nuclearData
    self%geom => geom

  end subroutine initMG

  !!
  !! Find the total cross-section from the nuclear data provided, continuous or multi-group
  !!
  function getSigmaT(self,p)result(sigmaT)
    class(transportOperator), intent(in) :: self
    class(particle), intent(in)          :: p
    real(defReal)                        :: sigmaT

    if (p % isMG) then
      !sigmaT = self % MGData % giveTotalXS(p)
      sigmaT = 0.6
      return
    else
      !sigmaT = self % CEData % giveTotalXS(p)
      sigmaT = 0.6
      return
    end if

  end function getSigmaT

  !!
  !! Find the majorant cross-section from the nuclear data provided, continuous or multi-group
  !!
  function getMajorant(self,p)result(majorant)
    class(transportOperator), intent(in) :: self
    class(particle), intent(in)          :: p
    real(defReal)                        :: majorant

    if (p % isMG) then
      !majorant = self % MGData % giveMajorantXS(p)
      majorant = 0.8
      return
    else
      !majorant = self % CEData % giveMajorantXS(p)
      majorant = 0.8
      return
    end if

  end function getMajorant

  !!
  !! Transport particle from r' to r given current direction and energy
  !!
  subroutine performTransport(self,p)
    class(transportOperator), intent(in) :: self
    class(particle), intent(inout)       :: p

    if (self % isDT) then
      call self % deltaTracking(p)
    else
      call self % surfaceTracking(p)
    end if

  end subroutine performTransport

  !!
  !! Performs delta tracking
  !!
  subroutine deltaTracking(self,p)
    class(transportOperator), intent(in) :: self
    class(particle), intent(inout)       :: p
    real(defReal)                        :: majorant, sigmaT, distance
    type(cell_ptr)                       :: currentCell

    majorant = self % getMajorant(p)
    DTLoop:do
      distance = -log(self%random%get())/majorant

      ! Move particle to new location in the global co-ordinate systems
      call p % moveGlobal(distance)

      ! Find the new cell which the particle occupies
      currentCell = self % geom % whichCell(p%coords)

      ! If the particle is outside the geometry, apply boundary conditions
      if (.not. currentCell % insideGeom()) then
        call self % applyBCsDT(p, currentCell)
        ! End the transport step if the particle is killed
        if (p % isDead) then
          exit DTLoop
        end if
      end if

      ! Obtain the local cross-section
      sigmaT = self % getSigmaT(p)

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
  subroutine applyBCsDT(self, p, currentCell)
    class(transportOperator), intent(in) :: self
    class(particle), intent(inout)       :: p
    class(cell_ptr), intent(inout)       :: currentCell
    type(surface_ptr)                    :: currentSurface

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
        call p % coords % resetNesting()
        call currentSurface % reflectiveTransform(p%coords%lvl(1)%r, p%coords%lvl(1)%dir)

        ! Identify which cell the particle now occupies
        currentCell = self % geom % whichCell(p%coords)

      else if (currentSurface % isPeriodic()) then

        ! Return particle to gloabl coordinates and apply the periodic translation associated with the surface
        call p % coords % resetNesting()
        p % coords % lvl(1) % r = p % rGlobal() + currentSurface % periodicTranslation()

        ! Identify which cell the particle now occupies
        currentCell = self % geom % whichCell(p%coords)

      else
        call fatalError('applyBCsDT','Could not identify correct boundary conditions')
      end if

    end do

    ! The particle should no longer be sat on a surface
    call currentSurface % kill()

  end subroutine applyBCsDT


  !!
  !! Performs surface tracking
  !!
  subroutine surfaceTracking(self,p)
    class(transportOperator), intent(in) :: self
    class(particle), intent(inout)       :: p
    real(defReal)                        :: sigmaT, distance, boundaryDistance, &
                                            testDistance, latDistance, lDist
    integer(shortInt)                    :: i, n, nMin, latIdx, iLat
    type(cell_ptr)                       :: c, currentCell
    type(lattice_ptr)                    :: lat
    logical(defBool)                     :: moveUp

    STLoop: do
      ! Obtain the local cross-section
      sigmaT = self % getSigmaT(p)

      ! Calculate boundary distance: descend the different co-ordinate levels starting from the highest
      ! Ensures bounds of parent cells are not exceeded
      n = p % coords % nesting
      nMin = 1
      boundaryDistance = INFINITY
      latDistance = INFINITY
      do i = 1,n
        c = self % geom % cells(p % coords % lvl(i) % cellIdx)
        testDistance = c % getDistance(p%rLocal(i), p%dirLocal(i))

        ! Check if the particle is in a lattice cell
        latIdx = p % coords % lvl(i) % latIdx
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

      ! Sample particle flight distance
      distance = -log(self%random%get())/sigmaT

      ! The particle escapes the cell and moves to the next
      !if (abs(boundaryDistance - distance)/boundaryDistance >= surface_tol) then
      if (boundaryDistance <= distance) then

        ! Move particle to the surface with a small nudge to move across the boundary
        call p % moveLocal(boundaryDistance + NUDGE, n)

        ! Find the new base cell which the particle occupies
        currentCell = self % geom % whichCell(p%coords,n)

        ! If the particle is outside the geometry, apply boundary conditions
        do while (.not. currentCell % insideGeom())
          call self % applyBCsST(p, currentCell)
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
  subroutine applyBCsST(self, p, currentCell)
    class(transportOperator), intent(in) :: self
    class(particle), intent(inout)       :: p
    class(cell_ptr), intent(inout)       :: currentCell
    type(surface_ptr)                    :: currentSurface

    ! Return to global coordinates - this may be a superfluous call!!
    call p % coords % resetNesting()

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
      p%coords%lvl(1)%r = p%rGlobal() - NUDGE * p%dirGlobal()

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
      p%coords%lvl(1)%r = p%rGlobal() + currentSurface % periodicTranslation() &
                          + NUDGE * p%dirGlobal()

      ! Identify which cell the particle now occupies
      currentCell = self % geom % whichCell(p%coords)

      ! The particle should no longer be sat on a surface
      call currentSurface % kill()

    else
      call fatalError('applyBCsST','Could not identify correct boundary conditions')
    end if

  end subroutine applyBCsST

  !!
  !! Walks the particle for a number of steps or until the particle is killed
  !! Does this by surface tracking without accounting for material properties
  !!
  subroutine walk(self,p,steps)
    class(transportOperator), intent(in) :: self
    class(particle), intent(inout)       :: p
    integer(shortInt), intent(in)        :: steps
    real(defReal)                        :: boundaryDistance, testDistance, latDistance, lDist
    integer(shortInt)                    :: i, n, nMin, step, latIdx, iLat
    type(cell_ptr)                       :: c, currentCell
    type(lattice_ptr)                    :: lat

    STLoop: do step =1,steps

      ! Calculate boundary distance: descend the different co-ordinate levels starting from the highest
      ! Ensures bounds of parent cells are not exceeded
      n = p % coords % nesting
      nMin = 1
      boundaryDistance = INFINITY
      latDistance = INFINITY
      do i = 1,n
        c = self % geom % cells(p % coords % lvl(i) % cellIdx)
        testDistance = c % getDistance(p%rLocal(i), p%dirLocal(i))

        ! Check if the particle is in a lattice cell
        latIdx = p % coords % lvl(i) % latIdx
        if (latIdx > 0) then
          lat = self % geom % lattices(latIdx)
          lDist = lat % getDistance(p%rLocal(i), p%dirLocal(i))
          if (latDistance >= lDist) then
            latDistance = lDist
            iLat = i - 1
          end if
        end if

        if (abs(boundaryDistance - testDistance)/boundaryDistance >= surface_tol) then
          boundaryDistance = testDistance
          nMin = i
          ! Must move up a level if crossing a lattice boundary
          if (boundaryDistance >= latDistance) then
            nMin = iLat
            boundaryDistance = latDistance
          end if
        end if
      end do
      call c % kill()
      call lat % kill()
      n = nMin

      ! Move particle to the surface with a small nudge to move across the boundary
      call p % moveLocal(boundaryDistance + NUDGE, n)

      ! Find the new base cell which the particle occupies
      currentCell = self % geom % whichCell(p%coords, n)

      ! If the particle is outside the geometry, apply boundary conditions
      do while (.not. currentCell % insideGeom())
        call self % applyBCsST(p, currentCell)
        ! End transport if the particle is killed
        if (p % isDead) then
          print *,'DEAD'
          exit STLoop
        end if
      end do

    end do STLoop

    call currentCell % kill()

  end subroutine walk

end module transportOperator_class
