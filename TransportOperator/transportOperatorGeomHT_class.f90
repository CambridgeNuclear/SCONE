!!
!! Hybrid transport operator using two geometries to switch between delta tracking and surface tracking
!!
module transportOperatorGeomHT_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle, P_PHOTON, P_MATERIAL
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use RNG_class,                  only : RNG

  ! Superclass
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,             only : geometry
  use geometryGrid_class,         only : geometryGrid
  use geometryReg_mod,            only : gr_geomPtr  => geomPtr, gr_addGeom => addGeom, &
                                         gr_geomIdx  => geomIdx
  use coord_class,                only : coordList

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase
  use mgIMCDatabase_inter,        only : mgIMCDatabase, mgIMCDatabase_CptrCast
  use nuclearDataReg_mod,         only : ndReg_get   => get

  use virtualMat_class,           only : virtualMat
  use simulationTime_class

  implicit none
  private

  !!
  !! Transport operator that moves a particle with using hybrid tracking, up to a time boundary
  !!
  type, public, extends(transportOperator)   :: transportOperatorGeomHT
    real(defReal)                            :: deltaT
    real(defReal)                            :: cutoff
    integer(shortInt)                        :: method
    class(geometry), pointer                 :: upperGeom
    integer(shortInt)                        :: upperGeomIdx
    integer(shortInt)                        :: thisTimeStep
    class(virtualMat), dimension(:), allocatable :: virtualMats
  contains
    procedure          :: transit => timeTracking
    procedure          :: init
    procedure, private :: surfaceTracking
    procedure, private :: deltaTracking
    procedure, private :: getMajInv
  end type transportOperatorGeomHT

contains

  subroutine timeTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorGeomHT), intent(inout) :: self
    class(particle), intent(inout)                :: p
    type(tallyAdmin), intent(inout)               :: tally
    class(particleDungeon), intent(inout)         :: thisCycle
    class(particleDungeon), intent(inout)         :: nextCycle
    integer(shortInt)                             :: i
    character(100), parameter :: Here = 'timeTracking (transportOperatorGeomHT_class.f90)' 

    ! Material transform not included in this class to avoid complicating any further. Can be
    ! essentially copied and pasted from transportOperatorTime_class if desired to use ISMC here
    if (p % type == P_MATERIAL) call fatalError(Here, 'No support for ISMC in this transOp')

    ! Update majorants if required - this would be better done at the end of time step in PP
    ! to avoid check for each particle but I wanted to keep this class self-contained
    if (self % thisTimeStep /= thisStep()) then
      self % thisTimeStep = thisStep()
      do i = 1, size(self % virtualMats)
        call self % virtualMats (i) % updateMajorant()
      end do
    end if

    call self % deltaTracking(p)

    ! Check for particle leakage
    if (p % matIdx() == OUTSIDE_FILL) then
      p % fate = LEAK_FATE
      p % isDead = .true.
    end if

    call tally % reportTrans(p)

  end subroutine timeTracking

  !!
  !! Perform delta tracking
  !!
  !! Note that the method used of dealing with upper geometry cell changes does lead to potential
  !! inaccuracy if a particle can change upperGeom cells and then return to the original one in the
  !! same move, for example if upperGeom is curved or if there is a reflection. For now this class
  !! requires upperGeom to be geometryGrid_class (no curves) and this sort of reflection should
  !! never happen for any vaguely sensible geometry/material choices. One could instead explicitly
  !! calculate distance to a new upperGeom cell and compare to dColl and dTime, at the cost of
  !! efficiency.
  !!
  subroutine deltaTracking(self, p)
    class(transportOperatorGeomHT), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(coordList), allocatable                 :: coords
    real(defReal)                                 :: dTime, dColl, sigmaT, majorant_inv, dist, ratio
    integer(shortInt)                             :: virtualMatIdx, testMat, uniqueID, event, i
    character(100), parameter :: Here = 'deltaTracking (transportOperatorGeomHT_class.f90)'

    ! Get majorant
    call self % getMajInv(p, majorant_inv, virtualMatIdx)

    ! Get initial local opacity
    sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

    DTLoop:do

      ! Switch to surface tracking if delta tracking is unsuitable
      ratio = sigmaT*majorant_inv
      if (ratio > ONE) call fatalError(Here, 'Local opacity greater than majorant')
      if (ratio < self % cutoff .or. majorant_inv == ZERO) then
        call self % surfaceTracking(p)
        return
      end if

      ! Find distance to time boundary
      dTime = lightSpeed * (p % timeMax - p % time)

      ! Sample distance to collision
      dColl = -log(p % pRNG % get()) * majorant_inv

      ! Move copy of coords by minimum distance without considering surface crossings
      coords = p % coords
      dist = min(dColl, dTime)
      call self % geom % teleport(coords, dist)

      ! Check for particle leakage
      if (coords % matIdx == OUTSIDE_FILL) then
        p % coords = coords
        return
      end if

      ! Check for change of upper geometry
      call self % upperGeom % whatIsAt(testMat, uniqueID, coords % lvl(1) % r, coords % lvl(1) % dir)
      if (testMat /= virtualMatIdx) then
        ! Move would take particle to a new cell
        call self % upperGeom % move(p % coords, dist, event)
        ! Get new majorant (particle already placed in upper geometry)
        virtualMatIdx = p % matIdx()
        majorant_inv = self % virtualMats(virtualMatIdx) % majorant_inv
        !call self % getMajInv(p, majorant_inv, virtualMatIdx)
        ! Update particle time and place back in lower geometry
        p % time = p % time + dist / lightSpeed
        call self % geom % placeCoord(p % coords)
        sigmaT = self % xsData % getTransMatXS(p, p % matIdx())
        ! Return to start of loop
        cycle DTLoop
      end if

      ! Move accepted, move p
      p % coords = coords

      ! Update particle time
      p % time = p % time + dist / lightSpeed

      ! Act based on distance moved
      if (dist == dTime) then
        ! Update particle fate and exit
        p % fate = AGED_FATE
        p % time = p % timeMax
        exit DTLoop

      else if (dist == dColl) then! Dist == dColl
        ! Get local opacity and check for real or virtual collision
        sigmaT = self % xsData % getTransMatXS(p, p % matIdx())
        if (p % pRNG % get() < sigmaT * majorant_inv) exit DTLoop

      else
        call fatalError(Here, 'aaa')

      end if

    end do DTLoop

  end subroutine deltaTracking

  !!
  !! Perform standard surface tracking using only the lower geometry.
  !! Once in this loop, delta tracking is not used again.
  !!
  subroutine surfaceTracking(self, p)
    class(transportOperatorGeomHT), intent(inout) :: self
    class(particle), intent(inout)                :: p
    real(defReal)                                 :: dTime, dColl, dist, sigmaT
    integer(shortInt)                             :: event
    character(100), parameter :: Here = 'surfaceTracking (transportOperatorGeomHT_class.f90)'

    STLoop:do

      ! Find distance to time boundary
      dTime = lightSpeed * (p % timeMax - p % time)

      ! Sample distance to collision
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())
      dColl = -log( p % pRNG % get() ) / sigmaT

      ! Ensure particle does not remain exactly on a boundary if dColl is close to 0
      if (event == CROSS_EV .and. dColl < SURF_TOL) then
        dColl = SURF_TOL
      end if

      ! Choose minimum distance
      dist = min(dTime, dColl)

      ! Move through geometry using minimum distance
      call self % geom % move(p % coords, dist, event)

      ! Check for particle leakage
      if (p % matIdx() == OUTSIDE_FILL) return

      ! Increase time based on distance moved
      p % time = p % time + dist / lightSpeed

      ! Check result of transport
      if (dist == dTime) then
        ! Time boundary
        p % fate = AGED_FATE
        p % time = p % timeMax
        exit STLoop

      else if (dist == dColl) then
        ! Collision
        exit STLoop

      end if

      if (event == COLL_EV) call fatalError(Here, 'Move outcome should be CROSS_EV or BOUNDARY_EV')

    end do STLoop

    if (event /= COLL_EV) call fatalError(Here, 'Move outcome should be COLL_EV')

  end subroutine surfaceTracking


  !!
  !! Return the inverse majorant opacity
  !! For DT or HT this will be constant, for GT this will be dependent on position
  !!
  !! Args:
  !!   p [in]  -> particle
  !!
  !! Result:
  !!   maj_inv -> 1 / majorant opacity
  !!
  subroutine getMajInv(self, p, majorant_inv, virtualMatIdx)
    class(transportOperatorGeomHT), intent(in) :: self
    class(particle), intent(in)                :: p
    real(defReal), intent(out)                 :: majorant_inv
    integer(shortInt), intent(out)             :: virtualMatIdx
    real(defReal), dimension(3)                :: r, dir
    integer(shortInt)                          :: uniqueID

    ! Get index of virtual material
    r = p % coords % lvl(1) % r
    dir = p % coords % lvl(1) % dir
    call self % upperGeom % whatIsAt(virtualMatIdx, uniqueID, r, dir)

    ! Get 1/majorant
    majorant_inv = self % virtualMats(virtualMatIdx) % majorant_inv

  end subroutine getMajInv

  !!
  !! Provide transport operator with delta tracking/surface tracking cutoff
  !!
  !! Cutoff of 1 gives exclusively delta tracking, cutoff of 0 gives exclusively surface tracking
  !!
  subroutine init(self, dict)
    class(transportOperatorGeomHT), intent(inout)    :: self
    class(dictionary), intent(in)                    :: dict
    character(nameLen)                               :: geomName
    class(dictionary),pointer                        :: tempdict
    class(geometry), pointer                         :: upperGeom
    integer(shortInt), dimension(:), allocatable     :: dimensions, searchN
    integer(shortInt), dimension(3)                  :: searchN3
    integer(shortInt)                                :: N, i, j, k
    integer(shortInt)                                :: realMatIdx, virtualMatIdx, uniqueID
    real(defReal), dimension(6)                      :: bounds
    real(defReal), dimension(3)                      :: corner, extent, searchRes, r
    character(100), parameter :: Here = "init (transportOperatorGeomHT_class.f90)"

    ! Initialise superclass
    call init_super(self, dict)

    ! Get cutoff value
    call dict % getOrDefault(self % cutoff, 'cutoff', 0.3_defReal)
    ! Flip to be consistent with transportOperatorHT_class
    self % cutoff = ONE - self % cutoff

    ! Build upper level geometry
    geomName = 'surfaceGeom'
    tempDict => dict % getDictPtr('geometry')
    call gr_addGeom(geomName, tempDict)
    self % upperGeomIdx =  gr_geomIdx(geomName)
    upperGeom           => gr_geomPtr(self % upperGeomIdx)
    self % upperGeom    => upperGeom

    ! Provide access to lower (standard) geometry
    ! TODO: This assumes that there is only 1 other defined geometry
    self % geom => gr_geomPtr(1)
    self % xsData => ndReg_get(P_PHOTON_MG)

    ! For now limited to grid geometry
    select type(upperGeom)
      class is(geometryGrid)
        ! Get some basic geometry info
        corner = upperGeom % corner
      class default
        call fatalError(Here, 'Geometry class should be geometryGrid')
    end select

    ! Initialise a virtual material object for each cell
    call tempDict % get(dimensions,'dimensions')
    N = dimensions(1) * dimensions(2) * dimensions(3)
    allocate(self % virtualMats(N))
    do i = 1, N
      call self % virtualMats(i) % init(self % xsData)
    end do

    ! Get resolution to search through grid
    call dict % get(searchN, 'searchN')
    if (size(searchN) == 3) then
      searchN3 = searchN
    else if (size(searchN) == 1) then
      do i = 1, 3
        searchN3(i) = searchN(1)
      end do
    else
      call fatalError(Here, 'searchN must be of size 1 or 3')
    end if
    if (any(searchN3 < 1)) call fatalError(Here, 'Invalid searchN')

    ! Search grid to assign real materials to virtual materials
    bounds = upperGeom % bounds()
    do i = 1, 3
      extent(i) = bounds(i+3) - bounds(i)
    end do
    searchRes = extent / (searchN3 + 1)
    do i = 1, searchN3(1)
      do j = 1, searchN3(2)
        do k = 1, searchN3(3)
          ! Find matIdx at search location
          r = corner + [i, j, k] * searchRes
          call self % geom % whatIsAt(realMatIdx, uniqueID, r)
          call upperGeom % whatIsAt(virtualMatIdx, uniqueID, r)
          call self % virtualMats(virtualMatIdx) % addRealMat(realMatIdx)
        end do
      end do
    end do

    do i = 1, size(self % virtualMats)
      call self % virtualMats (i) % updateMajorant()
    end do

  end subroutine init

end module transportOperatorGeomHT_class
