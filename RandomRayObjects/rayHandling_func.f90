module rayHandling_func

  use numPrecision
  use universalVariables
  use constantsRR
  use genericProcedures,              only : fatalError, numToChar, rotateVector
  use dictionary_class,               only : dictionary
  use rng_class,                      only : RNG

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : distCache
  use geometryStd_class,              only : geometryStd

  ! Random ray modules
  use arraysRR_class,                 only : arraysRR
  use dataRR_class,                   only : dataRR
  use mathsRR_func,                   only : F1

  ! Random ray - or a standard particle
  use particle_class,                 only : ray => particle

  implicit none
  private

  !!
  !! Set of functions and subroutines to handle everything to do with rays
  !! in random ray
  !!
  !! TODO: add uncollided sweep and volume tracing
  !!
  public  :: initialiseRay
  public  :: transportSweep
  private :: moveRay
  private :: checkRayLength
  private :: transportSweepFlatIso
  private :: transportSweepLinIso
  private :: transportSweepLIFA
  private :: transportSweepFlatAni

contains

  subroutine initialiseRay(r, arrays)
    type(ray), intent(inout)                   :: r
    class(arraysRR), pointer, intent(in)       :: arrays
    class(geometryStd), pointer                :: geom
    real(defReal)                              :: mu, phi
    real(defReal), dimension(6)                :: b
    real(defReal), dimension(3)                :: lb, ub, u, rand3, x
    integer(shortInt)                          :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (rayHandling_func.f90)'

    geom => arrays % getGeomPointer()
    
    i = 0
    mu = TWO * r % pRNG % get() - ONE
    phi = TWO_PI * r % pRNG % get()
    u = rotateVector([ONE, ZERO, ZERO], mu, phi)

    b = geom % bounds()
    lb = b(1:3)
    ub = b(4:6)

    rejection : do
      rand3(1) = r % pRNG % get()
      rand3(2) = r % pRNG % get()
      rand3(3) = r % pRNG % get()
      x = lb + (ub - lb) * rand3

      ! Exit if point is inside the geometry
      call geom % whatIsAt(matIdx, cIdx, x, u)
      if (matIdx /= OUTSIDE_MAT) exit rejection

      i = i + 1
      if (i > 5000) then
        call fatalError(Here, 'Infinite loop when searching for ray start in the geometry.')
      end if
    end do rejection

    ! Place in the geometry & process the ray
    call r % build(x, u, 1, ONE)
    call geom % placeCoord(r % coords)

    if (.not. arrays % found(cIdx)) call arrays % newFound(cIdx, x)

  end subroutine initialiseRay
  
  !!
  !! Move ray across a cell, into the next  
  !! Use distance caching or standard ray tracing
  !! Distance caching seems a little bit more unstable
  !! due to FP error accumulation, but is faster.
  !! This can be fixed by resetting the cache after X number
  !! of distance calculations.
  !!
  subroutine moveRay(r, doCache, ints, geom, length, event, cache, hitVacuum)
    type(ray), intent(inout)                 :: r
    logical(defBool), intent(in)             :: doCache
    integer(longInt), intent(inout)          :: ints
    class(geometryStd), pointer, intent(in)  :: geom 
    real(defReal), intent(inout)             :: length
    integer(shortInt), intent(out)           :: event  
    type(distCache), intent(inout)           :: cache
    logical(defBool), intent(out)            :: hitVacuum

    if (doCache) then
      if (mod(ints,20_longInt) == 0)  cache % lvl = 0
      call geom % moveRay_withCache(r % coords, length, event, cache, hitVacuum)
    else
      call geom % moveRay_noCache(r % coords, length, event, hitVacuum)
    end if
    ints = ints + 1

  end subroutine moveRay
  
  !!
  !! Set maximum flight distance and ensure ray is active
  !!
  subroutine checkRayLength(totalLength, dead, termination, activeRay, length)
    real(defReal), intent(in)       :: totalLength
    real(defReal), intent(in)       :: dead
    real(defReal), intent(in)       :: termination
    logical(defBool), intent(inout) :: activeRay
    real(defReal), intent(out)      :: length
      
    if (totalLength >= dead) then
      length = termination - totalLength 
      activeRay = .true.
    else
      length = dead - totalLength
    end if

  end subroutine checkRayLength
  
  !!
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume.
  !! Records the number of integrations/ray movements.
  !!
  subroutine transportSweep(r, ints, nG, doCache, dead, termination, arrays)
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt), intent(in)                         :: nG
    logical(defBool), intent(in)                          :: doCache
    real(defReal), intent(in)                             :: dead
    real(defReal), intent(in)                             :: termination
    class(arraysRR), pointer, intent(inout)               :: arrays
    integer(shortInt)                                     :: simType
    character(100), parameter :: Here = 'transportSweep (rayHandling_func.f90)'

    simType = arrays % getSimulationType()

    select case(simType)
      case (flatIso)
        call transportSweepFlatIso(r, ints, nG, doCache, dead, termination, arrays)
      case (linearIso)
        call transportSweepLinIso(r, ints, nG, doCache, dead, termination, arrays)
      case (flatAni)
        call transportSweepFlatAni(r, ints, nG, doCache, dead, termination, arrays)
      case (linearAni)
        call transportSweepLIFA(r, ints, nG, doCache, dead, termination, arrays)
      case default
        call fatalError(Here,'Unsupported simulation type')
    end select

  end subroutine transportSweep

  !!
  !! Transport sweep for flat isotropic sources
  !!
  subroutine transportSweepFlatIso(r, ints, nG, doCache, dead, termination, arrays)
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt), intent(in)                         :: nG
    logical(defBool), intent(in)                          :: doCache
    real(defReal), intent(in)                             :: dead
    real(defReal), intent(in)                             :: termination
    class(arraysRR), pointer, intent(in)                  :: arrays
    class(dataRR), pointer                                :: XSData
    class(geometryStd), pointer                           :: geom
    integer(shortInt)                                     :: matIdx, g, cIdx, event, matIdx0
    real(defReal)                                         :: totalLength, length
    logical(defBool)                                      :: activeRay, hitVacuum
    type(distCache)                                       :: cache
    real(defFlt)                                          :: lenFlt
    real(defFlt), dimension(nG)                           :: attenuate, delta, fluxVec, tau
    real(defFlt), pointer, dimension(:)                   :: scalarVec, sourceVec, totVec
    
    XSData => arrays % getDataPointer()
    geom => arrays % getGeomPointer()

    ! Set initial angular flux to angle average of cell source
    cIdx = r % coords % uniqueID
    matIdx  = r % coords % matIdx
    call XSData % getTotalPointer(matIdx, totVec)
    
    ! Catch for regions with voids
    ! Assumes these are defined as 'void'
    ! TODO: Use a more robust criterion, as for branching later
    if (matIdx <= XSData % getNMat()) then
      do g = 1, nG
        fluxVec(g) = arrays % getSource(cIdx,g) / totVec(g)
      end do
    else
      do g = 1, nG
        fluxVec(g) = arrays % getPrevFlux(cIdx,g)
      end do
    end if

    ints = 0
    matIdx0 = matIdx
    totalLength = ZERO
    activeRay = .false.
    do while (totalLength < termination)

      ! Get material and cell the ray is moving through
      matIdx = r % coords % matIdx
      cIdx   = r % coords % uniqueID
      if (matIdx0 /= matIdx) then
        matIdx0 = matIdx
        
        ! Cache total cross section
        call XSData % getTotalPointer(matIdx, totVec)
      end if

      ! Set maximum flight distance and ensure ray is active
      call checkRayLength(totalLength, dead, termination, activeRay, length)

      ! Move ray
      call moveRay(r, doCache, ints, geom, length, event, cache, hitVacuum)
      totalLength = totalLength + length
      
      ! Set new cell's position. Use half distance across cell
      ! to try and avoid FP error
      if (.not. arrays % found(cIdx)) then
        call arrays % newFound(cIdx, r % rGlobal() - length * HALF * r % dirGlobal())
      end if

      lenFlt = real(length,defFlt)
      call arrays % getSourcePointer(cIdx, sourceVec)

      ! Branch for voids etc
      ! TODO: Should use a better branching criterion. Maybe create it in data?
      ! Standard route
      if (matIdx <= XSData % getNMat()) then

        !$omp simd
        do g = 1, nG
          tau(g) = totVec(g) * lenFlt
          attenuate(g) = lenFlt * F1(tau(g))
          delta(g) = (totVec(g) * fluxVec(g) - sourceVec(g)) * attenuate(g)
          fluxVec(g) = fluxVec(g) - delta(g)
        end do

        ! Accumulate to scalar flux
        if (activeRay) then
      
          call arrays % setLock(cIdx)
            call arrays % getFluxPointer(cIdx, scalarVec)
            !$omp simd
            do g = 1, nG
              scalarVec(g) = scalarVec(g) + delta(g) 
            end do
            call arrays % incrementVolume(cIdx, length)
          call arrays % unsetLock(cIdx)
          if (arrays % hasHit(cIdx) == 0) call arrays % hitCell(cIdx)
      
        end if

      ! Route for void materials
      else
        
        ! Accumulate to scalar flux
        if (activeRay) then
      
          call arrays % setLock(cIdx)
            call arrays % getFluxPointer(cIdx, scalarVec)
            !$omp simd
            do g = 1, nG
              scalarVec(g) = scalarVec(g) + fluxVec(g) * lenFlt 
            end do
            call arrays % incrementVolume(cIdx, length)
            call arrays % incrementLengthSquared(cIdx, length)
          call arrays % unsetLock(cIdx)
          if (arrays % hasHit(cIdx) == 0) call arrays % hitCell(cIdx)
      
        end if

        !$omp simd
        do g = 1, nG
          fluxVec(g) = fluxVec(g) + sourceVec(g) * lenFlt
        end do

      end if

      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, nG
          fluxVec(g) = 0.0_defFlt
        end do
      end if

    end do

  end subroutine transportSweepFlatIso
  
  !!
  !! Transport sweep for flat isotropic sources
  !!
  subroutine transportSweepLinIso(r, ints, nG, doCache, dead, termination, arrays)
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt), intent(in)                         :: nG
    logical(defBool), intent(in)                          :: doCache
    real(defReal), intent(in)                             :: dead
    real(defReal), intent(in)                             :: termination
    class(arraysRR), pointer, intent(inout)               :: arrays

  end subroutine transportSweepLinIso
  
  !!
  !! Transport sweep for LIFA sources
  !!
  subroutine transportSweepLIFA(r, ints, nG, doCache, dead, termination, arrays)
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt), intent(in)                         :: nG
    logical(defBool), intent(in)                          :: doCache
    real(defReal), intent(in)                             :: dead
    real(defReal), intent(in)                             :: termination
    class(arraysRR), pointer, intent(inout)               :: arrays

  end subroutine transportSweepLIFA
  
  !!
  !! Transport sweep for flat aniisotropic sources
  !!
  subroutine transportSweepFlatAni(r, ints, nG, doCache, dead, termination, arrays)
    type(ray), intent(inout)                              :: r
    integer(longInt), intent(out)                         :: ints
    integer(shortInt), intent(in)                         :: nG
    logical(defBool), intent(in)                          :: doCache
    real(defReal), intent(in)                             :: dead
    real(defReal), intent(in)                             :: termination
    class(arraysRR), pointer, intent(inout)               :: arrays

  end subroutine transportSweepFlatAni
  

end module rayHandling_func
