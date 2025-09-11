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
  use mathsRR_func,                   only : expF1, expF1Tau, expG, expG2

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
  private :: transportSweepLinearIso
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

    if (.not. arrays % wasFound(cIdx)) call arrays % newFound(cIdx, x)

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
    type(ray), intent(inout)                :: r
    integer(longInt), intent(out)           :: ints
    integer(shortInt), intent(in)           :: nG
    logical(defBool), intent(in)            :: doCache
    real(defReal), intent(in)               :: dead
    real(defReal), intent(in)               :: termination
    class(arraysRR), pointer, intent(inout) :: arrays
    integer(shortInt)                       :: simType
    character(100), parameter :: Here = 'transportSweep (rayHandling_func.f90)'

    simType = arrays % getSimulationType()

    select case(simType)
      case (flatIso)
        call transportSweepFlatIso(r, ints, nG, doCache, dead, termination, arrays)
      case (linearIso)
        call transportSweepLinearIso(r, ints, nG, doCache, dead, termination, arrays)
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
    type(ray), intent(inout)             :: r
    integer(longInt), intent(out)        :: ints
    integer(shortInt), intent(in)        :: nG
    logical(defBool), intent(in)         :: doCache
    real(defReal), intent(in)            :: dead
    real(defReal), intent(in)            :: termination
    class(arraysRR), pointer, intent(in) :: arrays
    class(dataRR), pointer               :: XSData
    class(geometryStd), pointer          :: geom
    integer(shortInt)                    :: matIdx, g, cIdx, event, matIdx0
    real(defReal)                        :: totalLength, length
    logical(defBool)                     :: activeRay, hitVacuum
    type(distCache)                      :: cache
    real(defFlt)                         :: lenFlt, len_2
    real(defFlt), dimension(nG)          :: attenuate, delta, angular, tau, inc
    real(defFlt), pointer, dimension(:)  :: source, total
    real(defReal), pointer, dimension(:) :: scalar
    
    XSData => arrays % getDataPointer()
    geom => arrays % getGeomPointer()

    ! Set initial angular flux to angle average of cell source
    cIdx = r % coords % uniqueID
    matIdx  = r % coords % matIdx
    call XSData % getTotalPointer(matIdx, total)
    
    ! Catch for regions with voids
    ! Assumes these are defined as 'void'
    ! TODO: Use a more robust criterion, as for branching later
    if (matIdx <= XSData % getNMat()) then
      do g = 1, nG
        if (total(g) > 1.0E-6_defFlt) then
          angular(g) = arrays % getSource(cIdx,g) / total(g)
        else
          angular(g) = real(arrays % getPrevFlux(cIdx, g), defFlt)
        end if
      end do
    else
      do g = 1, nG
        angular(g) = real(arrays % getPrevFlux(cIdx, g), defFlt)
        !angular(g) = 0.0_defFlt
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
        call XSData % getTotalPointer(matIdx, total)
      end if

      ! Set maximum flight distance and ensure ray is active
      call checkRayLength(totalLength, dead, termination, activeRay, length)

      ! Move ray
      call moveRay(r, doCache, ints, geom, length, event, cache, hitVacuum)
      totalLength = totalLength + length
      
      ! Set new cell's position. Use half distance across cell
      ! to try and avoid FP error
      if (.not. arrays % wasFound(cIdx)) then
        call arrays % newFound(cIdx, r % rGlobal() - length * HALF * r % dirGlobal())
      end if
      
      lenFlt = real(length,defFlt)
      call arrays % getSourcePointer(cIdx, source)

      ! Branch for voids etc
      ! TODO: Should use a better branching criterion. Maybe create it in data?
      ! Standard route
      if (matIdx <= XSData % getNMat()) then
      
        !$omp simd
        do g = 1, nG
          !tau(g) = max(total(g) * lenFlt, 1.0e-6_defFlt) ! Not sure whether this does much for stability
          tau(g) = total(g) * lenFlt
        end do

        !$omp simd
        do g = 1, nG
          attenuate(g) = lenFlt * expF1(tau(g))
          delta(g) = (total(g) * angular(g) - source(g)) * attenuate(g)
          angular(g) = angular(g) - delta(g)
        end do
 
        ! Accumulate to scalar flux
        if (activeRay) then
      
          call arrays % setLock(cIdx)
            call arrays % getFluxPointer(cIdx, scalar)
            !$omp simd
            do g = 1, nG
              scalar(g) = scalar(g) + delta(g)
            end do
            call arrays % incrementVolume(cIdx, length)
            call arrays % hitCell(cIdx)
          call arrays % unsetLock(cIdx)
      
        end if

      ! Route for void materials
      else
        
        ! Accumulate to scalar flux
        if (activeRay) then
      
          len_2 = lenFlt * one_two
          !$omp simd
          do g = 1, nG
            inc(g) = lenFlt * (angular(g) + source(g) * len_2)
          end do

          call arrays % setLock(cIdx)
            call arrays % getFluxPointer(cIdx, scalar)
            !$omp simd
            do g = 1, nG
              scalar(g) = scalar(g) + inc(g)
            end do
            call arrays % incrementVolume(cIdx, length)
            call arrays % hitCell(cIdx)
          call arrays % unsetLock(cIdx)
      
        end if

        !$omp simd
        do g = 1, nG
          angular(g) = angular(g) + source(g) * lenFlt
        end do

      end if

      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, nG
          angular(g) = 0.0_defFlt
        end do
      end if

    end do

  end subroutine transportSweepFlatIso
  
  !!
  !! Transport sweep for linear isotropic sources
  !!
  subroutine transportSweepLinearIso(r, ints, nG, doCache, dead, termination, arrays)
    type(ray), intent(inout)                :: r
    integer(longInt), intent(out)           :: ints
    integer(shortInt), intent(in)           :: nG
    logical(defBool), intent(in)            :: doCache
    real(defReal), intent(in)               :: dead
    real(defReal), intent(in)               :: termination
    class(arraysRR), pointer, intent(inout) :: arrays
    class(dataRR), pointer                  :: XSData
    class(geometryStd), pointer             :: geom
    integer(shortInt)                       :: matIdx, g, cIdx, event, matIdx0
    real(defReal)                           :: totalLength, length, len2_12
    real(defReal), dimension(nDim)          :: mid, r0, rC, mu0, rNorm
    real(defReal), dimension(matSize)       :: matScore
    logical(defBool)                        :: activeRay, hitVacuum
    type(distCache)                         :: cache
    real(defFlt)                            :: lenFlt, lenFlt2_2, len_2
    real(defFlt), dimension(nDim)           :: muFlt, r0NormFlt, rNormFlt
    real(defFlt), dimension(nG)             :: delta, angular, tau, flatQ, gradQ, &
                                               F1, F2, angular0, G0, G1, G2, H, &
                                               xInc, yInc, zInc, inc
    real(defFlt), pointer, dimension(:)     :: source, total, sourceX, sourceY, sourceZ
    real(defReal), pointer, dimension(:)    :: scalar, scalarX, scalarY, scalarZ
    character(100), parameter :: Here = 'transportSweepLinearIso (rayHandling_func.f90)'
    
    XSData => arrays % getDataPointer()
    geom => arrays % getGeomPointer()

    ! Set initial angular flux to angle average of cell source
    cIdx = r % coords % uniqueID
    matIdx  = r % coords % matIdx
    call XSData % getTotalPointer(matIdx, total)
    
    ! Catch for regions with voids
    ! Assumes these are defined as 'void'
    ! TODO: Use a more robust criterion, as for branching later
    if (matIdx <= XSData % getNMat()) then
      do g = 1, nG
        angular(g) = arrays % getSource(cIdx,g) / total(g)
      end do
    else
      do g = 1, nG
        angular(g) = real(arrays % getPrevFlux(cIdx,g), defFlt)
      end do
    end if
      
    ints = 0
    matIdx0 = matIdx
    totalLength = ZERO
    activeRay = .false.
    do while (totalLength < termination)

      ! Get ray entry position and direction for LS calculations
      r0  = r % rGlobal()
      mu0 = r % dirGlobal()
      muFlt = real(mu0,defFlt)

      ! Get material and cell the ray is moving through
      matIdx  = r % coords % matIdx
      cIdx    = r % coords % uniqueID
      if (matIdx0 /= matIdx) then
        matIdx0 = matIdx
        
        ! Cache total cross section
        call XSData % getTotalPointer(matIdx, total)
      end if

      ! Set maximum flight distance and ensure ray is active
      call checkRayLength(totalLength, dead, termination, activeRay, length)

      ! Move ray
      call moveRay(r, doCache, ints, geom, length, event, cache, hitVacuum)
      totalLength = totalLength + length
      
      ! Calculate the track centre
      rC = r0 + length * HALF * mu0
      
      ! Set new cell's position
      if (.not. arrays % wasFound(cIdx)) call arrays % newFound(cIdx, rC)
      
      ! Compute the track centroid and entry point in local co-ordinates
      ! Convert to floats for speed
      ! If region is rarely visited, use ray's halfway point as centroid
      ! Prevents numerical trouble
      if (arrays % getVolume(cIdx) > ZERO) then
        mid = arrays % getCentroid(cIdx)
        rNorm = rC - mid
        rNormFlt = real(rNorm,defFlt)
        r0NormFlt = real(r0 - mid,defFlt)
      else
        rNorm = ZERO
        rNormFlt = 0.0_defFlt
        r0NormFlt = -real(HALF * mu0 * length,defFlt)
      end if

      call arrays % getSourcePointer(cIdx, source)
      call arrays % getSourceXYZPointers(cIdx, sourceX, sourceY, sourceZ)

      ! Calculate source terms
      !$omp simd aligned(sourceX, sourceY, sourceZ)
      do g = 1, nG
        flatQ(g) = rNormFlt(x) * sourceX(g)
        flatQ(g) = flatQ(g) + rNormFlt(y) * sourceY(g)
        flatQ(g) = flatQ(g) + rNormFlt(z) * sourceZ(g)
        flatQ(g) = flatQ(g) + source(g)

        gradQ(g) = muFlt(x) * sourceX(g)
        gradQ(g) = gradQ(g) + muFlt(y) * sourceY(g)
        gradQ(g) = gradQ(g) + muFlt(z) * sourceZ(g)
      end do

      lenFlt = real(length,defFlt)
      lenFlt2_2 = lenFlt * lenFlt * one_two
      
      ! Branch for voids etc
      ! TODO: Should use a better branching criterion. Maybe create it in data?
      ! Standard route
      if (matIdx <= XSData % getNMat()) then
      
        ! Compute exponentials necessary for angular flux update
        !$omp simd
        do g = 1, nG
          !tau(g) = max(total(g) * lenFlt, 1.0e-8_defFlt) ! This line worsens stability
          tau(g) = total(g) * lenFlt
        end do
      
        !$omp simd
        do g = 1, nG
          G0(g)  = expG(tau(g))
        end do
      
        !$omp simd
        do g = 1, nG
          F1(g)  = 1.0_defFlt - tau(g) * G0(g)
        end do
      
        !$omp simd
        do g = 1, nG
          F2(g)  = 2.0_defFlt * G0(g) - F1(g)
        end do
      
        !$omp simd
        do g = 1, nG
          delta(g) = (tau(g) * angular(g) - lenFlt * flatQ(g)) * F1(g) &
                     - gradQ(g) * F2(g) * lenFlt2_2
        end do
      
        ! Create an intermediate flux variable for use in LS scores
        !$omp simd
        do g = 1, nG
          angular0(g) = angular(g)
        end do
      
        !$omp simd
        do g = 1, nG
          angular(g) = angular(g) - delta(g)
        end do

        ! Accumulate to scalar flux
        if (activeRay) then
      
          ! Precompute geometric info to keep it out of the lock
          len2_12 = length * length * one_twelve
          matScore(xx) = length * (rNorm(x) * rNorm(x) + mu0(x) * mu0(x) * len2_12)
          matScore(xy) = length * (rNorm(x) * rNorm(y) + mu0(x) * mu0(y) * len2_12)
          matScore(xz) = length * (rNorm(x) * rNorm(z) + mu0(x) * mu0(z) * len2_12)
          matScore(yy) = length * (rNorm(y) * rNorm(y) + mu0(y) * mu0(y) * len2_12)
          matScore(yz) = length * (rNorm(y) * rNorm(z) + mu0(y) * mu0(z) * len2_12)
          matScore(zz) = length * (rNorm(z) * rNorm(z) + mu0(z) * mu0(z) * len2_12)
          rC = rC * length
        
          ! Compute necessary exponentials outside of the lock
          ! Follows those in Gunow
        
          !$omp simd
          do g = 1, nG
            H(g) = F1(g) - G0(g)
          end do
        
          !$omp simd
          do g = 1, nG
            G1(g) = one_two - H(g)
          end do
        
          !$omp simd
          do g = 1, nG
            G2(g) = expG2(tau(g)) 
          end do
     
          ! Make some more condensed variables to help vectorisation
          !$omp simd 
          do g = 1, nG
            G1(g) = G1(g) * flatQ(g) * lenFlt
            G2(g) = G2(g) * gradQ(g) * lenFlt2_2 
            H(g)  = H(g) * angular0(g) * tau(g)
            H(g) = (G1(g) + G2(g) + H(g)) * lenFlt
            flatQ(g) = flatQ(g) * lenFlt + delta(g)
          end do
        
          !$omp simd
          do g = 1, nG
            xInc(g) = r0NormFlt(x) * flatQ(g) + muFlt(x) * H(g) 
            yInc(g) = r0NormFlt(y) * flatQ(g) + muFlt(y) * H(g) 
            zInc(g) = r0NormFlt(z) * flatQ(g) + muFlt(z) * H(g) 
          end do

          call arrays % setLock(cIdx)
          
            call arrays % getFluxPointer(cIdx, scalar)
            call arrays % getFluxXYZPointers(cIdx, scalarX, scalarY, scalarZ)

            ! Update flux moments
            !$omp simd aligned(scalar, scalarX, scalarY, scalarZ)
            do g = 1, nG
              scalar(g) = scalar(g) + delta(g) 
              scalarX(g) = scalarX(g) + xInc(g) 
              scalarY(g) = scalarY(g) + yInc(g)
              scalarZ(g) = scalarZ(g) + zInc(g) 
            end do
            
            call arrays % incrementVolume(cIdx, length)
            call arrays % incrementCentroid(cIdx, rC)
            call arrays % incrementMoments(cIdx, matScore)
            call arrays % hitCell(cIdx)

          call arrays % unsetLock(cIdx)
      
        end if

      ! Handle void cells. Assume flat source.
      ! Does not accumulate geometric info.
      else
        
        ! Accumulate to scalar flux
        if (activeRay) then
      
          len_2 = lenFlt * one_two
          !$omp simd
          do g = 1, nG
            inc(g) = lenFlt * (angular(g) + source(g) * len_2)
          end do

          call arrays % setLock(cIdx)
            call arrays % getFluxPointer(cIdx, scalar)
            !$omp simd
            do g = 1, nG
              scalar(g) = scalar(g) + inc(g)
            end do
            call arrays % incrementVolume(cIdx, length)
            call arrays % hitCell(cIdx)
          call arrays % unsetLock(cIdx)
      
        end if

        !$omp simd
        do g = 1, nG
          angular(g) = angular(g) + source(g) * lenFlt
        end do

      end if

      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, nG
          angular(g) = 0.0_defFlt
        end do
      end if

    end do

  end subroutine transportSweepLinearIso
  
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
