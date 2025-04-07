module rayVolPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,    only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,   only : FNV_1
  use dictionary_class,     only : dictionary
  use rng_class,            only : RNG
  use physicsPackage_inter, only : physicsPackage

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : geometry, distCache
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx, &
                                             gr_kill    => kill
  use geometryFactory_func,           only : new_geometry

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat, mm_matName => matName
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_getMatNames => getMatNames, &
                                             ndReg_kill        => kill

  implicit none
  private

  ! Parameters
  integer(shortInt), parameter :: SCORE = 1, CSUM = 2, CSUM2 = 3

  !!
  !! Physics package to perform ray-tracking based volume calculation
  !!
  !! Calculates relative volume of diffrent materials in the problem by performing
  !! random ray tracing in the geometry. The volume is normalised that the total domain
  !! volume is 1.0.
  !!
  !! Rays travel by a random, exponentially distributed distance with a user-defined mean free
  !! path. After each segment, the ray scatters isotropically. At each scattering event
  !! ray may be terminated with provided probability `abs_prob`. Ray will also be killed
  !! if it has reached the OUTSIDE material.
  !!
  !! Simulation can be performed in the ROBUST mode, which is intended to be used for debugging.
  !! If it is enabled, for each ray segment the material idx at the middle of the segment is
  !! obtained and checked against what is stored in coords. If they do not match it means that
  !! for some reason, the geometry stopped to track material composition correctly for the ray.
  !!
  !! This Physics Package exists to serve as a geometry debugging and benchmarking tool.
  !! The calculation it preforms is unlikley to serve any practical purpose.
  !!
  !! Sample Input Dictionary:
  !!   PP {
  !!     type rayVolPhysicsPackage;
  !!     mfp 0.3;       // Mean length of ray segments
  !!     abs_prob 0.1;  // Ray absorbtion probability after each segment
  !!     pop 2000;      // Number of rays per cycle
  !!     cycles 100;    // Number of cycles
  !!     robust 1;      // 1 for true; 0 for false; Enable robust mode
  !!     cache  1;      // 1 for true; 0 for false; Enable distance caching
  !!     #seed 86868;#  // Optional RNG seed
  !!     geometry {<Geometry Definition>}
  !!     nuclearData {<Nuclear data definition. Requires material names only>}
  !!   }
  !!
  !! Private Members
  !!   geom      -> Pointer to the geometry
  !!   geomIdx   -> Index of the geometry in geometry Registry
  !!   rand      -> Random number generator
  !!   timerMain -> Index of the timer defined to measure calculation time
  !!   mfp       -> Mean length of the ray segment
  !!   abs_prob  -> Ray absorption probability after every segment
  !!   N_cycles  -> Number of cycles
  !!   robust    -> Flag to enable/disable robust mode
  !!   cache     -> Flag to enable/disable distance caching
  !!   res       -> Array to accumulate total track length in each material. Contains
  !!     score (SCORE), cumulative sum over cycles (CSUM) and cumulative sume of squares over
  !!     cycles (CSUM2)
  !!   totDist   -> Accumulator for total track distance in a single cycle
  !!   ray_speed -> Cumulative sume and cumulative sum of squares over cycles for the mean speed
  !!     of ray. Result in meters over CPU second.
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: rayVolPhysicsPackage
    private
    ! Components
    class(geometry), pointer :: geom
    integer(shortInt)        :: geomIdx   = 0
    type(RNG)                :: rand
    integer(shortInt)        :: timerMain = 0

    ! Settings
    real(defReal)      :: mfp       = ZERO
    real(defReal)      :: abs_prob  = ZERO
    integer(shortInt)  :: pop       = 0
    integer(shortInt)  :: N_cycles  = 0
    logical(defBool)   :: robust    = .false.
    logical(defBool)   :: cache     = .false.

    ! Results space
    real(defReal), dimension(:,:), allocatable :: res
    real(defReal)                              :: totDist   = ZERO
    real(defReal), dimension(2)                :: ray_speed = ZERO

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: run
    procedure :: kill

    ! Private procedures
    procedure, private :: cycles
    procedure, private :: trackRay
    procedure, private :: printResults
    procedure, private :: printSettings
  end type rayVolPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(rayVolPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)           :: dict
    integer(shortInt)                          :: seed_temp
    integer(longInt)                           :: seed
    character(10)                              :: time
    character(8)                               :: date
    character(:),allocatable                   :: string
    class(dictionary),pointer                  :: tempDict
    character(nameLen)                         :: geomName
    character(100), parameter :: Here = 'init (rayVolPhysicsPackage_class.f90)'

    ! Load settings
    call dict % get(self % mfp, 'mfp')
    call dict % get(self % abs_prob, 'abs_prob')
    call dict % get(self % pop, 'pop')
    call dict % get(self % N_cycles, 'cycles')
    call dict % get(self % robust, 'robust')
    call dict % get(self % cache, 'cache')

    ! Check settings
    if (self % mfp < ZERO) then
      call fatalError(Here, 'Was given -ve mean free path (mfp): '//numToChar(self % mfp))

    else if (self % abs_prob <= ZERO .or. self % abs_prob > ONE) then
      call fatalError(Here, 'Absorbtion probability is outside valid range &
                            &of 0.0-1.0: '//numToChar(self % abs_prob))
    end if

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Initialise RNG
    if( dict % isPresent('seed')) then
      call dict % get(seed_temp,'seed')

    else
      ! Obtain time string and hash it to obtain random seed
      call date_and_time(date, time)
      string = date // time
      call FNV_1(string,seed_temp)

    end if
    seed = seed_temp
    call self % rand % init(seed)

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'rayCalcGeom'
    call new_geometry(tempDict, geomName)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Allocate results space
    allocate(self % res(mm_nMat(), 3))
    self % res = ZERO
    self % totDist = ZERO
    self % ray_speed = ZERO

  end subroutine init

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(rayVolPhysicsPackage), intent(inout) :: self

    call self % printSettings()
    call self % cycles(self % rand)
    call self % printResults()

  end subroutine run

  !!
  !! Perform cycles of the stochatic volume calculation with ray tracing
  !!
  !! Randomly places the starting point based on uniform distribution.
  !!
  !! Args:
  !!   rand [inout] -> Initialised random number generator
  !!
  !! NOTE:
  !!   RNG needs to be given as an argument `class(RNG)` to prevent inlining. Compiler (gcc 8.3)
  !!   produced erroneous code withou it. Same random number would be produced for diffrent calls
  !!   of `get` function.
  !!
  subroutine cycles(self, rand)
    class(rayVolPhysicsPackage), intent(inout) :: self
    class(RNG), intent(inout)                  :: rand
    type(coordList)                     :: coords
    real(defReal), dimension(3)         :: rand3, bottom, top
    real(defReal), dimension(3)         :: r, u
    real(defReal)                       :: mu, phi
    integer(shortInt)                   :: gen, ray, matIdx, uniqueId, i
    type(RNG), save                     :: pRNG
    real(defReal)                       :: elapsed_T, end_T, T_toEnd, av_speed, cycle_T
    character(100), parameter :: Here = 'cycles (rayVolPhysicsPackage_class.f90)'
    !$omp threadprivate(pRNG)

    !$omp parallel
    pRNG = rand
    !$omp end parallel

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Get lower an upper corner of bounding box
    associate (aabb => self % geom % bounds())
      bottom = aabb(1:3)
      top    = aabb(4:6)
    end associate

    ! Perform clculation
    do gen = 1, self % N_cycles
      !$omp parallel do private(r, u, mu, phi, i, rand3, matIdx, uniqueID, coords)
      do ray = 1, self % pop

        ! Set seed
        call pRNG % stride( (gen-1) * self % pop + ray )

        ! Find starting point that is inside the geometry
        i = 0
        mu = TWO * rand % get() - ONE
        phi = TWO_PI * rand % get()
        u = rotateVector([ONE, ZERO, ZERO], mu, phi)

        rejection : do
          rand3(1) = pRNG % get()
          rand3(2) = pRNG % get()
          rand3(3) = pRNG % get()
          r = bottom + (top - bottom) * rand3

          ! Exit if point is inside the geometry
          call self % geom % whatIsAt(matIdx, uniqueId, r, u)
          if (matIdx /= OUTSIDE_MAT) exit rejection

          i = i + 1
          if (i > 1000) then
            call fatalError(Here, 'Infinite loop when searching ray start in the geometry.')
          end if
        end do rejection

        ! Place in the geometry & process the ray
        call coords % init(r, u)
        call self % geom % placeCoord(coords)
        call self % trackRay(coords, pRNG)

      end do
      !$omp end parallel do

      ! Calculate times
      call timerStop(self % timerMain)
      cycle_T = timerTime(self % timerMain) - elapsed_T
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(self % N_cycles, defReal) * elapsed_T / gen
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Calculate average tracking speed
      av_speed = self % totDist / cycle_T * 1.0E-3_defReal

      ! Display progress
      call printFishLineR(gen)
      print *
      print *, 'Cycle: ', numToChar(gen), ' of ', numToChar(self % N_cycles)
      print *, 'Pop: ', numToChar(self % pop)
      print '(A, ES12.5)', ' Av. Ray speed: [m/s]: ', av_speed
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))

      ! Process scores
      self % res(:, SCORE) = self % res(:, SCORE) / self % totDist
      self % res(:, CSUM)  = self % res(:, CSUM)  + self % res(:, SCORE)
      self % res(:, CSUM2) = self % res(:, CSUM2) + self % res(:, SCORE)**2
      self % res(:, SCORE) = ZERO
      self % totDist = ZERO

      ! Average ray speed
      self % ray_speed(1) = self % ray_speed(1) + av_speed
      self % ray_speed(2) = self % ray_speed(2) + av_speed**2

    end do

  end subroutine cycles

  !!
  !! Track single ray through the geometry
  !!
  !! Seperate this functionality in separate function to improve
  !! clarity of `cycles` procedure.
  !!
  !! Moves ray through the geometry untill it leaks or is absorbed.
  !! Accumulates track-length information along the way.
  !!
  !! Args:
  !!   coords [inout] -> Coordinates of ray placed in the geometry.
  !!   rand [inout] -> Initialised random number generator
  !!
  subroutine trackRay(self, coords, rand)
    class(rayVolPhysicsPackage), intent(inout) :: self
    class(RNG), intent(inout)                  :: rand
    type(coordList), intent(inout)             :: coords
    real(defReal)                              :: dist, mu, phi, maxDist, rn
    real(defReal), dimension(3)                :: r, r_pre, u_pre
    integer(shortInt)                          :: event, matIdx, uniqueId, mat_mid, unique_mid
    type(distCache)                            :: cache_space
    character(100), parameter :: Here = 'trackRay (rayVolPhysicsPackage_class.f90)'

    ! Keep compiler happy
    r_pre = ZERO
    r = ZERO
    u_pre = ZERO

    hist : do
      ! Sample distance
      dist = -log(rand % get()) * self % mfp

      event = LOST_EV
      do while (event /= COLL_EV)
        ! Save pre-movement state
        matIdx = coords % matIdx
        uniqueID = coords % uniqueID
        maxDist = dist
        if (self % robust) then
          r_pre = coords % lvl(1) % r
          u_pre = coords % lvl(1) % dir
        end if

        ! Move in geometry
        if (self % cache) then
          call self % geom % move_withCache(coords, dist, event, cache_space)

        else
          call self % geom % move(coords, dist, event)

        end if

        ! If robust verify matIdx in the mid point
        if (self % robust) then
          r = r_pre + u_pre * HALF * dist
          call self % geom % whatIsAt(mat_mid, unique_mid, r, u_pre)

          if (matIdx /= mat_mid ) then
            print *, "EVENT: ", event
            print *, "PRE MOVE: ", r_pre
            print *, "DIR", coords % lvl(1) % dir
            print *, "POST MOVE", coords % lvl(1) % r
            print *, "MAT CHECKED AT:", r
            print *, "WITH DIRECTION:", u_pre
            print *, "CORRECT MAT:", mat_mid
            print *, "HAS MAT:", matIdx
            call fatalError(Here, 'Ray has lost correct material.')
          end if

        end if

        ! Score result
        !$omp atomic
        self % totDist = self % totDist + dist
        if (matIdx /= VOID_MAT) then
          !$omp atomic
          self % res(matIdx, SCORE) = self % res(matIdx, SCORE) + dist
        end if

        ! Set to remaining distance
        dist = maxDist - dist
      end do

      ! Kill the ray
      rn = rand % get()
      if (self % abs_prob > rn .or. coords % matIdx == OUTSIDE_MAT) exit hist

      ! Scatter the ray
      mu = TWO * rand % get() - ONE
      phi = TWO_PI * rand % get()
      call coords % rotate(mu, phi)

    end do hist

  end subroutine trackRay

  !!
  !! Output calculation results to the console
  !!
  !! Convert cumulative sums to mean and absolute standard deviation and
  !! print them to the console.
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(rayVolPhysicsPackage), intent(in) :: self
    real(defReal)                           :: mean, SD, var
    real(defReal)                           :: V, V_SD
    integer(shortInt)                       :: i

    ! Calculate speed and its SD
    V = self % ray_speed(1) / self % N_cycles
    V_SD = self % ray_speed(2) / self % N_cycles - V**2
    V_SD = ONE/(self % N_cycles - 1) * sqrt(V_SD)

    print *
    print '(A, ES12.5, A, ES12.5)', " Ray speed [m/s]: ", V, " +/- ", V_SD
    print *, "RELATIVE VOLUME FOR MATERIALS: "
    do i = 1, mm_nMat()
      mean = self % res(i, CSUM) / self % N_cycles
      var = self % res(i, CSUM2) / self % N_cycles - mean**2
      SD = ONE/(self % N_cycles - 1) * sqrt(var)
      print '(A, A, A, ES12.5, A, ES12.5)', " Material: ", mm_matName(i), " Vol", mean, " +/-", SD
    end do

  end subroutine printResults

  !!
  !! Print settings of the Ray-tracking volume calculation
  !!
  !! Args:
  !!   None
  !!
  subroutine printSettings(self)
    class(rayVolPhysicsPackage), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ RAY-TRACING RELATIVE VOLUME CALCULATION /\/\"
    print *, "Total Cycles:    ", numToChar(self % N_cycles)
    print *, "Rays per cycle: ", numToChar(self % pop)
    print *, "Initial RNG Seed:   ", numToChar(self % rand % getSeed())
    print *
    print *, repeat("<>", MAX_COL/2)

  end subroutine printSettings

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(rayVolPhysicsPackage), intent(inout) :: self

    ! Clean Nuclear Data & Geometries
    call gr_kill()
    call ndreg_kill()

    ! Clean contents
    self % geom    => null()
    self % geomIdx = 0
    !call self % rand % kill()
    self % timerMain = 0

    self % mfp      = ZERO
    self % abs_prob = ZERO
    self % pop      = 0
    self % N_cycles = 0
    self % robust   = .false.

    if (allocated(self % res)) deallocate(self % res)
    self % totDist   = ZERO
    self % ray_speed = ZERO

  end subroutine kill

end module rayVolPhysicsPackage_class
