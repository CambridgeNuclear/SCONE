module randomRayPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use exponentialRA_func,             only : exponential
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile
  use rng_class,                      only : RNG
  use physicsPackage_inter,           only : physicsPackage

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : geometry, distCache
  use geometryStd_class,              only : geometryStd
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_addGeom => addGeom, &
                                             gr_geomIdx  => geomIdx, gr_kill    => kill

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat            => nMat, mm_matName => matName
  use nuclearDataReg_mod,             only : ndReg_init         => init, &
                                             ndReg_getMatNames  => getMatNames, &
                                             ndReg_activate     => activate, &
                                             ndReg_kill         => kill, &
                                             ndReg_getNeutronMG => getNeutronMG
  use materialHandle_inter,           only : materialHandle
  use mgNeutronDatabase_inter,        only : mgNeutronDatabase
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use baseMgNeutronMaterial_class,    only : baseMgNeutronMaterial, baseMgNeutronMaterial_CptrCast
  
  ! Visualisation
  use visualiser_class,               only : visualiser

  ! Random ray
  use ray_class,                      only : ray

  implicit none
  private

  ! Parameter for when to skip a tiny volume
  real(defReal), parameter :: volume_tolerance = 1.0E-10

  !!
  !! Physics package to perform The Random Ray Method (TRRM) eigenvalue calculations
  !!
  !! Tracks rays across the geometry, attenuating their flux. After some dead length,
  !! rays begin scoring to estimates of the scalar flux and volume. Each ray has a
  !! uniform termination length, after which it is stopped and the next ray is tracked.
  !! Once all rays have been tracked, a cycle concludes and fluxes, sources, and keff
  !! are updated.
  !!
  !! Both inactive and active cycles occur, as in Monte Carlo. These can be terminated
  !! after a specified number of iterations or on reaching some chosen convergence
  !! criterion.
  !!
  !! Calculates relative volume of diffrent materials in the problem by performing
  !! random ray tracing in the geometry. The volume is normalised that the total domain
  !! volume is 1.0.
  !!
  !! IMPORTANT N.B.: Geometry type must be extended! If shrunk, results may be dubious!
  !! This is because spatial discretisation is determined by the number of unique cells in the
  !! geometry.
  !! Also, this is obviously for multi-group calculations only.
  !!
  !! Sample Input Dictionary:
  !!   PP {
  !!     type randomRayPhysicsPackage;
  !!     dead 10;            // Dead length where rays do not score to scalar fluxes
  !!     termination 100;    // Length a ray travels before it is terminated
  !!     rays 1000;          // Number of rays to sample per iteration
  !!     inactive 100;       // Number of convergence cycles (would use accum and std otherwise)
  !!     active 200;         // Number of scoring cycles (would use eps otherwise)
  !!     !NOT YET SUPPORTED!
  !!     accum 15;           // Accumulation cycles for exiting inactive cycles and estimating keff std
  !!     std 30;             // Standard deviation below which inactive cycles are exited
  !!     eps 1e-5;           // RMS error on scalar fluxes below which simulation concludes
  !!     !NOT YET SUPPORTED!
  !!     #seed 86868;#       // Optional RNG seed
  !!     geometry {<Geometry Definition>}
  !!     nuclearData {<Nuclear data definition>}
  !!   }
  !!
  !! Private Members
  !!   geom        -> Pointer to the geometry.
  !!   geomIdx     -> Index of the geometry in geometry Registry.
  !!   top         -> Top co-ordinates of the geometry bounding box.
  !!   bottom      -> Bottom co-ordinates of the geometry bounding box.
  !!   rand        -> Random number generator.
  !!   timerMain   -> Index of the timer defined to measure calculation time.
  !!   mgData      -> MG database. Calculation obviously cannot be run in CE.
  !!   nG          -> Number of energy groups, kept for convenience.
  !!   nCells      -> Number of unique cells in the geometry, kept for convenience.
  !!
  !!   termination -> Distance a ray can travel before it is terminated
  !!   dead        -> Distance a ray must travel before it becomes active
  !!   pop         -> Number of rays to track per cycle
  !!   inactive    -> Number of inactive cycles to perform
  !!   active      -> Number of active cycles to perform
  !!   !NOT YET SUPPORTED!
  !!   criterionA  -> Should active cycles terminate on some convergence criterion?
  !!   criterionI  -> SHould inactive cycles terminate on some convergence criterion?
  !!   accum       -> Number of cycles over which to accumulate keff in a moving window
  !!   std         -> Standard deviation tolerance criterion for terminating inactive
  !!   eps         -> RMS tolerance criterion for terminating active
  !!   !NOT YET SUPPORTED!
  !!
  !!   keff        -> Estimated value of keff
  !!   keffScore   -> Vector holding cumulative keff score and keff^2 score
  !!   scalarFlux  -> Array of scalar flux values of length = nG * nCells
  !!   prevFlux    -> Array of previous scalar flux values of length = nG * nCells
  !!   fluxScore   -> Array of scalar flux values and squared values to be reported 
  !!                  in results, dimension =  [nG * nCells, 2]
  !!   source      -> Array of neutron source values of length = nG * nCells
  !!   volume      -> Array of stochastically estimated cell volumes of length = nCells
  !!   cellHit     -> Array tracking whether given cells have been hit during tracking
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: randomRayPhysicsPackage
    private
    ! Components
    class(geometryStd), pointer           :: geom
    integer(shortInt)                     :: geomIdx     = 0
    real(defReal), dimension(3)           :: top         = ZERO
    real(defReal), dimension(3)           :: bottom      = ZERO
    type(RNG)                             :: rand
    class(baseMgNeutronDatabase), pointer :: mgData      => null()
    integer(shortInt)                     :: nG          = 0
    integer(shortInt)                     :: nCells      = 0
    real(defReal)                         :: lengthPerIt = ZERO

    ! Settings
    real(defReal)      :: termination = ZERO
    real(defReal)      :: dead        = ZERO
    integer(shortInt)  :: pop         = 0
    integer(shortInt)  :: inactive    = 0
    integer(shortInt)  :: active      = 0
    !logical(defBool)   :: criterionA  = .FALSE.
    !logical(defBool)   :: criterionI  = .FALSE.
    !integer(shortInt)  :: accum       = 0
    !real(defReal)      :: std         = ZERO
    !real(defReal)      :: eps         = ZERO
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    logical(defBool)   :: plotResults = .FALSE.
    type(visualiser)   :: viz

    ! Results space
    real(defReal)                                :: keff
    real(defReal), dimension(2)                  :: keffScore
    real(defReal), dimension(:), allocatable     :: scalarFlux
    real(defReal), dimension(:), allocatable     :: prevFlux
    real(defReal), dimension(:,:), allocatable   :: fluxScores
    real(defReal), dimension(:), allocatable     :: source
    real(defReal), dimension(:), allocatable     :: volume
    real(defReal), dimension(:), allocatable     :: volumeTracks

    ! Tracking cell properites
    integer(shortInt), dimension(:), allocatable :: cellHit
    logical(defBool), dimension(:), allocatable  :: cellFound
    real(defReal), dimension(:,:), allocatable   :: cellPos

    ! Timer bins
    integer(shortInt) :: timerMain
    integer(shortInt) :: timerTransport
    real (defReal)    :: time_transport = ZERO
    real (defReal)    :: CPU_time_start
    real (defReal)    :: CPU_time_end

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: run
    procedure :: kill

    ! Private procedures
    procedure, private :: cycles
    procedure, private :: initialiseRay
    procedure, private :: moveRay
    procedure, private :: moveRayCache
    procedure, private :: transportSweep
    procedure, private :: calculateSources
    procedure, private :: calculateKeff
    procedure, private :: normaliseFluxAndVolume
    procedure, private :: resetFluxes
    procedure, private :: accumulateFluxAndKeffScores
    procedure, private :: finaliseFluxAndKeffScores
    procedure, private :: printResults
    procedure, private :: printSettings
  end type randomRayPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(randomRayPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)              :: dict
    integer(shortInt)                             :: seed_temp
    integer(longInt)                              :: seed
    character(10)                                 :: time
    character(8)                                  :: date
    character(:),allocatable                      :: string
    class(dictionary),pointer                     :: tempDict, graphDict
    class(mgNeutronDatabase),pointer              :: db
    character(nameLen)                            :: geomName, graphType, nucData
    class(geometry), pointer                      :: geom
    type(outputFile)                              :: test_out
    character(100), parameter :: Here = 'init (randomRayPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)
    
    ! Load settings
    call dict % get( nucData, 'XSdata')
    call dict % get(self % termination, 'termination')
    call dict % get(self % dead, 'dead')
    call dict % get(self % pop, 'pop')
    !call dict % get(self % accum, 'accum')
    !if (dict % isPresent('active')) then
    call dict % get(self % active, 'active')
    !  self % eps = ZERO
    !  self % criterionA = .FALSE.
    !else
    !  call dict % get(self % eps, 'eps')
    !end if
    !if (dict % isPresent('inactive')) then
    call dict % get(self % inactive, 'inactive')
    !  self % std = ZERO
    !  self % criterionI = .FALSE.
    !else
    !  call dict % get(self % std, 'std')
    !end if
    
    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be cought early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Check settings
    if (self % termination <= ZERO) call fatalError(Here, 'Ray termination distance (termination) is less than or equal to zero.')
    if (self % pop < 1) call fatalError(Here, 'Must have 1 or more rays (pop).')
    !if (self % accum < 1) call fatalError('Must have 1 or more accumulation cycles (accum).')
    !if (self % std < ZERO) call fatalError('Convergence standard deviation (std) cannot be less than zero.')
    !if (self % eps < ZERO) call fatalError('RMS tolerance (eps) cannot be less than zero.')

    ! Dead length can be less than zero but will be reset to zero if so
    if (self % dead < ZERO) then
      self % dead = ZERO
      print *,'Warning: Dead length of rays (dead) was negative. This has been set to 0 instead.'
    end if

    ! Ensure termination length is longer than dead length
    if (self % termination <= self % dead) call fatalError(Here,&
            'Ray termination length must be greater than ray dead length')

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
    geomName = 'randomRayGeom'
    call gr_addGeom(geomName, tempDict)
    self % geomIdx = gr_geomIdx(geomName)
    geom    => gr_geomPtr(self % geomIdx)

    ! Ensure geometry is geometryStd
    select type(geom)
      type is (geometryStd)
        self % geom => geom
      class default
        call fatalError(Here,'Unrecognised geometry type')
    end select

    ! Ensure that geometry graph is extended
    graphDict => tempDict % getDictPtr('graph')
    call graphDict % get(graphType,'type')
    if (graphType /= 'extended') call fatalError(Here,&
            'Geometry graph type must be "extended" for random ray calculations.')

    ! Activatee nuclear data
    call ndReg_activate(P_NEUTRON_MG, nucData, self % geom % activeMats())

    ! Ensure that nuclear data is multi-group
    db => ndReg_getNeutronMG()
    if (.NOT. associated(db)) call fatalError(Here,&
            'No MG nuclear database was constructed')

    ! Ensure nuclear data is baseMgNeutronDatabase
    select type(db)
      type is (baseMgNeutronDatabase)
        self % mgData => db
      class default
        call fatalError(Here,'Unrecognised MG database type')
    end select

    ! Store number of energy groups for convenience
    self % nG = self % mgData % nGroups()

    ! Get lower and upper corner of bounding box
    associate (aabb => self % geom % bounds())
      self % bottom = aabb(1:3)
      self % top    = aabb(4:6)
    end associate
    
    ! Call visualisation
    if (dict % isPresent('viz')) then
      print *, "Initialising visualiser"
      tempDict => dict % getDictPtr('viz')
      call self % viz % init(geom, tempDict)
      print *, "Constructing visualisation"
      call self % viz % makeViz()
      call self % viz % kill()
    endif
    
    ! Check for results plotting and initialise VTK
    call dict % getOrDefault(self % plotResults,'plot',.FALSE.)
    if (self % plotResults) then
      ! Initialise a visualiser to be used when results are available
      print *, "Initialising results visualiser"
      tempDict => dict % getDictPtr('viz')
      call self % viz % init(geom, tempDict)
      print *, "Constructing geometry visualisation"
      call self % viz % initVTK()
    end if

    ! Store number of cells in geometry for convenience
    self % nCells = self % geom % numberOfCells()

    ! Allocate results space
    allocate(self % scalarFlux(self % nCells * self % nG))
    allocate(self % prevFlux(self % nCells * self % nG))
    allocate(self % fluxScores(self % nCells * self % nG, 2))
    allocate(self % source(self % nCells * self % nG))
    allocate(self % volume(self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(self % nCells, 3))
    
    ! Set active length traveled per iteration
    self % lengthPerIt = (self % termination - self % dead) * self % pop

  end subroutine init

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(randomRayPhysicsPackage), intent(inout) :: self

    call self % printSettings()
    call self % cycles()
    call self % printResults()

  end subroutine run

  !!
  !! Perform cycles of The Random Ray Method.
  !!
  !! Randomly places the ray starting point and direction uniformly.
  !! Rays are tracked until they reach some specified termination length.
  !! During tracking, fluxes are attenuated (and adjusted according to BCs),
  !! scoring to fluxes and volume estimates when the ray has surpassed its 
  !! specified dead length.
  !!
  !! Inactive and active iterations occur, terminating subject either to 
  !! given criteria or when a fixed number of iterations has been passed.
  !!
  !! Args:
  !!   rand [inout] -> Initialised random number generator
  !!
  !! NOTE:
  !!   RNG needs to be given as an argument `class(RNG)` to prevent inlining. Compiler (gcc 8.3)
  !!   produced erroneous code withou it. Same random number would be produced for diffrent calls
  !!   of `get` function.
  !!
  subroutine cycles(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    type(ray), save                               :: r
    type(RNG), target, save                       :: pRNG
    real(defReal)                                 :: hitRate
    real(defReal)                                 :: elapsed_T, end_T, T_toEnd, transport_T
    logical(defBool)                              :: stoppingCriterion, isActive
    integer(shortInt)                             :: i, itInac, itAct, it
    integer(longInt), save                        :: ints
    integer(longInt)                              :: integrations
    !$omp threadprivate(pRNG, r, ints)

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Initialise fluxes 
    self % keff       = ONE
    self % scalarFlux = ZERO
    self % prevFlux   = ONE
    self % fluxScores = ZERO
    self % keffScore  = ZERO
    self % source     = ZERO

    ! Initialise other results
    self % cellHit      = 0
    self % volume       = ZERO
    self % volumeTracks = ZERO
    
    ! Initialise cell information
    self % cellFound = .FALSE.
    self % cellPos = -INFINITY

    ! Stopping criterion is initially on flux convergence or number of convergence iterations.
    ! Will be replaced by RMS error in flux or number of scoring iterations afterwards.
    itInac = 0
    itAct  = 0
    isActive = .FALSE.
    stoppingCriterion = .TRUE.

    ! Power iteration
    do while( stoppingCriterion )
      
      if (isActive) then
        itAct = itAct + 1
      else
        itInac = itInac + 1
      end if
      it = itInac + itAct

      call self % calculateSources()
    
      ! Reset and start transport timer
      call timerReset(self % timerTransport)
      call timerStart(self % timerTransport)
      integrations = 0
      
      !$omp parallel do schedule(dynamic) reduction(+: integrations)
      do i = 1, self % pop

        ! Set seed
        pRNG = self % rand
        call pRNG % stride(i)
        r % pRNG => pRNG 

        ! Set ray attributes
        call self % initialiseRay(r)

        ! Transport ray until termination criterion met
        call self % transportSweep(r,ints)
        integrations = integrations + ints

      end do
      !$omp end parallel do
      
      call timerStop(self % timerTransport)

      ! Update RNG on master thread
      call self % rand % stride(self % pop + 1)

      ! Normalise flux estimate and combines with source
      call self % normaliseFluxAndVolume(it)

      ! Calculate new k
      call self % calculateKeff()

      ! Accumulate flux scores
      if (isActive) call self % accumulateFluxAndKeffScores()

      ! Calculate proportion of cells that were hit
      hitRate = real(sum(self % cellHit),defReal) / self % nCells
      self % cellHit = 0

      ! Evaluate stopping criterion for active or inactive iterations
      if (isActive) then
        !if (self % criterionA) then
        !  rms = self % calculateRMS()
        !  stopping_criterion = rms > self % eps
        !else
        stoppingCriterion = (itAct < self % active)
        !end if
      else
        !if (self % criterionI) then
        !  isActive = ((it >= self % accum) .AND. (kStd < self % std))
        !else
        isActive = (itInac >= self % inactive)
        !end if
      end if

      ! Set previous iteration flux to scalar flux
      ! and zero scalar flux
      call self % resetFluxes()

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)
      transport_T = timerTime(self % timerTransport)
      self % time_transport = self % time_transport + transport_T

      ! Predict time to end
      end_T = real(self % active + self % inactive, defReal) * elapsed_T / it
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      call printFishLineR(it)
      print *
      print *, 'Iteration: ', numToChar(it), ' of ', numToChar(self % active + self % inactive)
      if(isActive) then
        print *,'Active iterations'
      else
        print *,'Inactive iterations'
      end if
      print *, 'Cell hit rate: ', trim(numToChar(hitRate))
      print *, 'keff: ', trim(numToChar(self % keff))
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      print *, 'Time per integration (ns): ', &
              trim(numToChar(transport_T*10**9/(self % nG * integrations)))

    end do

    ! Finalise flux scores
    call self % finaliseFluxAndKeffScores(itAct)

  end subroutine cycles

  !!
  !! Initialises rays: samples initial position and direction,
  !! performs the build operation, and sets the initial flux
  !!
  subroutine initialiseRay(self, r)
    class(randomRayPhysicsPackage), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    real(defReal)                                 :: mu, phi
    real(defReal), dimension(3)                   :: u, rand3, x
    integer(shortInt)                             :: i, matIdx, g, idx, cIdx
    character(100), parameter :: Here = 'initialiseRay (randomRayPhysicsPackage_class.f90)'

    i = 0
    mu = TWO * r % pRNG % get() - ONE
    phi = TWO_PI * r % pRNG % get()
    u = rotateVector([ONE, ZERO, ZERO], mu, phi)

    rejection : do
      rand3(1) = r % pRNG % get()
      rand3(2) = r % pRNG % get()
      rand3(3) = r % pRNG % get()
      x = self % bottom + (self % top - self % bottom) * rand3

      ! Exit if point is inside the geometry
      call self % geom % whatIsAt(matIdx, cIdx, x, u)
      if (matIdx /= OUTSIDE_MAT) exit rejection

      i = i + 1
      if (i > 5000) then
        call fatalError(Here, 'Infinite loop when searching ray start in the geometry.')
      end if
    end do rejection

    ! Place in the geometry & process the ray
    call r % build(x, u, self % nG)
    call self % geom % placeCoord(r % coords)

    ! Set angular flux to angle average of cell source
    do g = 1, self % nG
      idx = (cIdx - 1) * self % nG + g
      r % flux(g) = self % source(idx)
    end do

    if (.NOT. self % cellFound(cIdx)) then
      !$omp critical 
      self % cellFound(cIdx) = .TRUE.
      self % cellPos(cIdx,:) = x
      !$omp end critical
    end if

  end subroutine initialiseRay

  !!
  !! Geometry handling of the ray: moves the ray across a FSR,
  !! returning the distance travelled and applying boundary conditions
  !! where necessary.
  !!
  !! Accounts for biases that may be introduced by rays moving
  !! significantly past their dead/termination length during an FSR crossing.
  !! The ray is only moved as far as the point where the dead/termination 
  !! length is exceeded and this is the distance which is travelled during 
  !! the operation.
  !!
  !! Reports the vacuum hit to be dealt with after attenuation
  !! has taken place.
  !!
  subroutine moveRay(self, r, length, hitVacuum)
    class(randomRayPhysicsPackage), intent(in) :: self
    type(ray), intent(inout)                   :: r
    real(defReal), intent(out)                 :: length
    logical(defBool), intent(out)              :: hitVacuum
    integer(shortInt)                          :: event

    ! Set maximum flight distance and ensure ray is active
    if (r % length >= self % dead) then
      length = self % termination - r % length 
      r % isActive = .TRUE.
    else
      length = self % dead - r % length
    end if

    ! Perform the movement
    call self % geom % moveRay_noCache(r % coords, length, event, hitVacuum)

    r % length = r % length + length

  end subroutine moveRay

  !!
  !! Geometry handling of the ray: moves the ray across a FSR,
  !! returning the distance travelled and applying boundary conditions
  !! where necessary.
  !!
  !! Accounts for biases that may be introduced by rays moving
  !! significantly past their dead/termination length during an FSR crossing.
  !! The ray is only moved as far as the point where the dead/termination 
  !! length is exceeded and this is the distance which is travelled during 
  !! the operation.
  !!
  !! Reports the vacuum hit to be dealt with after attenuation
  !! has taken place.
  !!
  subroutine moveRayCache(self, r, length, hitVacuum, cache)
    class(randomRayPhysicsPackage), intent(in) :: self
    type(ray), intent(inout)                   :: r
    real(defReal), intent(out)                 :: length
    logical(defBool), intent(out)              :: hitVacuum
    type(distCache), intent(inout)             :: cache
    integer(shortInt)                          :: event

    ! Set maximum flight distance and ensure ray is active
    if (r % length >= self % dead) then
      length = self % termination - r % length 
      r % isActive = .TRUE.
    else
      length = self % dead - r % length
    end if

    ! Perform the movement
    call self % geom % moveRay_withCache(r % coords, length, event, cache, hitVacuum)

    r % length = r % length + length

  end subroutine moveRayCache

  !!
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume.
  !! Records the number of integrations/ray movements.
  !!
  subroutine transportSweep(self, r, ints)
    class(randomRayPhysicsPackage), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    integer(longInt), intent(out)                 :: ints
    integer(shortInt)                             :: matIdx, g, cIdx, idx
    real(defReal)                                 :: attenuate, length, delta, total
    logical(defBool)                              :: hitVacuum
    type(distCache)                               :: cache
    class(baseMgNeutronMaterial), pointer         :: mat
    class(materialHandle), pointer                :: matPtr
    
    ints = 0
    do while (r % length < self % termination)

      ! Get material and cell the ray is moving through
      matIdx  = r % coords % matIdx
      cIdx    = r % coords % uniqueID
      matPtr  => self % mgData % getMaterial(matIdx)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)

      ! Remember new cell positions
      if (.NOT. self % cellFound(cIdx)) then
        !$omp critical 
        self % cellFound(cIdx) = .TRUE.
        self % cellPos(cIdx,:) = r % rGlobal()
        !$omp end critical
      end if
          
      ! Move ray
      call self % moveRayCache(r, length, hitVacuum, cache)
      !call self % moveRay(r, length, hitVacuum) ! For debugging

      ints = ints + 1
 
      do g = 1, self % nG
        
        total = mat % getTotalXS(g, self % rand)
            
        ! Calculate delta
        idx = (cIdx - 1) * self % nG + g
        attenuate = exponential(total * length)
        !attenuate = ONE - exp(-SigmaT * length) ! For debugging
        delta = (r % flux(g) - self % source(idx)) * attenuate

        ! Accumulate scalar flux
        if (r % isActive) then
          !$omp atomic
          self % scalarFlux(idx) = self % scalarFlux(idx) + delta 
        end if

        ! Update flux
        r % flux(g) = r % flux(g) - delta

      end do

      ! Accumulate cell volume estimates
      if (r % isActive) then
        self % cellHit(cIdx) = 1

        !$omp atomic
        self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length
      end if

      ! Check for a vacuum hit
      if (hitVacuum) r % flux = ZERO

    end do

  end subroutine transportSweep

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(randomRayPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defReal)                                 :: norm, normVol
    real(defReal), save                           :: total
    integer(shortInt), save                       :: g, matIdx, idx
    integer(shortInt)                             :: cIdx
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(mat, matPtr, total, idx, g, matIdx)

    norm = ONE / self % lengthPerIt
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      matPtr => self % mgData % getMaterial(matIdx)
      mat    => baseMgNeutronMaterial_CptrCast(matPtr)
      
      ! Update volume due to additional rays
      self % volume(cIdx) = self % volumeTracks(cIdx) * normVol

      do g = 1, self % nG

        total =  mat % getTotalXS(g, self % rand)
        idx   = self % nG * (cIdx - 1) + g

        if (self % volume(cIdx) > volume_tolerance) then
          self % scalarFlux(idx) = self % scalarFlux(idx) * norm / ( total * self % volume(cIdx))
        end if
        self % scalarFlux(idx) = self % scalarFlux(idx) + self % source(idx)

      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume
  
  !!
  !! Calculate sources in all cells and energy groups
  !!
  subroutine calculateSources(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    real(defReal)                                 :: ONE_KEFF
    real(defReal), save                           :: scatter, fission
    real(defReal), save                           :: nuFission, total, chi, scatterXS
    integer(shortInt), save                       :: matIdx, g, gIn, idx, idx0
    integer(shortInt)                             :: cIdx
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(mat, matPtr, scatter, fission, scatterXS, nuFission, total, chi, matIdx, g, gIn, idx, idx0)

    ONE_KEFF = ONE / self % keff

    !$omp parallel do schedule(static) 
    do cIdx = 1, self % nCells

      ! Identify material
      matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
      matPtr  => self % mgData % getMaterial(matIdx)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)

      do g = 1, self % nG

        ! Calculate fission and scattering source
        scatter = ZERO
        fission = ZERO

        ! Sum contributions from all energies
        do gIn = 1, self % nG
      
          ! Input index
          idx0 = self % nG * (cIdx - 1) + gIn

          scatterXS = mat % getScatterXS(gIn, g, self % rand)
          nuFission = mat % getNuFissionXS(gIn, self % rand)

          fission = fission + self % prevFlux(idx0) * nuFission
          scatter = scatter + self % prevFlux(idx0) * scatterXS
          
        end do

        ! Output index
        idx = self % nG * (cIdx - 1) + g

        total = mat % getTotalXS(g, self % rand)
        chi   = mat % getChi(g, self % rand)
    
        self % source(idx) = chi * fission * ONE_KEFF + scatter
        self % source(idx) = self % source(idx) / total

      end do

    end do
    !$omp end parallel do

  end subroutine calculateSources

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    real(defReal)                                 :: fissionRate, prevFissionRate
    real(defReal), save                           :: fissLocal, prevFissLocal, nuFission, vol
    integer(shortInt), save                       :: matIdx, g, idx
    integer(shortInt)                             :: cIdx
    class(baseMgNeutronMaterial), pointer, save   :: mat
    class(materialHandle), pointer, save          :: matPtr
    !$omp threadprivate(mat, matPtr, fissLocal, prevFissLocal, nuFission, matIdx, g, idx, vol)

    fissionRate     = ZERO
    prevFissionRate = ZERO
    !$omp parallel do schedule(static) reduction(+: fissionRate, prevFissionRate)
    do cIdx = 1, self % nCells

      ! Identify material
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
      matPtr => self % mgData % getMaterial(matIdx)
      mat    => baseMgNeutronMaterial_CptrCast(matPtr)

      vol = self % volume(cIdx)

      if (vol <= volume_tolerance) cycle

      fissLocal = ZERO
      prevFissLocal = ZERO
      do g = 1, self % nG
      
        nuFission = mat % getNuFissionXS(g, self % rand)
        
        ! Source index
        idx = self % nG * (cIdx - 1) + g
        fissLocal     = fissLocal     + self % scalarFlux(idx) * nuFission
        prevFissLocal = prevFissLocal + self % prevFlux(idx) * nuFission

      end do

      fissionRate     = fissionRate     + fissLocal * vol
      prevFissionRate = prevFissionRate + prevFissLocal * vol

    end do
    !$omp end parallel do

    ! Update k
    self % keff = self % keff * fissionRate / prevFissionRate

  end subroutine calculateKeff

  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !!
  subroutine resetFluxes(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    integer(shortInt)                             :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = ZERO
    end do
    !$omp end parallel do

  end subroutine resetFluxes

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxAndKeffScores(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    real(defReal), save                           :: flux
    integer(shortInt)                             :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = self % scalarFlux(idx)
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) + flux
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) + flux*flux
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) + self % keff
    self % keffScore(2) = self % keffScore(2) + self % keff * self % keff

  end subroutine accumulateFluxAndKeffScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxAndKeffScores(self,it)
    class(randomRayPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    integer(shortInt)                             :: idx
    real(defReal)                                 :: N1, Nm1

    if (it /= 1) then
      Nm1 = ONE/(it - 1)
    else
      Nm1 = ONE
    end if
    N1 = ONE/it

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) * N1
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) * N1
      self % fluxScores(idx,2) = Nm1 *(self % fluxScores(idx,2) - &
            self % fluxScores(idx,1) * self % fluxScores(idx,1)) 
      if (self % fluxScores(idx,2) <= ZERO) then
        self % fluxScores(idx,2) = ZERO
      else
        self % fluxScores(idx,2) = sqrt(self % fluxScores(idx,2))
      end if
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) * N1
    self % keffScore(2) = self % keffScore(2) * N1
    self % keffScore(2) = sqrt(Nm1*(self % keffScore(2) - &
            self % keffScore(1) * self % keffScore(1))) 

  end subroutine finaliseFluxAndKeffScores
  
  !!
  !! Output calculation results to a file
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    type(outputFile)                              :: out
    character(nameLen)                            :: name
    integer(shortInt)                             :: g, cIdx
    integer(shortInt), save                       :: idx
    integer(shortInt),dimension(:),allocatable    :: resArrayShape
    real(defReal), dimension(:), allocatable      :: groupFlux
    !$omp threadprivate(idx)

    call out % init(self % outputFormat)
    
    name = 'seed'
    call out % printValue(self % rand % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Inactive_Cycles'
    call out % printValue(self % inactive,name)

    name = 'Active_Cycles'
    call out % printValue(self % active,name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Total_Transport_Time'
    call out % printValue(self % time_transport,name)
    
    ! Print keff
    name = 'keff'
    call out % startBlock(name)
    call out % printResult(self % keffScore(1), self % keffScore(2), name)
    call out % endBlock()

    ! Print cell volumes
    name = 'volume'
    call out % startBlock(name)
    resArrayShape = [size(self % volume)]
    call out % startArray(name, resArrayShape)
    do cIdx = 1, self % nCells
      call out % addResult(self % volume(cIdx), ZERO)
    end do
    call out % endArray()
    call out % endBlock()

    ! Print cell positions
    name = 'position'
    call out % startBlock(name)
    resArrayShape = [size(self % cellPos)]
    call out % startArray(name, resArrayShape)
    do cIdx = 1, self % nCells
      call out % addResult(self % cellPos(cIdx,1), ZERO)
      call out % addResult(self % cellPos(cIdx,2), ZERO)
      call out % addResult(self % cellPos(cIdx,3), ZERO)
    end do
    call out % endArray()
    call out % endBlock()

    ! Print fluxes
    resArrayShape = [size(self % volume)]
    do g = 1, self % nG
      name = 'flux_g'//numToChar(g)
      call out % startBlock(name)
      call out % startArray(name, resArrayShape)
      do cIdx = 1, self % nCells
        idx = (cIdx - 1)* self % nG + g
        call out % addResult(self % fluxScores(idx,1), self % fluxScores(idx,2))
      end do
      call out % endArray()
      call out % endBlock()
    end do

    call out % writeToFile(self % outputFile)

    ! Send all fluxes and stds to VTK
    if (self % plotResults) then
      allocate(groupFlux(self % nCells))
      do g = 1, self % nG
        name = 'flux_g'//numToChar(g)
        !$omp parallel do
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g
          groupFlux(cIdx) = self % fluxScores(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g = 1, self % nG
        name = 'std_g'//numToChar(g)
        !$omp parallel do
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g
          groupFlux(cIdx) = self % fluxScores(idx,2)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      call self % viz % finaliseVTK
    end if

  end subroutine printResults

  !!
  !! Print settings of the random ray calculation
  !!
  !! Args:
  !!   None
  !!
  subroutine printSettings(self)
    class(randomRayPhysicsPackage), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ RANDOM RAY EIGENVALUE CALCULATION /\/\"
    !if (self % criterionI) then
    !  print *, "Using convergence criterion for inactive"
    !  print *, "Accumulates k over ",numToChar(self % accum),&
    !          " cycles to reach a std. dev. of ",numToChar(self % std)
    !else
    print *, "Using "//numToChar(self % inactive)// " iterations for "&
              //"the inactive cycles"
    !end if
    !if (self % criterionA) then
    !  print *, "Using convergence criterion for active"
    !  print *, "Compares RMS within a tolerance of ",numToChar(self % eps)
    !else
    print *, "Using "//numToChar(self % active)// " iterations for "&
              //"the active cycles"
    !end if
    print *, 
    print *, "Rays per cycle: "// numToChar(self % pop)
    print *, "Ray dead length: "//numToChar(self % dead)
    print *, "Ray termination length: "//numToChar(self % termination)
    print *, "Initial RNG Seed:   "// numToChar(self % rand % getSeed())
    print *,
    print *, "Number of cells in the geometry: "// numToChar(self % nCells)
    print *, "Number of energy groups: "// numToChar(self % nG)
    print *, repeat("<>", MAX_COL/2)

  end subroutine printSettings

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(randomRayPhysicsPackage), intent(inout) :: self

    ! Clean Nuclear Data, Geometry and visualisation
    call gr_kill()
    call ndreg_kill()
    call self % viz % kill()

    ! Clean contents
    self % geom    => null()
    self % geomIdx = 0
    self % timerMain = 0

    self % top       = ZERO
    self % bottom    = ZERO
    self % mgData    => null()
    self % nG        = 0
    self % nCells    = 0

    self % termination = ZERO
    self % dead        = ZERO
    self % pop         = 0
    self % inactive    = 0
    self % active      = 0
    !self % criterionA  = .FALSE.
    !self % criterionI  = .FALSE.
    !self % accum       = 0
    !self % std         = ZERO
    !self % eps         = ZERO

    self % keff        = ZERO
    self % keffScore   = ZERO
    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % source)) deallocate(self % source)
    if(allocated(self % volume)) deallocate(self % volume)
    if(allocated(self % volumeTracks)) deallocate(self % volumeTracks)
    if(allocated(self % cellHit)) deallocate(self % cellHit)

  end subroutine kill

end module randomRayPhysicsPackage_class
