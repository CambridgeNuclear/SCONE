module randomRayPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,         only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,        only : FNV_1
  use exponentialRationalApprox, only : exponential
  use dictionary_class,          only : dictionary
  use rng_class,                 only : RNG
  use physicsPackage_inter,      only : physicsPackage

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : geometry, distCache
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_addGeom => addGeom, &
                                             gr_geomIdx  => geomIdx, gr_kill    => kill

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat, mm_matName => matName
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_getMatNames => getMatNames, &
                                             ndReg_kill        => kill
  use mgNeutronDatabase_class         only : mgNeutronDatabase
  use baseMgNeutronDatabase_class     only : baseMgNeutronDatabase
  use baseMgNeutronMaterial_class     only : baseMgNeutronMaterial

  implicit none
  private

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
    class(geometry), pointer              :: geom
    integer(shortInt)                     :: geomIdx     = 0
    real(defReal), dimension(3)           :: top         = ZERO
    real(defReal), dimension(3)           :: bottom      = ZERO
    type(RNG)                             :: rand
    integer(shortInt)                     :: timerMain   = 0
    class(baseMgNeutronDatabase), pointer :: mgData      => NULL
    integer(shortInt)                     :: nG          = 0
    integer(shortInt)                     :: nCells      = 0
    real(defReal)                         :: lengthPerIt = ZERO

    ! Settings
    real(defReal)      :: termination = ZERO
    real(defReal)      :: dead        = ZERO
    integer(shortInt)  :: pop         = 0
    integer(shortInt)  :: inactive    = 0
    integer(shortInt)  :: active      = 0
    !logical(defBool)   :: criterionA  = .TRUE.
    !logical(defBool)   :: criterionI  = .TRUE.
    !integer(shortInt)  :: accum       = 0
    !real(defReal)      :: std         = ZERO
    !real(defReal)      :: eps         = ZERO

    ! Results space
    real(defReal)                                :: keff
    real(defReal), dimension(2)                  :: keffScore
    real(defReal), dimension(:), allocatable     :: scalarFlux
    real(defReal), dimension(:), allocatable     :: prevFlux
    real(defReal), dimension(:,:), allocatable   :: fluxScores
    real(defReal), dimension(:), allocatable     :: source
    real(defReal), dimension(:), allocatable     :: volume
    integer(shortInt), dimension(:), allocatable :: cellHit

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: run
    procedure :: kill

    ! Private procedures
    procedure, private :: cycles
    procedure, private :: initialiseRay
    procedure, private :: moveRay
    procedure, private :: transportSweep
    procedure, private :: calculateSources
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
    character(nameLen)                            :: geomName, graphType
    character(100), parameter :: Here = 'init (randomRayPhysicsPackage_class.f90)'

    ! Load settings
    call dict % get(self % termination, 'termination')
    call dict % get(self % dead, 'dead')
    call dict % get(self % pop, 'pop')
    call dict % get(self % accum, 'accum')
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

    ! Check settings
    if (self % termination =< ZERO) call fatalError(Here, 'Ray termination distance (termination) is less than or equal to zero.')
    if (self % pop < 1) call fatalError('Must have 1 or more rays (pop).')
    !if (self % accum < 1) call fatalError('Must have 1 or more accumulation cycles (accum).')
    !if (self % std < ZERO) call fatalError('Convergence standard deviation (std) cannot be less than zero.')
    !if (self % eps < ZERO) call fatalError('RMS tolerance (eps) cannot be less than zero.')

    ! Dead length can be less than zero but will be reset to zero if so
    if (self % dead < ZERO) then
      self % dead = ZERO
      print *,'Warning: Dead length of rays (dead) was negative. This has been set to 0 instead.'
    end if

    ! Ensure termination length is longer than dead length
    if (self % termination <= self % dead) call fataError(Here,&
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

    ! Ensure that nuclear data is multi-group
    self % mgData => getNeutronMG()
    if (.NOT. associated(self % mgData)) call fatalError(Here,&
            'No MG nuclear database was constructed')

    ! Store number of energy groups for convenience
    self % nG = self % mgData % nGroups()

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'randomRayGeom'
    call gr_addGeom(geomName, tempDict)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Ensure that geometry graph is extended
    graphDict => tempDict % getDictPtr('graph')
    call graphDict % get(graphType,'type')
    if (graphType /= 'extended') call fatalError(Here,&
            'Geometry graph type must be "extended" for random ray calculations.')

    ! Get lower and upper corner of bounding box
    associate (aabb => self % geom % bounds())
      self % bottom = aabb(1:3)
      self % top    = aabb(4:6)
    end associate

    ! Store number of cells in geometry for convenience
    self % nCells = self % geom % numberOfCells()

    ! Allocate results space
    allocate(self % scalarFlux(self % nCells * self % nG))
    allocate(self % prevFlux(self % nCells * self % nG))
    allocate(self % fluxScore(self % nCells * self % nG, 2))
    allocate(self % source(self % nCells * self % nG))
    allocate(self % volume(self % nCells))
    allocate(self % cellHit(self % nCells))

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
    call self % cycles(self % rand)
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
  subroutine cycles(self, rand)
    class(randomRayPhysicsPackage), intent(inout) :: self
    class(RNG), intent(inout)                     :: rand
    type(ray)                                     :: r
    type(RNG), save                               :: pRNG
    real(defReal)                                 :: hitRate
    real(defReal)                       :: elapsed_T, end_T, T_toEnd, av_speed, cycle_T
   
    !$omp threadprivate(pRNG)

    !$omp parallel
    pRNG = rand
    !$omp end parallel

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Initialise fluxes 
    self % keff       = ONE
    self % keffVec(1) = self % keff
    self % scalarFlux = ONE
    self % prevFlux   = self % scalarFlux

    ! Initialise other results
    self % cellHit = 0
    self % volume  = ZERO 

    ! Stopping criterion is initially on flux convergence or number of convergence iterations.
    ! Will be replaced by RMS error in flux or number of scoring iterations afterwards.
    itInac = 0
    itAct  = 0
    isActive = .FALSE.

    ! Power iteration
    do while( stopping_criterion )
      
      if (isActive) then
        itAct = itAct + 1
      else
        itInac = itInac + 1
      end if
      it = itInac + itAct

      call self % calculateSource()
      
      !$omp parallel do schedule(dynamic) private(r)
      do i = 1, self % pop

        ! Set seed
        call pRNG % stride( (it-1) * self % pop + i )
        r % pRNG => pRNG 

        ! Set ray attributes
        call self % initialiseRay(r)

        ! Transport ray until termination criterion met
        call self % transportSweep(r)

      end do
      !$omp end parallel do

      ! Normalise flux estimate and combines with source
      call self % normaliseFluxAndVolume(it)

      ! Calculate new k
      call self % calculateKeff()

      ! Accumulate flux scores
      if (isActive) call self % accumulateFluxAndKeffScores()

      ! Calculate proportion of cells that were hit
      hitRate = sum(self % cellHit) / self % nCells
      self % cellHit = 0
      if (hitRate <= 0.99) print *,&
              "Warning: rays are missing ",numToChar(100*(ONE-hitRate)),"% of cells"

      ! Evaluate stopping criterion for active or inactive iterations
      if (isActive) then
        !if (self % criterionA) then
        !  rms = self % calculateRMS()
        !  stopping_criterion = rms > self % eps
        !else
        stopping_criterion = itAct < self % active
        !end if
      else
        if (self % criterionI) then
          isActive = ((it >= self % accum) .AND. (kStd < self % std))
        else
          isActive = (itInac >= self % inactive)
        end if
      end if

      ! Set fluxes 
      self % prevFlux = self % scalarFlux
      self % scalarFlux = ZERO

      ! Calculate times
      call timerStop(self % timerMain)
      cycle_T = timerTime(self % timerMain) - elapsed_T
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(self % N_cycles, defReal) * elapsed_T / gen
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      call printFishLineR(gen)
      print *
      print *, 'Cycle: ', numToChar(gen), ' of ', numToChar(self % N_cycles)
      print *, 'Pop: ', numToChar(self % pop)
      print '(A, ES12.5)', ' Av. Ray speed: [m/s]: ', av_speed
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))

    end do

    ! Finalise flux scores
    call self % finaliseFluxAndKeffScores(itAct)

      ! Process scores
      self % res(:, SCORE) = self % res(:, SCORE) / self % totDist
      self % res(:, CSUM)  = self % res(:, CSUM)  + self % res(:, SCORE)
      self % res(:, CSUM2) = self % res(:, CSUM2) + self % res(:, SCORE)**2
      self % res(:, SCORE) = ZERO
      self % totDist = ZERO


  end subroutine cycles

  !!
  !! Initialises rays: samples initial position and direction,
  !! performs the build operation, and sets the initial flux
  !!
  subroutine initialiseRay(self, r)
    class(randomRayPhysicsPackage), intent(in) :: self
    type(ray), intent(inout)                   :: r
    real(defReal)                              :: mu, phi
    real(defReal), dimension(3)                :: u, rand3, x
    integer(shortInt)                          :: i, g, cIdx, idx
    character(100), parameter :: Here = 'sampleInitial (randomRayPhysicsPackage_class.f90)'

    i = 0
    mu = TWO * r % rand % get() - ONE
    phi = TWO_PI * r % rand % get()
    u = rotateVector([ONE, ZERO, ZERO], mu, phi)

    rejection : do
      rand3(1) = pRNG % get()
      rand3(2) = pRNG % get()
      rand3(3) = pRNG % get()
      x = bottom + (self % top - self % bottom) * rand3

      ! Exit if point is inside the geometry
      call self % geom % whatIsAt(matIdx, uniqueId, r, u)
      if (matIdx /= OUTSIDE_MAT) exit rejection

      i = i + 1
      if (i > 5000) then
        call fatalError(Here, 'Infinite loop when searching ray start in the geometry.')
      end if
    end do rejection

    ! Place in the geometry & process the ray
    call r % build(x, u, self % nG)
    call self % geom % placeCoord(r % coords)

    ! Set angular flux to angle average of cell scalar flux
    cIdx = r % getCellIdx()
    do g = 1, self % nG
      idx = (cIdx - 1) * self % nG + g
      r % flux(g) = self % source(idx)
    end do

  end subroutine initialiseRay

  !!
  !! Geometry handling of the ray: moves the ray across a FSR,
  !! returning the distance travelled and applying boundary conditions
  !! where necessary.
  !!
  !! Accounts for biases that may be introduced by rays moving
  !! significantly past their dead length during an FSR crossing.
  !! If this occurs, the ray is moved back to the point where the
  !! dead length is exceeded and this is the distance which is 
  !! travelled during the operation.
  !!
  !! Reports the vacuum hit to be dealt with after attenuation
  !! has taken place.
  !!
  subroutine moveRay(self, r, length, hitVacuum, cache)
    class(randomRayPhysicsPackage), intent(in) :: self
    type(ray), intent(inout)                   :: r
    real(defReal), intent(out)                 :: length
    logical(defBool), intent(out)              :: hitVacuum
    type(distCache), intent(inout)             :: cache
    integer(shortInt)                          :: event

    ! Set maximum flight distance and ensure ray is active
    if (r % length >= self % dead) then
      length = INFINITY 
      r % isActive = .TRUE.
    else
      length = self % dead - r % length
    end if

    ! Perform the movement
    call self % geom % moveRay_withCache(r % coords, length, event, cache, hitVacuum)

    r % length = r % length + length

  end subroutine moveRay

  !!
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume
  !!
  subroutine transportSweep(self, r)
    class(randomRayPhysicsPackage), intent(inout) :: self
    type(ray), intent(inout)                      :: r
    integer(shortInt)                             :: matIdx, cIdx, idx, g
    real(defReal)                                 :: attenuate, length, delta
    real(defReal), dimension(:)                   :: SigmaT
    logical(defBool)                              :: hitVacuum
    type(distCache)                               :: cache

    do while (r % length < self % termination)

      ! Get material and cell the ray is moving through
      matIdx = r % getMatIdx()
      cIdx   = r % getCellIdx()
      mat    => self % mgData % getMaterial(matIdx)
      SigmaT = mat % getAllTotal(matIdx)
          
      ! Move ray
      call moveRay(r,length, hitVacuum, cache)

      do g = 1, self % nG
            
        ! Calculate delta
        idx = (cIdx - 1) * self % nG + g
        attenuate = exponential(SigmaT(g) * length)
        delta = (r % flux(g) - self % source(idx)) * attenuate

        ! Accumulate scalar flux
        if (r % isActive) then
          !$omp atomic
          self % scalarFlux(idx) = self % scalarFlux(idx) + delta
        end if

        ! Update flux
        r % flux(g) = r % flux(g) - delta

      end do

      ! Accumulatve cell volume estimates
      if (r % isActive) then
        self % cellHit(cIdx) = 1

        !$omp atomic
        self % volume(cIdx) = self % volume(cIdx) + length
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
    class(randomRayPHysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    real(defReal)                                 :: norm, normVol
    integer(shortInt)                             :: i, g, cIdx, matIdx, idx

    norm = ONE / self % lengthPerIt
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do
    do i = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(i) 
      mat    => self % mgData % getMaterial(matIdx)
      SigmaT =  mat % getAllTotalXS()
      
      ! Scale volume due to additional rays
      self % volume(i) = self % volume(i) * normVol

      do g = 1, self % nG

        idx = self % nG * (i - 1) + g
        if (self % volume(i) /= ZERO) then
          self % flux(idx) = flux(idx) * norm / ( SigmaT(g) * self % volume(i))
        end if
        self % flux(idx) = self % flux(idx) + self % source(idx)

      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFlux

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxAndKeffScores(self)
    class(randomRayPHysicsPackage), intent(inout) :: self
    integer(shortInt)                             :: i, g, idx

    !$omp parallel do
    do i = 1, self % nCells
      do g = 1, self % nG
        idx = self % nG * (i - 1) + g
        self % fluxScores(idx,1) = self % fluxScores(idx, 1) + self % scalarFlux(idx)
        self % fluxScores(idx,2) = self % fluxScores(idx, 2) + &
                self % scalarFlux(idx) * self % scalarFlux(idx)
      end do
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) + self % keff
    self % keffScore(2) = self % keffScore(2) + self % keff * self % keff

  end subroutine accumulateFluxScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxScores(self,it)
    class(randomRayPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                 :: it
    integer(shortInt)                             :: i, g, idx

    !$omp parallel do
    do i = 1, self % nCells
      do g = 1, self % nG
        idx = self % nG * (i - 1) + g
        self % fluxScores(idx,1) = self % fluxScores(idx, 1) / it
        self % fluxScores(idx,2) = self % fluxScores(idx, 2) / it
        self % fluxScores(idx,2) = &
                self % fluxScores(idx,1) * self % fluxScores(idx,1) - &
                self % fluxScores(idx,2)
      end do
    end do
    !$omp end parallel do

    self % keffScore(1) = self % keffScore(1) / it
    self % keffScore(2) = self % keffScore(2) / it
    self % keffScore(2) = &
            self % keffScore(1) * self % keffScore(1) - &
            self % keffScore(2)

  end subroutine finaliseFluxScores
  
  !!
  !! Calculate sources in all cells and energy groups
  !!
  subroutine calculateSource(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    real(defReal)                                 :: fission, scatter, ONE_KEFF
    real(defReal), dimension(:)                   :: nuFission, total, chi
    real(defReal), dimension(:,:)                 :: scatterMatrix
    integer(shortInt)                             :: i, idx, idx0, matIdx, g, gIn

    ONE_KEFF = ONE / self % keff

    !$omp parallel do
    do i = 1, self % nCells

      ! Identify material
      matIdx =  self % geom % geom % graph % getMatFromUID(i) 
      mat    => self % mgData % getMaterial(matIdx)

      ! Get XS info
      chi           = mat % getAllChi()
      nuFission     = mat % getAllNuFissionXS()
      scatterMatrix = mat % getScatterMatrix()
      total         = mat % getAllTotalXS()

      do g = 1, self % nG
      
        ! Source index
        idx = self % nG * (i - 1) + g

        ! Calculate fission and scattering source
        scatter = ZERO
        fission = ZERO

        ! Sum contributions from all energies
        do gIn = 1, self % nG

          idx0 = self % nG * (i - 1) + gIn
          fission = fission + self % scalarFlux(idx0) * nuSigmaF(gIn)
          scatter = scatter + self % scalarFlux(idx0) * scatterMatrix(g,gIn)
          
        end do
        
        self % source(idx) = chi(g) * fission * ONE_KEFF + scatter
        self % source(idx) = self % source(idx) / total(g)

      end do

    end do
    !$omp end parallel do

  end subroutine calculateSources

  !!
  !! Calculate keff
  !!
  subroutine calculateKeff(self)
    class(randomRayPhysicsPackage), intent(inout) :: self
    real(defReal)                                 :: fissionRate, prevFissionRate, vol
    real(defReal), dimension(:)                   :: nuFission
    integer(shortInt)                             :: i, idx, idx0, matIdx, g, gIn
    class(mgMaterial), pointer                    :: mat

    fissionRate     = ZERO
    prevFissionRate = ZERO
    !$omp parallel do reduction(+: fissionRate, prevFissionRate)
    do i = 1, self % nCells

      ! Identify material
      matIdx =  self % geom % geom % graph % getMatFromUID(i) 
      mat    => self % mgData % getMaterial(matIdx)

      ! Get XS info
      nuFission = mat % getAllNuFissionXS()
      vol = self % volume(i)

      if (vol == ZERO) continue

      do g = 1, self % nG
      
        ! Source index
        idx = self % nG * (i - 1) + g

        fissionRate     = fissionRate     + self % scalarFlux(idx) * nuFission(g) * vol
        prevFissionRate = prevFissionRate + self % prevFlux(idx) * nuFission(g) * vol

      end do

    end do
    !$omp end parallel do reduction

    ! Update k
    self % keff = self % keff * fissionRate / prevFissionRate

  end subroutine calculateKeff

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
    class(randomRayPhysicsPackage), intent(in) :: self
    real(defReal)                           :: mean, SD, var
    real(defReal)                           :: V, V_SD
    integer(shortInt)                       :: i


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
  !! Print settings of the random ray calculation
  !!
  !! Args:
  !!   None
  !!
  subroutine printSettings(self)
    class(randomRayPhysicsPackage), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ RANDOM RAY EIGENVALUE CALCULATION /\/\"
    if (self % criterionI) then
      print *, "Using convergence criterion for inactive"
      print *, "Accumulates k over ",numToChar(self % accum),&
              " cycles to reach a std. dev. of ",numToChar(self % std)
    else
      print *, "Using ",numToChar(self % inactive), " iterations for "&
              //"the inactive cycles"
    end if
    if (self % criterionA) then
      print *, "Using convergence criterion for active"
      print *, "Compares RMS within a tolerance of ",numToChar(self % eps)
    else
      print *, "Using ",numToChar(self % active), " iterations for "&
              //"the active cycles"
    end if
    print *, 
    print *, "Rays per cycle: ", numToChar(self % pop)
    print *, "Ray dead length: ",numToChar(self % dead)
    print *, "Ray termination length: ",numToChar(self % termination)
    print *, "Initial RNG Seed:   ", numToChar(self % rand % getSeed())
    print *,
    print *, "Number of cells in the geometry: ", self % nCells
    print *, "Number of energy groups: ", self % nG
    print *, repeat("<>", MAX_COL/2)

  end subroutine printSettings

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(randomRayPhysicsPackage), intent(inout) :: self

    ! Clean Nuclear Data & Geometries
    call gr_kill()
    call ndreg_kill()

    ! Clean contents
    self % geom    => null()
    self % geomIdx = 0
    self % timerMain = 0

    self % top       = ZERO
    self % bottom    = ZERO
    self % mgData    => NULL
    self % nG        = 0
    self % nCells    = 0

    self % termination = ZERO
    self % dead        = ZERO
    self % pop         = 0
    self % criterionA  = .TRUE.
    self % criterionI  = .TRUE.
    self % accum       = 0
    self % std         = ZERO
    self % eps         = ZERO
    self % inactive    = 0
    self % active      = 0

    self % keff        = ZERO
    if(allocated(keffVec)) deallocate(keffVec)
    if(allocated(scalarFlux)) deallocate(scalarFlux)
    if(allocated(prevFlux)) deallocate(prevFlux)
    if(allocated(source)) deallocate(source)
    if(allocated(volume)) deallocate(volume)
    if(allocated(cellHit)) deallocate(cellHit)

  end subroutine kill

end module randomRayPhysicsPackage_class
