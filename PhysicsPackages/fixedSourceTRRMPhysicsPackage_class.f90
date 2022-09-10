module fixedSourceTRRMPhysicsPackage_class

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

  ! Tally map for fission rate
  use tallyMap_inter,                 only : tallyMap
  use tallyMapFactory_func,           only : new_tallyMap

  ! Random ray
  use ray_class,                      only : ray

  ! Particle state (just for easier output)
  use particle_class,                 only : particleState

  implicit none
  private

  ! Parameter for when to skip a tiny volume
  real(defReal), parameter :: volume_tolerance = 1.0E-15

  !!
  !! Physics package to perform The Random Ray Method (TRRM) fixed source calculations
  !!
  !! Tracks rays across the geometry, attenuating their flux. After some dead length,
  !! rays begin scoring to estimates of the scalar flux and volume. Each ray has a
  !! uniform termination length, after which it is stopped and the next ray is tracked.
  !! Once all rays have been tracked, a cycle concludes and fluxes and sources are updated.
  !!
  !! Both inactive and active cycles occur, as in eigenvalue calculations. These can be terminated
  !! after a specified number of iterations or on reaching some chosen convergence
  !! criterion (though the latter hasn't been implemented yet).
  !!
  !! Calculates relative volume of diffrent materials in the problem by performing
  !! random ray tracing in the geometry. The volume is normalised such that the total domain
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
  !!     dead 10;              // Dead length where rays do not score to scalar fluxes
  !!     termination 100;      // Length a ray travels before it is terminated
  !!     rays 1000;            // Number of rays to sample per iteration
  !!     inactive 100;         // Number of convergence cycles (would use accum and std otherwise)
  !!     active 200;           // Number of scoring cycles (would use eps otherwise)
  !!     #seed 86868;#         // Optional RNG seed
  !!     #cache 1;#            // Optionally use distance caching to accelerate ray tracing
  !!     #fissionMap {<map>}#  // Optionally output fission rates according to a given map
  !!     #plot 1;#             // Optionally make VTK viewable plot of fluxes and uncertainties
  !!
  !!     #source {             // Fixed sources for named materials and their intensities n/cm3/s
  !!                           // Intensities are in each energy group, from 1 to G
  !!         material_name1 ( s_g1 s_g2 ... s_gG );
  !!         ...
  !!      } #
  !!
  !!     geometry {<Geometry definition>}
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
  !!   lengthPerIt -> Distance all rays travel in a single iteration - for convenience.
  !!
  !!   termination -> Distance a ray can travel before it is terminated
  !!   dead        -> Distance a ray must travel before it becomes active
  !!   pop         -> Number of rays to track per cycle
  !!   inactive    -> Number of inactive cycles to perform
  !!   active      -> Number of active cycles to perform
  !!   cache       -> Perform distance caching?
  !!   outputFile  -> Output file name
  !!   outputFormat-> Output file format
  !!   plotResults -> Plot results?
  !!   printFluxes -> Print fluxes?
  !!   printVolume -> Print volumes?
  !!   printCells  -> Print cell positions?
  !!   viz         -> Output visualiser
  !!   mapFission  -> Output fission rates across a given map?
  !!   resultsMap  -> The map across which to output fission rate results
  !!
  !!   scalarFlux   -> Array of scalar flux values of length = nG * nCells
  !!   prevFlux     -> Array of previous scalar flux values of length = nG * nCells
  !!   fluxScore    -> Array of scalar flux values and squared values to be reported 
  !!                   in results, dimension =  [nG * nCells, 2]
  !!   source       -> Array of neutron source values of length = nG * nCells
  !!   fixedSource  -> Array of fixed source values of length = nG * nCells
  !!   volume       -> Array of stochastically estimated cell volumes of length = nCells
  !!   cellHit      -> Array tracking whether given cells have been hit during tracking
  !!   cellFound    -> Array tracking whether a cell was ever found
  !!   cellPos      -> Array of cell positions, populated once they are found
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: fixedSourceTRRMPhysicsPackage
    private
    ! Components
    class(geometryStd), pointer           :: geom
    integer(shortInt)                     :: geomIdx     = 0
    real(defReal), dimension(3)           :: top         = ZERO
    real(defReal), dimension(3)           :: bottom      = ZERO
    class(baseMgNeutronDatabase), pointer :: mgData      => null()
    type(RNG)                             :: rand
    integer(shortInt)                     :: nG          = 0
    integer(shortInt)                     :: nCells      = 0
    real(defReal)                         :: lengthPerIt = ZERO

    ! Settings
    real(defReal)      :: termination = ZERO
    real(defReal)      :: dead        = ZERO
    integer(shortInt)  :: pop         = 0
    integer(shortInt)  :: inactive    = 0
    integer(shortInt)  :: active      = 0
    logical(defBool)   :: cache       = .FALSE.
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    logical(defBool)   :: plotResults = .FALSE.
    logical(defBool)   :: printFlux   = .FALSE.
    logical(defBool)   :: printVolume = .FALSE.
    logical(defBool)   :: printCells  = .FALSE.
    type(visualiser)   :: viz
    logical(defBool)   :: mapFission  = .FALSE.
    class(tallyMap), allocatable :: resultsMap

    ! Results space
    real(defReal), dimension(:), allocatable     :: scalarFlux
    real(defReal), dimension(:), allocatable     :: prevFlux
    real(defReal), dimension(:,:), allocatable   :: fluxScores
    real(defReal), dimension(:), allocatable     :: source
    real(defReal), dimension(:), allocatable     :: fixedSource
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
    procedure, private :: initialiseSource
    procedure, private :: initialiseRay
    procedure, private :: transportSweep
    procedure, private :: sourceUpdateKernel
    procedure, private :: normaliseFluxAndVolume
    procedure, private :: resetFluxes
    procedure, private :: accumulateFluxScores
    procedure, private :: finaliseFluxScores
    procedure, private :: printResults
    procedure, private :: printSettings

  end type fixedSourceTRRMPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)                    :: dict
    integer(shortInt)                                   :: seed_temp
    integer(longInt)                                    :: seed
    character(10)                                       :: time
    character(8)                                        :: date
    character(:),allocatable                            :: string
    class(dictionary),pointer                           :: tempDict, graphDict, sourceDict
    class(mgNeutronDatabase),pointer                    :: db
    character(nameLen)                                  :: geomName, graphType, nucData
    class(geometry), pointer                            :: geom
    type(outputFile)                                    :: test_out
    character(100), parameter :: Here = 'init (fixedSourceTRRMPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)
    
    ! Load settings
    call dict % get( nucData, 'XSdata')
    call dict % get(self % termination, 'termination')
    call dict % get(self % dead, 'dead')
    call dict % get(self % pop, 'pop')
    call dict % get(self % active, 'active')
    call dict % get(self % inactive, 'inactive')
    
    ! Perform distance caching?
    call dict % getOrDefault(self % cache, 'cache', .FALSE.)

    ! Print fluxes?
    call dict % getOrDefault(self % printFlux, 'printFlux', .FALSE.)

    ! Print volumes?
    call dict % getOrDefault(self % printVolume, 'printVolume', .FALSE.)

    ! Print cell positions?
    call dict % getOrDefault(self % printCells, 'printCells', .FALSE.)

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be cought early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Check settings
    if (self % termination <= ZERO) call fatalError(Here, 'Ray termination distance (termination) is less than or equal to zero.')
    if (self % pop < 1) call fatalError(Here, 'Must have 1 or more rays (pop).')

    ! Dead length can be less than zero but will be reset to zero if so
    if (self % dead < ZERO) then
      self % dead = ZERO
      print *,'Warning: Dead length of rays (dead) was negative. This has been set to 0 instead.'
    end if

    ! Ensure termination length is longer than dead length
    if (self % termination <= self % dead) call fatalError(Here,&
            'Ray termination length must be greater than ray dead length')

    ! Check whether there is a map for outputting fission rates
    ! If so, read and initialise the map to be used
    if (dict % isPresent('fissionMap')) then
      self % mapFission = .TRUE.
      tempDict => dict % getDictPtr('fissionMap')
      call new_tallyMap(self % resultsMap, tempDict)
    else
      self % mapFission = .FALSE.
    end if

    ! Register timer
    self % timerMain = registerTimer('simulationTime')
    self % timerTransport = registerTimer('transportTime')

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

    ! Activate nuclear data
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
    allocate(self % fixedSource(self % nCells * self % nG))
    allocate(self % volume(self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(self % nCells, 3))
    
    ! Read and initialise the fixed sources
    sourceDict => dict % getDictPtr('source')
    call self % initialiseSource(sourceDict)

    ! Set active length traveled per iteration
    self % lengthPerIt = (self % termination - self % dead) * self % pop

  end subroutine init

  !!
  !! Initialises the fixed source to be used in the simulation
  !! Takes a dictionary containing names of materials in the geometry and
  !! source strengths in each energy group and places these in the appropriate
  !! elements of the fixed source vector
  !!
  subroutine initialiseSource(self, dict)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)                    :: dict
    character(nameLen),dimension(:), allocatable        :: names
    real(defReal), dimension(:), allocatable            :: sourceStrength
    integer(shortInt)                                   :: cIdx, i
    integer(shortInt), save                             :: g, matIdx, idx
    logical(defBool)                                    :: found
    character(nameLen)                                  :: sourceName 
    character(nameLen), save                            :: localName
    character(100), parameter :: Here = 'initialiseSource (fixedSourceTRRMPhysicsPackage_class.f90)'
    !$omp threadprivate(matIdx, localName, idx, g)

    self % fixedSource = ZERO

    call dict % keys(names)

    ! Cycle through entries of the dictionary
    do i = 1, size(names)

      sourceName = names(i)
      call dict % get(sourceStrength, sourceName)

      ! Ensure correct number of energy groups
      if (size(sourceStrength) /= self % nG) call fatalError(Here,'Source '//sourceName//&
              ' has '//numToChar(size(sourceStrength))//' groups rather than '//numToChar(self % nG))
      
      ! Make sure that the source corresponds to a material present in the geometry
      found = .FALSE.
      !$omp parallel do schedule(static) 
      do cIdx = 1, self % nCells

        matIdx    = self % geom % geom % graph % getMatFromUID(cIdx)
        localName = mm_matName(matIdx)

        if (localName == sourceName) then

          found = .TRUE.
          do g = 1, self % nG
            
            idx = (cIdx - 1) * self % nG + g
            self % fixedSource(idx) = sourceStrength(g)

          end do

        end if

      end do
      !$omp end parallel do

      if (.NOT. found) call fatalError(Here,'The source '//trim(sourceName)//' does not correspond to '//&
              'any material found in the geometry.')

    end do

  end subroutine initialiseSource

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self

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
  subroutine cycles(self)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self
    type(ray), save                                     :: r
    type(RNG), target, save                             :: pRNG
    real(defReal)                                       :: hitRate
    real(defReal)                                       :: elapsed_T, end_T, T_toEnd, transport_T
    logical(defBool)                                    :: stoppingCriterion, isActive
    integer(shortInt)                                   :: i, itInac, itAct, it
    integer(longInt), save                              :: ints
    integer(longInt)                                    :: intersections
    !$omp threadprivate(pRNG, r, ints)

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Initialise fluxes 
    self % scalarFlux = ZERO
    self % prevFlux   = ZERO
    self % fluxScores = ZERO
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

    ! Source iteration
    do while( stoppingCriterion )
      
      if (isActive) then
        itAct = itAct + 1
      else
        itInac = itInac + 1
      end if
      it = itInac + itAct

      !$omp parallel do schedule(static)
      do i = 1, self % nCells
        call self % sourceUpdateKernel(i)
      end do
      !$omp end parallel do
    
      ! Reset and start transport timer
      call timerReset(self % timerTransport)
      call timerStart(self % timerTransport)
      intersections = 0
      
      !$omp parallel do schedule(dynamic) reduction(+: intersections)
      do i = 1, self % pop

        ! Set seed
        pRNG = self % rand
        call pRNG % stride(i)
        r % pRNG => pRNG 

        ! Set ray attributes
        call self % initialiseRay(r)

        ! Transport ray until termination criterion met
        call self % transportSweep(r,ints)
        intersections = intersections + ints

      end do
      !$omp end parallel do
      
      call timerStop(self % timerTransport)

      ! Update RNG 
      call self % rand % stride(self % pop + 1)

      ! Normalise flux estimate and combines with source
      call self % normaliseFluxAndVolume(it)

      ! Accumulate flux scores
      if (isActive) call self % accumulateFluxScores()

      ! Calculate proportion of cells that were hit
      hitRate = real(sum(self % cellHit),defReal) / self % nCells
      self % cellHit = 0

      ! Evaluate stopping criterion for active or inactive iterations
      if (isActive) then
        stoppingCriterion = (itAct < self % active)
      else
        isActive = (itInac >= self % inactive)
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
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      print *, 'Time per integration (ns): ', &
              trim(numToChar(transport_T*10**9/(self % nG * intersections)))

    end do

    ! Finalise flux scores
    call self % finaliseFluxScores(itAct)

  end subroutine cycles

  !!
  !! Initialises rays: samples initial position and direction,
  !! and performs the build operation
  !!
  subroutine initialiseRay(self, r)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self
    type(ray), intent(inout)                            :: r
    real(defReal)                                       :: mu, phi
    real(defReal), dimension(3)                         :: u, rand3, x
    integer(shortInt)                                   :: i, matIdx, cIdx
    character(100), parameter :: Here = 'initialiseRay (fixedSourceTRRMPhysicsPackage_class.f90)'

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
        call fatalError(Here, 'Infinite loop when searching for ray start in the geometry.')
      end if
    end do rejection

    ! Place in the geometry & process the ray
    call r % build(x, u)
    call self % geom % placeCoord(r % coords)

    if (.NOT. self % cellFound(cIdx)) then
      !$omp critical 
      self % cellFound(cIdx) = .TRUE.
      self % cellPos(cIdx,:) = x
      !$omp end critical
    end if

  end subroutine initialiseRay

  !!
  !! Moves ray through geometry, updating angular flux and
  !! scoring scalar flux and volume.
  !! Records the number of integrations/ray movements.
  !!
  subroutine transportSweep(self, r, ints)
    class(fixedSourceTRRMPhysicsPackage), target, intent(inout) :: self
    type(ray), intent(inout)                                    :: r
    integer(longInt), intent(out)                               :: ints
    integer(shortInt)                                           :: matIdx, g, cIdx, idx, event, matIdx0, baseIdx
    real(defReal)                                               :: totalLength, length
    logical(defBool)                                            :: activeRay, hitVacuum
    type(distCache)                                             :: cache
    class(baseMgNeutronMaterial), pointer                       :: mat
    class(materialHandle), pointer                              :: matPtr
    real(defReal), dimension(self % nG)                         :: attenuate, delta, fluxVec
    real(defReal), pointer, dimension(:)                        :: scalarVec, sourceVec, totVec
    
    ! Set initial angular flux to angle average of cell source
    cIdx = r % coords % uniqueID
    do g = 1, self % nG
      idx = (cIdx - 1) * self % nG + g
      fluxVec(g) = self % source(idx)
    end do
    
    ints = 0
    matIdx0 = 0
    totalLength = ZERO
    activeRay = .FALSE.
    do while (totalLength < self % termination)

      ! Get material and cell the ray is moving through
      matIdx  = r % coords % matIdx
      cIdx    = r % coords % uniqueID
      if (matIdx0 /= matIdx) then
        matPtr  => self % mgData % getMaterial(matIdx)
        mat     => baseMgNeutronMaterial_CptrCast(matPtr)
        matIdx0 = matIdx
        
        ! Cache total cross section
        totVec => mat % getTotalPtr()
      end if

      ! Remember new cell positions
      if (.NOT. self % cellFound(cIdx)) then
        !$omp critical 
        self % cellFound(cIdx) = .TRUE.
        self % cellPos(cIdx,:) = r % rGlobal()
        !$omp end critical
      end if
          
      ! Set maximum flight distance and ensure ray is active
      if (totalLength >= self % dead) then
        length = self % termination - totalLength 
        activeRay = .TRUE.
      else
        length = self % dead - totalLength
      end if

      ! Move ray
      ! Use distance caching or standard ray tracing
      ! Distance caching seems a little bit more unstable
      ! due to FP error accumulation, but is faster.
      ! This can be fixed by resetting the cache after X number
      ! of distance calculations.
      if (self % cache) then
        if (mod(ints,100_longInt) == 0)  cache % lvl = 0
        call self % geom % moveRay_withCache(r % coords, length, event, cache, hitVacuum)
      else
        call self % geom % moveRay_noCache(r % coords, length, event, hitVacuum)
      end if
      totalLength = totalLength + length

      ints = ints + 1

      baseIdx = (cIdx - 1) * self % nG
      sourceVec => self % source(baseIdx + 1 : baseIdx + self % nG)
      scalarVec => self % scalarFlux(baseIdx + 1 : baseIdx + self % nG)

      ! Accumulate to scalar flux
      if (activeRay) then
        !$omp simd
        do g = 1, self % nG
          attenuate(g) = exponential(totVec(g) * length)
          delta(g) = (fluxVec(g) - sourceVec(g)) * attenuate(g)
          fluxVec(g) = fluxVec(g) - delta(g)
        end do

        ! Accumulate scalar flux
        !$omp simd
        do g = 1, self % nG
          !$omp atomic
          scalarVec(g) = scalarVec(g) + delta(g) 
        end do

        ! Accumulate cell volume estimates
        self % cellHit(cIdx) = 1
        !$omp atomic
        self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length

      ! Don't accumulate to scalar flux
      else
        !$omp simd
        do g = 1, self % nG
          attenuate(g) = exponential(totVec(g) * length)
          delta(g) = (fluxVec(g) - sourceVec(g)) * attenuate(g)
          fluxVec(g) = fluxVec(g) - delta(g)
        end do
      end if

      ! Check for a vacuum hit
      if (hitVacuum) then
        !$omp simd
        do g = 1, self % nG
          fluxVec(g) = ZERO
        end do
      end if

    end do

  end subroutine transportSweep

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                       :: it
    real(defReal)                                       :: norm, normVol
    real(defReal), save                                 :: total
    integer(shortInt), save                             :: g, matIdx, idx
    integer(shortInt)                                   :: cIdx
    class(baseMgNeutronMaterial), pointer, save         :: mat
    class(materialHandle), pointer, save                :: matPtr
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

        if (self % volume(cIdx) > ZERO) then
          self % scalarFlux(idx) = self % scalarFlux(idx) * norm / ( total * self % volume(cIdx))
        end if
        self % scalarFlux(idx) = self % scalarFlux(idx) + self % source(idx)

      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolume
  
  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernel(self, cIdx)
    class(fixedSourceTRRMPhysicsPackage), target, intent(inout) :: self
    integer(shortInt), intent(in)                               :: cIdx
    real(defReal)                                               :: scatter, fission
    real(defReal), dimension(:), pointer                        :: nuFission, total, chi 
    real(defReal), dimension(:,:), pointer                      :: scatterXS
    integer(shortInt)                                           :: matIdx, g, gIn 
    integer(shortInt)                                           :: baseIdx, idx
    class(baseMgNeutronMaterial), pointer                       :: mat
    class(materialHandle), pointer                              :: matPtr
    logical(defBool)                                            :: isFiss
    real(defReal), pointer, dimension(:)                        :: fluxVec

    ! Identify material
    matIdx  =  self % geom % geom % graph % getMatFromUID(cIdx) 
    matPtr  => self % mgData % getMaterial(matIdx)
    mat     => baseMgNeutronMaterial_CptrCast(matPtr)

    ! Obtain XSs
    scatterXS =>  mat % getScatterPtr()
    total     =>  mat % getTotalPtr()
    isFiss    = mat % isFissile()
    if (isFiss) then
      nuFission =>  mat % getNuFissionPtr()
      chi       =>  mat % getChiPtr()
    end if

    baseIdx = self % ng * (cIdx - 1)
    fluxVec => self % prevFlux(baseIdx+1:baseIdx+self % nG)

    ! There is a scattering and fission source
    if (isFiss) then
      do g = 1, self % nG

        ! Calculate fission and scattering source
        scatter = ZERO
        fission = ZERO

        ! Sum contributions from all energies
        !$omp simd reduction(+:fission, scatter)
        do gIn = 1, self % nG
          fission = fission + fluxVec(gIn) * nuFission(gIn)
          scatter = scatter + fluxVec(gIn) * scatterXS(g,gIn)
        end do

        ! Output index
        idx = baseIdx + g

        self % source(idx) = chi(g) * fission + scatter + self % fixedSource(idx)
        self % source(idx) = self % source(idx) / total(g)

      end do

    ! There is only a scattering source
    else
      do g = 1, self % nG

        ! Calculate scattering source
        scatter = ZERO

        ! Sum contributions from all energies
        !$omp simd reduction(+:scatter)
        do gIn = 1, self % nG
          scatter = scatter + fluxVec(gIn) * scatterXS(g,gIn)
        end do

        ! Output index
        idx = baseIdx + g

        self % source(idx) = (self % fixedSource(idx) + scatter) / total(g)

      end do

    end if

  end subroutine sourceUpdateKernel

  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !!
  subroutine resetFluxes(self)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt)                                   :: idx

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
  subroutine accumulateFluxScores(self)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self
    real(defReal), save                                 :: flux
    integer(shortInt)                                   :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = self % scalarFlux(idx)
      self % fluxScores(idx,1) = self % fluxScores(idx, 1) + flux
      self % fluxScores(idx,2) = self % fluxScores(idx, 2) + flux*flux
    end do
    !$omp end parallel do

  end subroutine accumulateFluxScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxScores(self,it)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                       :: it
    integer(shortInt)                                   :: idx
    real(defReal)                                       :: N1, Nm1

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

  end subroutine finaliseFluxScores
  
  !!
  !! Output calculation results to a file
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self
    type(outputFile)                                    :: out
    character(nameLen)                                  :: name
    integer(shortInt)                                   :: cIdx, g1
    integer(shortInt), save                             :: idx, matIdx, i, g
    real(defReal), save                                 :: vol, SigmaF
    type(particleState), save                           :: s
    integer(shortInt),dimension(:),allocatable          :: resArrayShape
    real(defReal), dimension(:), allocatable            :: groupFlux, fiss, fissSTD
    class(baseMgNeutronMaterial), pointer, save         :: mat
    class(materialHandle), pointer, save                :: matPtr
    !$omp threadprivate(idx, matIdx, i, mat, matPtr, vol, s, SigmaF, g)

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
    
    name = 'Clock_Time'
    call out % printValue(timerTime(self % timerMain),name)
    
    ! Print cell volumes
    if (self % printVolume) then
      name = 'volume'
      call out % startBlock(name)
      resArrayShape = [size(self % volume)]
      call out % startArray(name, resArrayShape)
      do cIdx = 1, self % nCells
        call out % addResult(self % volume(cIdx), ZERO)
      end do
      call out % endArray()
      call out % endBlock()
    end if

    ! Print cell positions
    if (self % printCells) then
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
    end if

    ! Print fluxes
    if (self % printFlux) then
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
    end if

    ! Send fission rates to map output
    if (self % mapFission) then
      resArrayShape = self % resultsMap % binArrayShape()
      allocate(fiss(self % resultsMap % bins(0)))
      allocate(fissSTD(self % resultsMap % bins(0)))
      fiss    = ZERO
      fissSTD = ZERO

      ! Find whether cells are in map and sum their contributions
      !$omp parallel do reduction(+: fiss, fissSTD)
      do cIdx = 1, self % nCells
        
        ! Identify material
        matIdx =  self % geom % geom % graph % getMatFromUID(cIdx) 
        matPtr => self % mgData % getMaterial(matIdx)
        mat    => baseMgNeutronMaterial_CptrCast(matPtr)
        vol    =  self % volume(cIdx)

        if (vol < volume_tolerance) cycle

        ! Fudge a particle state to search tally map
        s % r = self % cellPos(cIdx,:)
        i = self % resultsMap % map(s)

        if (i > 0) then
          do g = 1, self % nG
            SigmaF = mat % getFissionXS(g, self % rand)
            idx = (cIdx - 1)* self % nG + g
            fiss(i) = fiss(i) + vol * self % fluxScores(idx,1) * SigmaF
            ! Is this correct? Also neglects uncertainty in volume - assumed small.
            fissSTD(i) = fissSTD(i) + &
                    vol * self % fluxScores(idx,2)*self % fluxScores(idx,2) * SigmaF
          end do
        end if

      end do
      !$omp end parallel do

      do i = 1,size(fissSTD)
        fissSTD(i) = sqrt(fissSTD(i))
        if (fiss(i) > 0) fissSTD(i) = fissSTD(i) / fiss(i)
      end do

      name = 'fissionRate'
      call out % startBlock(name)
      call out % startArray(name, resArrayShape)
      ! Add all map elements to results
      do idx = 1, self % resultsMap % bins(0)
        call out % addResult(fiss(idx), fissSTD(idx))
      end do
      call out % endArray()
      call out % endBlock()
      
      ! Output tally map
      call self % resultsMap % print(out)
    end if

    call out % writeToFile(self % outputFile)

    ! Send all fluxes and stds to VTK
    if (self % plotResults) then
      allocate(groupFlux(self % nCells))
      do g1 = 1, self % nG
        name = 'flux_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % fluxScores(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'std_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % fluxScores(idx,2) / self % fluxScores(idx,1)
        end do
        !$omp end parallel do
        call self % viz % addVTKData(groupFlux,name)
      end do
      do g1 = 1, self % nG
        name = 'source_g'//numToChar(g1)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          idx = (cIdx - 1)* self % nG + g1
          groupFlux(cIdx) = self % source(idx)
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
    class(fixedSourceTRRMPhysicsPackage), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ RANDOM RAY FIXED SOURCE CALCULATION /\/\"
    print *, "Using "//numToChar(self % inactive)// " iterations for "&
              //"the inactive cycles"
    print *, "Using "//numToChar(self % active)// " iterations for "&
              //"the active cycles"
    print * 
    print *, "Rays per cycle: "// numToChar(self % pop)
    print *, "Ray dead length: "//numToChar(self % dead)
    print *, "Ray termination length: "//numToChar(self % termination)
    print *, "Initial RNG Seed:   "// numToChar(self % rand % getSeed())
    print *
    print *, "Number of cells in the geometry: "// numToChar(self % nCells)
    print *, "Number of energy groups: "// numToChar(self % nG)
    if (self % cache) print *, "Accelerated with distance caching"
    print *, repeat("<>", MAX_COL/2)

  end subroutine printSettings

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(fixedSourceTRRMPhysicsPackage), intent(inout) :: self

    ! Clean Nuclear Data, Geometry and visualisation
    call gr_kill()
    call ndreg_kill()
    call self % viz % kill()

    ! Clean contents
    self % geom    => null()
    self % geomIdx = 0
    self % timerMain = 0
    self % timerTransport = 0

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
    self % cache       = .FALSE.
    self % mapFission  = .FALSE.
    self % plotResults = .FALSE.
    self % printFlux   = .FALSE.
    self % printVolume = .FALSE.
    self % printCells  = .FALSE.

    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % source)) deallocate(self % source)
    if(allocated(self % fixedSource)) deallocate(self % fixedSource)
    if(allocated(self % volume)) deallocate(self % volume)
    if(allocated(self % volumeTracks)) deallocate(self % volumeTracks)
    if(allocated(self % cellHit)) deallocate(self % cellHit)
    if(allocated(self % resultsMap)) then
      call self % resultsMap % kill()
      deallocate(self % resultsMap)
    end if

  end subroutine kill

end module fixedSourceTRRMPhysicsPackage_class
