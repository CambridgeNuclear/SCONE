module randomRayPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar, printFishLineR
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile
  use rng_class,                      only : RNG
  use physicsPackage_inter,           only : physicsPackage

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryStd_class,              only : geometryStd
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx, &
                                             gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use geometryFactory_func,           only : new_geometry

  ! Nuclear Data
  use nuclearDataReg_mod,             only : ndReg_init         => init, &
                                             ndReg_getMatNames  => getMatNames, &
                                             ndReg_activate     => activate, &
                                             ndReg_kill         => kill, &
                                             ndReg_getNeutronMG => getNeutronMG
  use mgNeutronDatabase_inter,        only : mgNeutronDatabase
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use baseMgNeutronMaterial_class,    only : baseMgNeutronMaterial, baseMgNeutronMaterial_CptrCast
  
  ! Visualisation
  use visualiser_class,               only : visualiser

  ! Tally map for fission rate
  use tallyAdmin_class,               only : tallyAdmin

  ! Random ray specific modules
  use dataRR_class,                   only : dataRR
  use arraysRR_class,                 only : arraysRR
  use rayHandling_func,               only : transportSweep, initialiseRay

  ! Random ray - or a standard particle
  use particle_class,                 only : ray => particle
  use particleDungeon_class,          only : particleDungeon

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
  !! criterion (though the latter hasn't been implemented yet).
  !!
  !! Calculates relative volume of different materials in the problem by performing
  !! random ray tracing in the geometry. The volume is normalised such that the total domain
  !! volume is 1.0.
  !!
  !! IMPORTANT N.B.: Geometry type must be extended! Won't run if shrunk.
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
  !!     inactive 100;         // Number of convergence cycles 
  !!     active 200;           // Number of scoring cycles 
  !!     #seed 86868;#         // Optional RNG seed
  !!     #cache 1;#            // Optionally use distance caching to accelerate ray tracing
  !!     #fissionMap {<map>}#  // Optionally output fission rates according to a given map
  !!     #fluxMap {<map>}#     // Optionally output one-group fluxes according to a given map
  !!     #plot 1;#             // Optionally make VTK viewable plot of fluxes and uncertainties
  !!     #rho 0;#              // Optional stabilisation for negative in-group scattering XSs
  !!     #lin 0;#              // Optionally use linear (rather than flat) sources
  !!     #2d 0;#               // Optional input to stablise linear sources in 2D problems
  !!     #volPolicy 1;#        // Optional input to specify how volumes should be handled
  !!     #missPolicy 1;#       // Optional input to specify how misses should be handled
  !!
  !!     geometry {<Geometry Definition>}
  !!     nuclearData {<Nuclear data definition>}
  !!     tally {<Tally definition>}
  !!   }
  !!
  !! Private Members
  !!   geom        -> Pointer to the geometry.
  !!   geomIdx     -> Index of the geometry in geometry Registry.
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
  !!   cache       -> Logical check whether to use distance caching
  !!   outputFile  -> Output file name
  !!   outputFormat-> Output file format
  !!   plotResults -> Plot results?
  !!   viz         -> Output visualiser
  !!   tally       -> Tally admin for outputting results
  !!
  !!   keff        -> Estimated value of keff
  !!   keffScore   -> Vector holding cumulative keff score and keff^2 score
  !!
  !!   tally       -> Tally admin to output results
  !!
  !!   intersectionsTotal -> Total number of ray traces for the calculation
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: randomRayPhysicsPackage
    private
    ! Components
    class(geometryStd), pointer           :: geom
    integer(shortInt)                     :: geomIdx  = 0
    type(RNG)                             :: rand
    type(arraysRR)                        :: arrays
    type(dataRR)                          :: XSData
    class(baseMgNeutronDatabase), pointer :: mgData   => null()
    integer(shortInt)                     :: nG       = 0
    integer(shortInt)                     :: nCells   = 0
    type(tallyAdmin),pointer              :: tally    => null()

    ! Settings
    real(defReal)      :: termination = ZERO
    real(defReal)      :: dead        = ZERO
    integer(shortInt)  :: pop         = 0
    integer(shortInt)  :: inactive    = 0
    integer(shortInt)  :: active      = 0
    logical(defBool)   :: cache       = .false.
    real(defReal)      :: rho         = ZERO
    logical(defBool)   :: lin         = .false.
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    logical(defBool)   :: plotResults = .false.
    logical(defBool)   :: printFlux   = .false.
    logical(defBool)   :: printVolume = .false.
    logical(defBool)   :: printCells  = .false.
    type(visualiser)   :: viz

    ! Results space
    ! keffScore is public for integration testing
    real(defReal)                        :: keff = ONE
    real(defReal), dimension(2), public  :: keffScore = ZERO
    integer(longInt)                     :: intersectionsTotal = 0

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
    procedure, private :: printResults
    procedure, private :: printSettings

  end type randomRayPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self, dict, loud)
    class(randomRayPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)              :: dict
    logical(defBool), intent(in), optional        :: loud
    integer(shortInt)                             :: seed_temp, volP, missP
    integer(longInt)                              :: seed
    character(10)                                 :: time
    character(8)                                  :: date
    character(:),allocatable                      :: string
    class(dictionary),pointer                     :: tempDict, graphDict
    class(mgNeutronDatabase),pointer              :: db
    character(nameLen)                            :: geomName, graphType, nucData
    class(geometry), pointer                      :: geom
    type(outputFile)                              :: test_out
    logical(defBool)                              :: set2D
    character(100), parameter :: Here = 'init (randomRayPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    if (present(loud)) then
      self % loud = loud
    else
      self % loud = .true.
    end if
    
    ! Load settings
    call dict % get( nucData, 'XSdata')
    call dict % get(self % termination, 'termination')
    call dict % get(self % dead, 'dead')
    call dict % get(self % pop, 'pop')
    call dict % get(self % active, 'active')
    call dict % get(self % inactive, 'inactive')
    call dict % getOrDefault(self % keff, 'keff', ONE)
    
    call dict % getOrDefault(volP, 'volPolicy', 1)
    call dict % getOrDefault(missP, 'missPolicy', 1)
    
    ! Perform distance caching?
    call dict % getOrDefault(self % cache, 'cache', .false.)
    
    ! Stabilisation factor for negative in-group scattering
    call dict % getOrDefault(self % rho, 'rho', ZERO)
    
    ! Use linear sources?
    call dict % getOrDefault(self % lin, 'lin', .false.)

    ! Is the problem 2D?
    call dict % getOrDefault(set2D, '2d', .false.)
    
    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be cought early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Check settings
    if (self % termination <= ZERO) call fatalError(Here, &
            'Ray termination distance (termination) is less than or equal to zero.')
    if (self % pop < 1) call fatalError(Here, 'Must have 1 or more rays (pop).')
    if (self % dead < ZERO) call fatalError(Here, 'Dead length must be positive.')
    if (self % termination <= self % dead) call fatalError(Here,&
            'Ray termination length must be greater than ray dead length')

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
    call new_geometry(tempDict, geomName, silent = .not. self % loud)
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
    
    ! Initialise tally Admin
    if (dict % isPresent('tally')) then
      tempDict => dict % getDictPtr('tally')
      allocate(self % tally)
      call self % tally % init(tempDict)
    end if

    ! Activatee nuclear data
    call ndReg_activate(P_NEUTRON_MG, nucData, self % geom % activeMats(), silent = .not. self % loud)

    ! Ensure that nuclear data is multi-group
    db => ndReg_getNeutronMG()
    if (.not. associated(db)) call fatalError(Here,&
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

    ! Call visualisation
    if (dict % isPresent('viz')) then
      if (self % loud) print *, "Initialising visualiser"
      tempDict => dict % getDictPtr('viz')
      call self % viz % init(geom, tempDict)
      if (self % loud) print *, "Constructing visualisation"
      call self % viz % makeViz()
      call self % viz % kill()
    endif
    
    ! Check for results plotting and initialise VTK
    call dict % getOrDefault(self % plotResults,'plot',.false.)
    if (self % plotResults) then
      ! Initialise a visualiser to be used when results are available
      if (self % loud) print *, "Initialising results visualiser"
      tempDict => dict % getDictPtr('viz')
      call self % viz % init(geom, tempDict)
      if (self % loud) print *, "Constructing geometry visualisation"
      call self % viz % initVTK()
    end if

    ! Store number of cells in geometry for convenience
    self % nCells = self % geom % numberOfCells()

    ! Initialise RR arrays and nuclear data
    call self % arrays % init(self % mgData, self % geom, &
            self % pop * (self % termination - self % dead), self % rho, self % lin, &
            .false., self % loud, volPolicy = volP, missPolicy = missP, set2D = set2D)
    
  end subroutine init

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(randomRayPhysicsPackage), intent(inout) :: self

    if (self % loud) call self % printSettings()
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
    class(randomRayPhysicsPackage), target, intent(inout) :: self
    type(ray), save            :: r
    type(RNG), target, save    :: pRNG
    real(defReal)              :: hitRate 
    real(defReal)              :: ONE_KEFF, elapsed_T, end_T, T_toEnd, transport_T, &
                                  N1, Nm1
    logical(defBool)           :: keepRunning, isActive
    integer(shortInt)          :: i, itInac, itAct, it
    integer(longInt), save     :: ints
    integer(longInt)           :: intersections
    class(arraysRR), pointer   :: arrayPtr
    type(particleDungeon)      :: dummyDungeon
    !$omp threadprivate(pRNG, r, ints)

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    arrayPtr => self % arrays

    ! Stopping criterion is on number of convergence iterations.
    ! TODO: Make this on, e.g., entropy during inactive, followed by stochastic noise during active!
    itInac = 0
    itAct  = 0
    isActive = .false.
    keepRunning = .true.
    
    ! Power iteration
    do while( keepRunning )
      
      if (isActive) then
        itAct = itAct + 1
      else
        itInac = itInac + 1
      end if
      it = itInac + itAct
      
      ONE_KEFF = ONE / self % keff
      call arrayPtr % updateSource(ONE_KEFF)

      ! Reset and start transport timer
      call timerReset(self % timerTransport)
      call timerStart(self % timerTransport)
      intersections = 0
      
      !$omp parallel do schedule(dynamic) reduction(+:intersections)
      do i = 1, self % pop

        ! Set seed
        pRNG = self % rand
        call pRNG % stride(i)
        r % pRNG => pRNG 

        ! Set ray attributes
        call initialiseRay(r, arrayPtr)

        ! Transport ray until termination criterion met
        call transportSweep(r, ints, self % nG, self % cache, self % dead, &
                self % termination, arrayPtr)
        intersections = intersections + ints

      end do
      !$omp end parallel do

      self % intersectionsTotal = self % intersectionsTotal + intersections
      
      call timerStop(self % timerTransport)

      ! Update RNG on master thread
      call self % rand % stride(self % pop + 1)

      ! Normalise flux estimate and combines with source
      call arrayPtr % normaliseFluxAndVolume(it)

      ! Calculate new k and accumulate stats
      self % keff = arrayPtr % calculateKeff(self % keff)
      if (isActive) then
        self % keffScore(1) = self % keffScore(1) + self % keff
        self % keffScore(2) = self % keffScore(2) + self % keff * self % keff
      end if

      ! Accumulate flux scores and tally results
      if (isActive) then
        call arrayPtr % accumulateFluxScores()
        if (associated(self % tally)) then
          call arrayPtr % tallyResults(self % tally)
          call self % tally % reportCycleEnd(dummyDungeon)
        end if
      end if

      ! Calculate proportion of cells that were hit
      hitRate = arrayPtr % getCellHitRate(it)
      call arrayPtr % wipeCellHits()

      ! Evaluate stopping criterion for active or inactive iterations
      if (isActive) then
        keepRunning = (itAct < self % active)
      else
        isActive = (itInac >= self % inactive)
      end if

      ! Set previous iteration flux to scalar flux
      ! and zero scalar flux
      call arrayPtr % resetFluxes()

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)
      transport_T = timerTime(self % timerTransport)
      self % time_transport = self % time_transport + transport_T

      ! Predict time to end
      end_T = real(self % active + self % inactive, defReal) * elapsed_T / it
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      if (self % loud) then
        call printFishLineR(it)
        print *
        print *, 'Iteration: ', numToChar(it), ' of ', numToChar(self % active + self % inactive)
        if(isActive) then
          print *,'Active iterations'
        else
          print *,'Inactive iterations'
        end if
        print *, 'Cell hit rate: ', trim(numToChar(real(hitRate,defReal)))
        print *, 'keff: ', trim(numToChar(real(self % keff,defReal)))
        print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
        print *, 'End time:     ', trim(secToChar(end_T))
        print *, 'Time to end:  ', trim(secToChar(T_toEnd))
        print *, 'Time per integration (ns): ', &
                trim(numToChar(transport_T*10**9/(self % nG * intersections)))
      end if

    end do

    ! Finalise flux and keff scores
    call arrayPtr % finaliseFluxScores(itAct)
    if (itAct /= 1) then
      Nm1 = 1.0_defReal/(itAct - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/itAct
    self % keffScore(1) = self % keffScore(1) * N1
    self % keffScore(2) = self % keffScore(2) * N1
    self % keffScore(2) = sqrt(Nm1*(self % keffScore(2) - &
            self % keffScore(1) * self % keffScore(1))) 

  end subroutine cycles

  !!
  !! Output calculation results to a file
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(randomRayPhysicsPackage), target, intent(inout) :: self
    type(outputFile), target                              :: out
    character(nameLen)                                    :: name
    class(visualiser), pointer                            :: vizPtr

    call out % init(self % outputFormat, filename = self % outputFile)
    
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
    
    name = 'Time_Per_Integration'
    call out % printValue(self % time_transport/(self % intersectionsTotal * self % nG),name)
    
    name = 'Clock_Time'
    call out % printValue(timerTime(self % timerMain),name)
    
    name = 'Hit_rate'
    call out % printValue(self % arrays % getAverageHitRate(),name)

    ! Print keff
    name = 'keff'
    call out % startBlock(name)
    call out % printResult(self % keffScore(1), self % keffScore(2), name)
    call out % endBlock()
    
    ! Print tally
    if (associated(self % tally)) then
      name = 'tally'
      call out % startBlock(name)
      call self % tally % print(out)
      call out % endBlock()
    end if
    
    ! Send all fluxes and SDs to VTK
    vizPtr => self % viz
    if (self % plotResults) call self % arrays % outputToVTK(vizPtr)

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
    if (self % lin) print *, "Using linear source"
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
    class(randomRayPhysicsPackage), intent(inout) :: self

    ! Clean Nuclear Data, Geometry and visualisation
    call ndreg_kill()
    call self % viz % kill()

    ! Clean contents
    self % geom    => null()
    self % geomIdx = 0
    self % timerMain = 0
    self % timerTransport = 0
    self % mgData    => null()
    self % nG        = 0
    self % nCells    = 0
    self % termination = ZERO
    self % dead        = ZERO
    self % pop         = 0
    self % inactive    = 0
    self % active      = 0
    self % cache       = .false.
    self % lin         = .false.
    self % plotResults = .false.
    self % keff        = ONE
    self % keffScore   = ZERO
    call self % arrays % kill()
    call self % XSData % kill()
    if (associated(self % tally)) then
      call self % tally % kill()
      self % tally => null()
    end if

  end subroutine kill

end module randomRayPhysicsPackage_class
