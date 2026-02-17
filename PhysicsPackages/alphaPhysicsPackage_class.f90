module alphaPhysicsPackage_class

  use numPrecision
  use universalVariables
  use endfConstants
  use genericProcedures,              only : fatalError, numToChar, rotateVector
  use display_func,                   only : printFishLineR, statusMsg, printSectionStart, &
                                             printSectionEnd, printSeparatorLine
  use mpi_func,                       only : isMPIMaster, getWorkshare, getOffset, getMPIRank
#ifdef MPI
  use mpi_func,                       only : mpi_reduce, mpi_bcast, MPI_COMM_WORLD, &
                                             MASTER_RANK, MPI_SUM, MPI_DEFREAL, MPI_SHORTINT
#endif
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Particle classes and Random number generator
  use particle_class,                 only : particle, P_NEUTRON
  use particleDungeon_class,          only : particleDungeon
  use RNG_class,                      only : RNG

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx, &
                                             gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use geometryFactory_func,           only : new_geometry

  ! Fields
  use field_inter,                    only : field
  use uniFissSitesField_class,        only : uniFissSitesField, uniFissSitesField_TptrCast
  use fieldFactory_func,              only : new_field

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_activate    => activate ,&
                                             ndReg_display     => display, &
                                             ndReg_kill        => kill, &
                                             ndReg_get         => get ,&
                                             ndReg_getMatNames => getMatNames
  use nuclearDatabase_inter,          only : nuclearDatabase

  ! Sources
  use source_inter,                   only : source
  use sourceFactory_func,             only : new_source

  ! Operators
  use collisionOperator_class,        only : collisionOperator
  use transportOperator_inter,        only : transportOperator

  ! Tallies
  use tallyCodes
  use tallyAdmin_class,               only : tallyAdmin
  use tallyResult_class,              only : tallyResult
  use kAlphaAnalogClerk_class,        only : kAlphaResult

  ! Factories
  use transportOperatorFactory_func,  only : new_transportOperator

  ! Visualisation
  use visualiser_class,               only : visualiser

  implicit none
  private

  !!
  !! Physics Package for alpha eigenvalue calculations
  !!
  type, public,extends(physicsPackage) :: alphaPhysicsPackage
    private
    ! Building blocks
    class(nuclearDatabase), pointer        :: nucData       => null()
    class(geometry), pointer               :: geom          => null()
    integer(shortInt)                      :: geomIdx       = 0
    type(collisionOperator)                :: collOp
    class(transportOperator), allocatable  :: transOp
    class(source), allocatable             :: initSource
    class(RNG), pointer                    :: pRNG          => null()
    type(tallyAdmin),pointer               :: inactiveTally => null()
    type(tallyAdmin),pointer               :: activeTally   => null()
    type(tallyAdmin),pointer               :: inactiveAtch  => null()
    type(tallyAdmin),pointer               :: activeAtch    => null()
    class(uniFissSitesField),pointer       :: ufsField      => null()


    ! Settings
    integer(shortInt)  :: N_inactive
    integer(shortInt)  :: N_active
    integer(shortInt)  :: pop
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: particleType
    real(defReal)      :: keff_0
    real(defReal)      :: alpha_0
    real(defReal)      :: eta
    integer(shortInt)  :: bufferSize
    logical(defBool)   :: UFS = .false.

    ! Calculation components
    type(particleDungeon), pointer :: thisCycle    => null()
    type(particleDungeon), pointer :: nextCycle    => null()
    type(particleDungeon), pointer :: temp_dungeon => null()

    ! Timer bins
    integer(shortInt) :: timerMain
    real (defReal)    :: time_transport = 0.0
    real (defReal)    :: CPU_time_start
    real (defReal)    :: CPU_time_end

  contains
    procedure :: init
    procedure :: printSettings
    procedure :: cycles
    procedure :: generateInitialState
    procedure :: collectResults
    procedure :: run
    procedure :: buildActiveAtch
    procedure :: kill

  end type alphaPhysicsPackage

contains

  subroutine run(self)
    class(alphaPhysicsPackage), intent(inout) :: self

    call printSeparatorLine()
    call printSectionStart("ALPHA EIGENVALUE CALCULATION")

    ! Skip RNG state forward based on the process rank
    call self % pRNG % stride(getOffset(self % pop))

    call self % generateInitialState()
    call self % cycles(self % inactiveTally, self % inactiveAtch, self % N_inactive)
    call self % buildActiveAtch()
    call self % cycles(self % activeTally, self % activeAtch, self % N_active)
    
    ! Collect results from other processes
    call self % inactiveTally % collectDistributed()
    call self % activeTally % collectDistributed()

    if (isMpiMaster()) call self % collectResults()

    call statusMsg("")
    call printSectionEnd("END OF ALPHA EIGENVALUE CALCULATION")
    call statusMsg("")

  end subroutine run

  !!
  !!
  !!
  subroutine cycles(self, tally, tallyAtch, N_cycles)
    class(alphaPhysicsPackage), intent(inout) :: self
    type(tallyAdmin), pointer,intent(inout)   :: tally
    type(tallyAdmin), pointer,intent(inout)   :: tallyAtch
    integer(shortInt), intent(in)             :: N_cycles
    type(particleDungeon), save               :: buffer
    integer(shortInt)                         :: i, n, Nend, nParticles
    class(tallyResult),allocatable            :: res
    type(collisionOperator), save             :: collOp
    class(transportOperator),allocatable,save :: transOp
    type(RNG), target, save                   :: pRNG
    type(particle), save                      :: neutron
    real(defReal)                             :: k_new, alpha_new
    real(defReal)                             :: elapsed_T, end_T, T_toEnd
#ifdef MPI
    integer(shortInt)                         :: error, nTemp
#endif
    character(100),parameter :: Here ='cycles (alphaPhysicsPackage_class.f90)'
    !$omp threadprivate(neutron, buffer, collOp, transOp, pRNG)

    !$omp parallel
    ! Create particle buffer
    call buffer % init(self % bufferSize)

    ! Initialise neutron
    neutron % geomIdx = self % geomIdx

    ! Create a collision + transport operator which can be made thread private
    collOp = self % collOp
    transOp = self % transOp
    !$omp end parallel

    ! Set initial k-eff and alpha
    k_new = self % keff_0
    alpha_new = self % alpha_0

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    do i = 1, N_cycles

      ! Send start of cycle report
      call tally % reportCycleStart(self % thisCycle)

      nParticles = self % thisCycle % popSize()

      !$omp parallel do schedule(dynamic)
      gen: do n = 1, nParticles

        ! TODO: Further work to ensure reproducibility!
        ! Create RNG which can be thread private
        pRNG = self % pRNG
        neutron % pRNG => pRNG
        call neutron % pRNG % stride(n)

        ! Obtain particle current cycle dungeon
        call self % thisCycle % copy(neutron, n)

        bufferLoop: do
          call self % geom % placeCoord(neutron % coords)

          ! Set values for normalisation in the particle
          neutron % k_eff = k_new
          neutron % alpha = alpha_new
          neutron % eta = self % eta

          ! Save state
          call neutron % savePreHistory()

          ! Transport particle until its death
          history: do
            call transOp % transport(neutron, tally, buffer, self % nextCycle)
            if(neutron % isDead) exit history

            call collOp % collide(neutron, tally, buffer, self % nextCycle)
            if(neutron % isDead) exit history
          end do history

          ! Clear out buffer
          if (buffer % isEmpty()) then
            exit bufferLoop
          else
            call buffer % release(neutron)
          end if

        end do bufferLoop

      end do gen
      !$omp end parallel do

      call self % thisCycle % cleanPop()

      ! Update RNG
      call self % pRNG % stride(self % pop + 1)

      ! Send end of cycle report
      Nend = self % nextCycle % popSize()
      call tally % reportCycleEnd(self % nextCycle)

      if (self % UFS) then
        call self % ufsField % updateMap()
      end if

      ! Normalise population
      call self % nextCycle % combing(self % pop, self % pRNG)

      if (self % printSource /= 0) then
        call self % nextCycle % printToFile(trim(self % outputFile) // '_source' // numToChar(i) // &
                                            '_rank' // numToChar(getMPIRank()), self % printSource == BINARY_FILE)
      end if

      ! Flip cycle dungeons
      self % temp_dungeon => self % nextCycle
      self % nextCycle    => self % thisCycle
      self % thisCycle    => self % temp_dungeon

      ! Obtain estimate of k_eff
      call tallyAtch % getResult(res,'kAlpha')

      select type(res)
        class is(kAlphaResult)
          k_new = res % k
          alpha_new = res % alpha

        class default
          call fatalError(Here, 'Invalid result has been returned')

      end select

#ifdef MPI
      ! Broadcast k_eff and alpha obtained in the master to all processes
      call mpi_bcast([k_new, alpha_new], 2, MPI_DEFREAL, MASTER_RANK, MPI_COMM_WORLD)
#endif

      ! Load new k-eff estimate into next cycle dungeon
      self % nextCycle % k_eff = k_new

      ! Used to normalise fission source of the first active cycle
      ! Should hopefully be ~1
      self % keff_0 = k_new
      self % alpha_0 = alpha_new

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(N_cycles,defReal) * elapsed_T / i
      T_toEnd = max(ZERO, end_T - elapsed_T)

#ifdef MPI
      ! Print the population numbers referred to all processes to screen
      call mpi_reduce(nParticles, nTemp, 1, MPI_SHORTINT, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD, error)
      nParticles = nTemp
      call mpi_reduce(nEnd, nTemp, 1, MPI_SHORTINT, MPI_SUM, MASTER_RANK, MPI_COMM_WORLD, error)
      nEnd = nTemp
#endif

      ! Display progress
      call printFishLineR(i)
      call statusMsg("")
      call statusMsg("Cycle: " // numToChar(i) // " of " // numToChar(N_cycles))
      call statusMsg("Pop: " // numToChar(nParticles) // " -> " // numToChar(nEnd))
      call statusMsg("Elapsed time: " // trim(secToChar(elapsed_T)))
      call statusMsg("End time:     " // trim(secToChar(end_T)))
      call statusMsg("Time to end:  " // trim(secToChar(T_toEnd)))
      call tally % display()

    end do

    ! Load elapsed time
    self % time_transport = self % time_transport + elapsed_T

  end subroutine cycles

  !!
  !!
  !!
  subroutine generateInitialState(self)
    class(alphaPhysicsPackage), intent(inout) :: self
    character(100), parameter :: Here =' generateInitialState (alphaPhysicsPackage_class.f90)'

    ! Allocate and initialise particle Dungeons
    allocate(self % thisCycle)
    allocate(self % nextCycle)
    call self % thisCycle % init(2 * self % pop)
    call self % nextCycle % init(2 * self % pop)

    ! Generate initial source
    print *, "GENERATING INITIAL SOURCE"
    call self % initSource % generate(self % thisCycle, self % pop, self % pRNG)
    print *, "DONE!"
    
    ! Update RNG after source generation
    call self % pRNG % stride(self % pop)

  end subroutine generateInitialState

  !!
  !! Print calculation results to file
  !!
  subroutine collectResults(self)
    class(alphaPhysicsPackage), intent(inout) :: self
    type(outputFile)                          :: out
    character(nameLen)                        :: name

    call out % init(self % outputFormat, filename=self % outputFile)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Inactive_Cycles'
    call out % printValue(self % N_inactive,name)

    name = 'Active_Cycles'
    call out % printValue(self % N_active,name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Total_Transport_Time'
    call out % printValue(self % time_transport,name)

    ! Print Inactive tally
    name = 'inactive'
    call out % startBlock(name)
    call self % inactiveTally % print(out)
    call out % endBlock()
    
    ! Print Active attachment
    ! Is printed into the root block
    call self % activeAtch % print(out)

    name = 'active'
    call out % startBlock(name)
    call self % activeTally % print(out)
    call out % endBlock()

  end subroutine collectResults


  !!
  !! Initialise from individual components and dictionaries for inactive and active tally
  !!
  subroutine init(self, dict)
    class(alphaPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)          :: dict
    class(dictionary),pointer                 :: tempDict
    type(dictionary)                          :: locDict1, locDict2
    integer(shortInt)                         :: seed_temp
    integer(longInt)                          :: seed
    character(10)                             :: time
    character(8)                              :: date
    character(:),allocatable                  :: string
    character(nameLen)                        :: nucData, energy, geomName
    type(outputFile)                          :: test_out
    type(visualiser)                          :: viz
    class(field), pointer                     :: field
    character(100), parameter :: Here ='init (alphaPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_inactive,'inactive')
    call dict % get( self % N_active,'active')
    call dict % get( nucData, 'XSdata')
    call dict % get( energy, 'dataType')

    ! Parallel buffer size
    call dict % getOrDefault( self % bufferSize, 'buffer', 1000)

    ! Process type of data
    select case(energy)
      case('mg')
        self % particleType = P_NEUTRON_MG
      case('ce')
        self % particleType = P_NEUTRON_CE
      case default
        call fatalError(Here,"dataType must be 'mg' or 'ce'.")
    end select

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be caught early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Initialise RNG
    allocate(self % pRNG)

    ! *** It is a bit silly but dictionary cannot store longInt for now
    !     so seeds are limited to 32 bits (can be -ve)
    if( dict % isPresent('seed')) then
      call dict % get(seed_temp,'seed')

    else
      ! Obtain time string and hash it to obtain random seed
      call date_and_time(date, time)
      string = date // time
      call FNV_1(string,seed_temp)

    end if
    seed = seed_temp
    call self % pRNG % init(seed)

    ! Initial k_effective guess
    call dict % getOrDefault(self % keff_0,'keff_0', ONE)
    
    ! Initial alpha guess
    call dict % getOrDefault(self % alpha_0,'alpha_0', 0.1_defReal)
    
    ! Alpha stabilising factor
    call dict % getOrDefault(self % eta, 'eta', ONE)
    
    ! Read whether to print particle source per cycle
    call dict % getOrDefault(self % printSource, 'printSource', 0)
    if (self % printSource < NO_PRINTING .or. self % printSource > BINARY_FILE) then
      call fatalError(Here, 'printSource must be 0 (No printing), 1 (ASCII) or 2 (BINARY)')
    end if

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'alphaGeom'
    call new_geometry(tempDict, geomName)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Activate Nuclear Data *** All materials are active
    call ndReg_activate(self % particleType, nucData, self % geom % activeMats())
    self % nucData => ndReg_get(self % particleType)

    ! Call visualisation
    if (dict % isPresent('viz')) then
      print *, "Initialising visualiser"
      tempDict => dict % getDictPtr('viz')
      call viz % init(self % geom, tempDict)
      print *, "Constructing visualisation"
      call viz % makeViz()
      call viz % kill()
    endif

    ! Read uniform fission site option as a geometry field
    ! WARNING: Not sure how this will behave in alpha mode
    if (dict % isPresent('uniformFissionSites')) then
      self % ufs = .true.
      ! Build and initialise
      tempDict => dict % getDictPtr('uniformFissionSites')
      call new_field(tempDict, nameUFS)
      ! Save UFS field
      field => gr_fieldPtr(gr_fieldIdx(nameUFS))
      self % ufsField => uniFissSitesField_TptrCast(field)
      ! Initialise
      call self % ufsField % estimateVol(self % geom, self % pRNG, self % particleType)
    end if

    ! Read variance reduction option as a geometry field
    if (dict % isPresent('varianceReduction')) then
      ! Build and initialise
      tempDict => dict % getDictPtr('varianceReduction')
      call new_field(tempDict, nameWW)
    end if

    ! Build collision operator
    tempDict => dict % getDictPtr('collisionOperator')
    call self % collOp % init(tempDict)

    ! Build transport operator
    tempDict => dict % getDictPtr('transportOperator')
    call new_transportOperator(self % transOp, tempDict)

    ! Initialise active & inactive tally Admins
    tempDict => dict % getDictPtr('inactiveTally')
    allocate(self % inactiveTally)
    call self % inactiveTally % init(tempDict)

    tempDict => dict % getDictPtr('activeTally')
    allocate(self % activeTally)
    call self % activeTally % init(tempDict)

    ! Load Initial source
    ! Requires source specification for non-multiplying systems
    if (dict % isPresent('source')) then ! Load definition from file
      call new_source(self % initSource, dict % getDictPtr('source'), self % geom)

    else
      call locDict1 % init(3)
      call locDict1 % store('type', 'fissionSource')
      call locDict1 % store('data', trim(energy))
      call new_source(self % initSource, locDict1, self % geom)
      call locDict1 % kill()

    end if

    ! Initialise active and inactive tally attachments
    ! Inactive tally attachment
    ! Note: mpiSync ensures that alpha/k is synchronised between all processes each cycle
    call locDict1 % init(2)
    call locDict2 % init(4)

    call locDict2 % store('type','kAlphaAnalogClerk')
    call locDict2 % store('alpha_0',self % alpha_0)
    call locDict1 % store('kAlpha', locDict2)
    call locDict1 % store('display',['kAlpha'])
    call locDict1 % store('mpiSync', 1)

    allocate(self % inactiveAtch)
    call self % inactiveAtch % init(locDict1)

    call locDict2 % kill()
    call locDict1 % kill()

    ! Attach attachment to result tallies
    call self % inactiveTally % push(self % inactiveAtch)

    ! Don't create the active tally yet as the guess value
    ! strongly affects the convergence rate and may bias
    ! the active cycle estimate

    call self % printSettings()

  end subroutine init

  !!
  !! Build active attachment given better initial guess of alpha
  !! to prevent biasing the active cycle estimate
  !! Note: mpiSync ensures that alpha/k is synchronised between all processes each cycle
  !!
  subroutine buildActiveAtch(self)
    class(alphaPhysicsPackage), intent(inout) :: self
    type(dictionary)                          :: locDict1, locDict2

    call locDict1 % init(2)
    call locDict2 % init(4)

    call locDict2 % store('type','kAlphaAnalogClerk')
    call locDict2 % store('alpha_0',self % alpha_0)
    call locDict1 % store('kAlpha', locDict2)
    call locDict1 % store('display',['kAlpha'])
    call locDict1 % store('mpiSync', 1)

    allocate(self % activeAtch)
    call self % activeAtch % init(locDict1)

    call locDict2 % kill()
    call locDict1 % kill()

    call self % activeTally % push(self % activeAtch)

  end subroutine buildActiveAtch

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(alphaPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(alphaPhysicsPackage), intent(in) :: self

    call printSeparatorLine()
    call printSectionStart("ALPHA EIGENVALUE CALCULATION WITH K-ALPHA ITERATION METHOD")
    call statusMsg("Inactive Cycles:    " // numToChar(self % N_inactive))
    call statusMsg("Active Cycles:      " // numToChar(self % N_active))
    call statusMsg("Neutron Population: " // numToChar(self % pop))
    call statusMsg("Initial RNG Seed:   " // numToChar(self % pRNG % getSeed()))
    call statusMsg("")
    call printSeparatorLine()

  end subroutine printSettings

end module alphaPhysicsPackage_class
