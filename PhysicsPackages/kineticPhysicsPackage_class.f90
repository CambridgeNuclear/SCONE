module kineticPhysicsPackage_class

  use numPrecision
  use universalVariables
  use endfConstants
  use genericProcedures,              only : fatalError, printFishLineR, numToChar, rotateVector
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Particle classes and Random number generator
  use particle_class,                 only : particle, P_NEUTRON
  use particleDungeon_class,          only : particleDungeon
  use source_inter,                   only : source
  use RNG_class,                      only : RNG

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx
  use geometryFactory_func,           only : new_geometry

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_activate    => activate ,&
                                             ndReg_display     => display, &
                                             ndReg_kill        => kill, &
                                             ndReg_get         => get ,&
                                             ndReg_getMatNames => getMatNames
  use nuclearDatabase_inter,          only : nuclearDatabase
  use neutronMaterial_inter,          only : neutronMaterial, neutronMaterial_CptrCast
  use ceNeutronMaterial_class,        only : ceNeutronMaterial
  use mgNeutronMaterial_inter,        only : mgNeutronMaterial
  use fissionCE_class,                only : fissionCE, fissionCE_TptrCast
  use fissionMG_class,                only : fissionMG, fissionMG_TptrCast
  use ceNeutronDatabase_inter,        only : ceNeutronDatabase, ceNeutronDatabase_CptrCast

  ! Operators
  use collisionOperator_class,        only : collisionOperator
  use transportOperator_inter,        only : transportOperator

  ! Tallies
  use tallyCodes
  use tallyAdmin_class,               only : tallyAdmin

  ! Factories
  use transportOperatorFactory_func,  only : new_transportOperator
  use sourceFactory_func,             only : new_source

  implicit none
  private

  !!
  !! Physics Package for time dependent calculations
  !!
  type, public,extends(physicsPackage) :: kineticPhysicsPackage
    private
    ! Building blocks
    class(nuclearDatabase), pointer        :: nucData   => null()
    class(geometry), pointer               :: geom      => null()
    integer(shortInt)                      :: geomIdx   = 0
    type(collisionOperator)                :: collOp
    class(transportOperator), allocatable  :: transOp
    class(RNG), pointer                    :: pRNG      => null()
    type(tallyAdmin),pointer               :: tally     => null()
    type(tallyAdmin),pointer               :: tallyAtch => null()

    ! Settings
    integer(shortInt)  :: pop
    integer(shortInt)  :: N_cycles
    real(defReal)      :: k0
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: particleType
    integer(shortInt)  :: bufferSize
    integer(shortInt)  :: bufferShift = 0
    logical(defBool)   :: forcedPrecursorDecay
    real(defReal), dimension(:), allocatable :: timeBounds

    ! Calculation components
    type(particleDungeon), pointer :: thisTime       => null()
    type(particleDungeon), pointer :: nextTime       => null()
    type(particleDungeon), pointer :: tempTime       => null()
    type(particleDungeon), pointer :: thisPrecursors => null() 
    type(particleDungeon), pointer :: nextPrecursors => null() 
    type(particleDungeon), pointer :: tempPrecursors => null() 
    type(particleDungeon), pointer :: commonBuffer   => null()
    class(source), allocatable     :: fixedSource

    ! Timer bins
    integer(shortInt) :: timerMain
    real (defReal)    :: CPU_time_start
    real (defReal)    :: CPU_time_end

  contains
    procedure :: init
    procedure :: printSettings
    procedure :: cycles
    procedure :: collectResults
    procedure :: run
    procedure :: kill

  end type kineticPhysicsPackage

contains

  subroutine run(self)
    class(kineticPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ KINETIC CALCULATION /\/\"

    call self % cycles()
    call self % collectResults()

    print *
    print *, "\/\/ END OF KINETIC CALCULATION \/\/"
    print *
  end subroutine

  !!
  !!
  !!
  subroutine cycles(self)
    class(kineticPhysicsPackage), intent(inout)  :: self
    integer(shortInt)                            :: b, t, n
    type(particle), save                         :: p, pNew
    type(particleDungeon), save                  :: buffer
    type(collisionOperator), save                :: collOp
    class(transportOperator), allocatable, save  :: transOp
    type(RNG), target, save                      :: pRNG
    real(defReal) , save                         :: decayT, dT
    real(defReal)                                :: elapsed_T, end_T, T_toEnd, tPrev, tNext
    integer(shortInt), save                      :: j, bufferExtra
    character(100),parameter :: Here ='cycles (kineticPhysicsPackage_class.f90)'
    !$omp threadprivate(p, pNew, buffer, collOp, transOp, pRNG, decayT, dT, j, bufferExtra)

    !$omp parallel
    ! Create particle buffer
    call buffer % init(self % bufferSize)

    ! Initialise neutron
    p % geomIdx = self % geomIdx
    p % k_eff = self % k0
    pNew = p

    ! Create thread private collision + transport operator
    collOp = self % collOp
    transOp = self % transOp
    !$omp end parallel

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)
    
    batch: do b = 1, self % N_cycles

      ! Produce initial source    
      call self % fixedSource % generate(self % thisTime, self % pop, self % pRNG)
      
      call self % tally % reportCycleStart(self % thisTime)

      time: do t = 1, size(self % timeBounds)
    
        ! Set start-of-step and end-of-step time
        if (t == 1) then
          tPrev = ZERO
        else
          tPrev = self % timeBounds(t-1)
        end if
        tNext = self % timeBounds(t)

        ! Set all particles in the dungeon to having the start-of-step time
        ! TODO: allow sources which are distributed in time
        call self % thisTime % setTime(tPrev)

        ! Add delayed neutrons to the sampled particles
        if (self % thisPrecursors % popSize() > 0) then

          ! Perform population control on precursors
          if (self % thisPrecursors % popSize() > self % pop) then
            call self % thisPrecursors % precursorCombing(self % pop, self % pRNG, tPrev, tNext)
          end if

          ! Sample which delayed neutrons to simulate
          ! Implicit precursor sampling
          if (self % forcedPrecursorDecay) then

            ! Force precursors to decay during the step
            !$omp parallel do
            do n = 1, self % thisPrecursors % popSize()
        
              pRNG = self % pRNG
              call pRNG % stride(n)
              
              call self % thisPrecursors % copy(p, n)
              dT = tNext - tPrev
              decayT = pRNG % get() * dT + tPrev

              ! Weight adjustment
              call p % forcedPrecursorDecay(decayT, dT, pNew)
              
              ! Add particle to the cycle
              call self % thisTime % detain(pNew)

            end do
            !$omp end parallel do

          ! Analog precursor sampling
          else

            call self % nextPrecursors % cleanPop()

            !$omp parallel do
            do n = 1, self % thisPrecursors % popSize()
              call self % thisPrecursors % copy(p, n)
              
              pRNG = self % pRNG
              call pRNG % stride(n)
          
              ! Sample time to decay
              p % time = tPrev - log(pRNG % get()) / p % lambda

              if (p % time < tNext) then
                call p % emitDelayedNeutron()
                call self % thisTime % detain(p)
              else
                call self % nextPrecursors % detain(p)
              end if

            end do
            !$omp end parallel do
            
            ! Flip precursor dungeons
            self % tempPrecursors => self % nextPrecursors
            self % nextPrecursors => self % thisPrecursors
            self % thisPrecursors => self % tempPrecursors

          end if
            
          ! Update RNG
          call self % pRNG % stride(self % thisPrecursors % popSize() + 1)

        end if

        ! Comb the population to keep the number of simulated particles constant
        call self % thisTime % combing(self % pop, self % pRNG)

        !$omp parallel do schedule(dynamic)
        gen: do n = 1, self % pop
          
          call self % thisTime % copy(p, n)
          
          pRNG = self % pRNG
          p % pRNG => pRNG
          call p % pRNG % stride(n)

          bufferLoop: do

            p % fate = no_FATE
            p % timeMax = tNext
          
            call self % geom % placeCoord(p % coords)
            
            call p % savePreHistory()

            ! Transport particle until its death
            history: do

              ! Catch for undecayed precursors
              if (p % isDead) exit history
              
              call transOp % transport(p, self % tally, buffer, buffer)
              if(p % isDead) exit history
         
              ! Particle hit the time boundary     
              if(p % fate == AGED_FATE) then
                call self % nextTime % detain(p)
                exit history
              end if

              call collOp % collide(p, self % tally, buffer, buffer)
              if(p % isDead) exit history
            end do history

            ! If buffer is quite full, shift some particles to the commonBuffer
            if (associated(self % commonBuffer) .and. (buffer % popSize() > self % bufferShift)) then
              bufferExtra = buffer % popSize() - self % bufferShift
              do j = 1, bufferExtra
                call buffer % release(p)
                call self % commonBuffer % detainCritical(p)
              end do
            end if
          
            ! Clear out buffer
            if ((.not. buffer % isEmpty()) .or. associated(self % commonBuffer)) then

              if (.not. buffer % isEmpty()) then
                call buffer % release(p)
              
              ! Clear out common queue
              ! Note the apparently redundant critical sections (one here in PP, one in the dungeon).
              ! This is to prevent the situation where two threads both enter the conditional and compete
              ! for the final particle in the dungeon. The first thread would pop the particle while the
              ! second would try to pop from an empty dungeon.
              elseif (associated(self % commonBuffer)) then
                p = p
                p % isDead = .true.
                !$omp critical
                if (.not. self % commonBuffer % isEmpty()) then
                  call self % commonBuffer % releaseCritical(p)
                end if
                !$omp end critical
                if (p % isDead) exit bufferLoop
              else
                exit bufferLoop
              end if


              ! Set RNG for the broodID?

              ! Is the particle a precursor?
              if (p % isPrecursor()) then

                ! If forced decay, obtain the particle weight and place the
                ! remainder in the dungeon
                if (self % forcedPrecursorDecay) then
            
                  dT = tNext - p % time
                  decayT = pRNG % get() * dT + p % time

                  ! Weight adjustment
                  call p % forcedPrecursorDecay(decayT, dT, pNew)
                  call self % thisPrecursors % detain(p)

                  ! Swap particles and simulate
                  p = pNew

                ! If analog sampling, does the precursor decay in this step?
                ! If not, store in the precursor dungeon and 'kill' the particle
                else

                  ! Precursor won't decay until later.
                  ! Sample decay time
                  decayT = p % time - log(pRNG % get()) / p % lambda

                  ! Set the particle as dead, skipping the history loop.
                  if (decayT >= tNext) then
                    call self % thisPrecursors % detain(p)
                    p % isDead = .true.

                  ! Precursor will decay in this step
                  else
                    call p % emitDelayedNeutron()
                  end if

                end if 
              end if
            else
              exit bufferLoop

            end if

          end do bufferLoop
        end do gen
        !$omp end parallel do

        ! Update RNG
        call self % pRNG % stride(self % pop + 1)

        call self % thisTime % cleanPop()
        self % tempTime => self % nextTime
        self % nextTime => self % thisTime
        self % thisTime => self % tempTime

      end do time
        
      call self % tally % reportCycleEnd(self % thisTime)
      
      ! Clean for the next batch
      call self % thisTime % cleanPop()
      call self % thisPrecursors % cleanPop()
      
      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(self % N_cycles,defReal) * elapsed_T / b
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      call printFishLineR(t)
      print *
      print *, 'Batch: ', numToChar(b), ' of ', numToChar(self % N_cycles)
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      call self % tally % display()

    end do batch

  end subroutine cycles

  !!
  !! Print calculation results to file
  !!
  subroutine collectResults(self)
    class(kineticPhysicsPackage), intent(inout) :: self
    type(outputFile)                            :: out
    character(nameLen)                          :: name

    call out % init(self % outputFormat, filename=self % outputFile)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Source_batches'
    call out % printValue(self % N_cycles,name)

    name = 'Time_increment'
    call out % printValue(self % timeBounds(1),name)

    name = 'Time_bins'
    call out % printValue(size(self % timeBounds),name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Transport_time'
    call out % printValue(timerTime(self % timerMain),name)

    ! Print tally
    call self % tally % print(out)

  end subroutine collectResults

  !!
  !! Initialise from individual components and dictionaries for source and tally
  !!
  subroutine init(self, dict)
    class(kineticPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)            :: dict
    class(dictionary),pointer                   :: tempDict
    type(dictionary)                            :: locDict1, locDict2
    integer(shortInt)                           :: seed_temp, i, nSteps, commonBufferSize
    integer(longInt)                            :: seed
    real(defReal)                               :: dt, t0
    character(10)                               :: time
    character(8)                                :: date
    character(:),allocatable                    :: string
    character(nameLen)                          :: nucData, energy, geomName
    type(outputFile)                            :: test_out
    character(100), parameter :: Here ='init (kineticPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_cycles,'cycles')
    call dict % get( nSteps,'timeSteps')
    call dict % get( nucData, 'XSdata')
    call dict % get( energy, 'dataType')
    call dict % get( dt, 'dt')
    call dict % getOrDefault(self % k0, 'k0', ONE)

    if (self % pop < 1) then
      call fatalError(Here, 'Population must be greater than zero: '//numToChar(self % pop))
    end if
    if (self % N_cycles < 1) then
      call fatalError(Here, 'Number of cycles must be greater than zero: '//numToChar(self % N_cycles))
    end if
    if (self % k0 < ZERO) then
      call fatalError(Here, 'Unphysical k0 provided: '//numToChar(self % k0))
    end if
    if (dt < ZERO) then
      call fatalError(Here, 'Timestep width must be positive: '//numToChar(dt))
    end if
    if (nSteps < 1) then
      call fatalError(Here, 'Must have at least one timestep: '//numToChar(nSteps))
    end if

    allocate(self % timeBounds(nSteps))
    t0 = ZERO
    do i = 1, nSteps
      self % timeBounds(i) = t0 + dt
      t0 = self % timeBounds(i)
    end do

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
    ! Initialise output file before calculation (so mistake in format will be cought early)
    call dict % getOrDefault(self % outputFormat, 'outputFormat', 'asciiMATLAB')
    call test_out % init(self % outputFormat)

    ! Parallel buffer size
    call dict % getOrDefault( self % bufferSize, 'buffer', 50)

    ! Whether to use analog or implicit precursors treatment (default = implicit)
    call dict % getOrDefault(self % forcedPrecursorDecay, 'forcedPrecursorDecay', .true.)
    
    ! Is the common buffer turned on? Set the size if so
    if (dict % isPresent('commonBufferSize')) then
      call dict % get(commonBufferSize,'commonBufferSize')
      allocate(self % commonBuffer)
      call self % commonBuffer % init(commonBufferSize)

      ! Set threshold at which to shift particles from private buffer
      ! to common buffer
      call dict % getOrDefault(self % bufferShift, 'bufferShift', 10)
      if (self % bufferShift > self % bufferSize) call fatalError(Here, &
              'Buffer size should be greater than the shift threshold')
    end if

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

    ! Read whether to print particle source per cycle
    call dict % getOrDefault(self % printSource, 'printSource', 0)

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'kineticGeom'
    call new_geometry(tempDict, geomName)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Activate Nuclear Data *** All materials are active
    call ndReg_activate(self % particleType, nucData, self % geom % activeMats())
    self % nucData => ndReg_get(self % particleType)

    ! Read particle source definition
    tempDict => dict % getDictPtr('source')
    call new_source(self % fixedSource, tempDict, self % geom)

    ! Build collision operator
    tempDict => dict % getDictPtr('collisionOperator')
    call self % collOp % init(tempDict)

    ! Build transport operator
    tempDict => dict % getDictPtr('transportOperator')
    call new_transportOperator(self % transOp, tempDict)

    ! Initialise tally Admin
    tempDict => dict % getDictPtr('tally')
    allocate(self % tally)
    call self % tally % init(tempDict)

    ! Size particle dungeon
    allocate(self % thisTime)
    allocate(self % nextTime)
    call self % thisTime % init(3 * self % pop)
    call self % nextTime % init(3 * self % pop)

    ! Size precursor dungeon
    allocate(self % thisPrecursors)
    allocate(self % nextPrecursors)
    call self % thisPrecursors % init(3 * self % pop)
    call self % nextPrecursors % init(3 * self % pop)

    call self % printSettings()
    
    ! Initialise tally for useful diagnostics
    call locDict1 % init(2)
    call locDict2 % init(2)

    call locDict2 % store('type','keffAnalogClerk')
    call locDict1 % store('keff', locDict2)
    call locDict1 % store('display',['keff'])

    allocate(self % tallyAtch)
    call self % tallyAtch % init(locDict1)
    call self % tally % push(self % tallyAtch)

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(kineticPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(kineticPhysicsPackage), intent(in) :: self
    real(defReal)                                  :: TStart, Tstop, Tincrement

    TStart = ZERO
    Tstop = self % timeBounds(size(self % timeBounds))
    Tincrement = self % timeBounds(1)
    print *, repeat("<>",50)
    print *, "/\/\ KINETIC CALCULATION /\/\"
    print *, "Time grid [start, stop, increment]: ", numToChar(TStart), numToChar(Tstop), numToChar(Tincrement)
    print *, "Initial Population:                 ", numToChar(self % pop)
    print *, "Initial RNG Seed:                   ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings

end module kineticPhysicsPackage_class
