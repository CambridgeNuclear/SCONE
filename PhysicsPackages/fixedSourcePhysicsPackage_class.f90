module fixedSourcePhysicsPackage_class

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

  ! Fields
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

  ! Visualisation
  use visualiser_class,               only : visualiser

  implicit none
  private

  !!
  !! Physics Package for fixed source calculations
  !!
  type, public,extends(physicsPackage) :: fixedSourcePhysicsPackage
    private
    ! Building blocks
    class(nuclearDatabase), pointer        :: nucData => null()
    class(geometry), pointer               :: geom    => null()
    integer(shortInt)                      :: geomIdx = 0
    type(collisionOperator)                :: collOp
    class(transportOperator), allocatable  :: transOp
    class(RNG), pointer                    :: pRNG    => null()
    type(tallyAdmin),pointer               :: tally   => null()

    ! Settings
    integer(shortInt)  :: N_cycles
    integer(shortInt)  :: pop
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: particleType
    integer(shortInt)  :: bufferSize
    integer(shortInt)  :: bufferShift

    ! Calculation components
    type(particleDungeon), pointer :: thisCycle       => null()
    type(particleDungeon), pointer :: commonBuffer    => null()
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

  end type fixedSourcePhysicsPackage

contains

  subroutine run(self)
    class(fixedSourcePhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ FIXED SOURCE CALCULATION /\/\"

    call self % cycles(self % tally, self % N_cycles)
    call self % collectResults()

    print *
    print *, "\/\/ END OF FIXED SOURCE CALCULATION \/\/"
    print *
  end subroutine

  !!
  !!
  !!
  subroutine cycles(self, tally, N_cycles)
    class(fixedSourcePhysicsPackage), intent(inout) :: self
    type(tallyAdmin), pointer,intent(inout)         :: tally
    integer(shortInt), intent(in)                   :: N_cycles
    integer(shortInt)                               :: i, n, nParticles
    integer(shortInt), save                         :: j, bufferExtra
    type(particle), save                            :: p, transferP
    type(particleDungeon), save                     :: buffer
    type(collisionOperator), save                   :: collOp
    class(transportOperator), allocatable, save     :: transOp
    type(RNG), target, save                         :: pRNG
    real(defReal)                                   :: elapsed_T, end_T, T_toEnd
    character(100),parameter :: Here ='cycles (fixedSourcePhysicsPackage_class.f90)'
    !$omp threadprivate(p, buffer, collOp, transOp, pRNG, j, bufferExtra, transferP)

    !$omp parallel
    ! Create particle buffer
    call buffer % init(self % bufferSize)

    ! Initialise neutron
    p % geomIdx = self % geomIdx
    p % k_eff = ONE

    ! Create a collision + transport operator which can be made thread private
    collOp = self % collOp
    transOp = self % transOp
    !$omp end parallel

    nParticles = self % pop

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    do i=1,N_cycles

      ! Send start of cycle report
      call self % fixedSource % generate(self % thisCycle, nParticles, self % pRNG)
      if(self % printSource == 1) then
        call self % thisCycle % printToFile(trim(self % outputFile)//'_source'//numToChar(i))
      end if

      call tally % reportCycleStart(self % thisCycle)

      !$omp parallel do schedule(dynamic)
      gen: do n = 1, nParticles

        ! TODO: Further work to ensure reproducibility!
        ! Create RNG which can be thread private
        pRNG = self % pRNG
        p % pRNG => pRNG
        call p % pRNG % stride(n)

        ! Obtain particle from dungeon
        call self % thisCycle % copy(p, n)
        p % showerID = n

        bufferLoop: do

          call self % geom % placeCoord(p % coords)

          ! Save state
          call p % savePreHistory()

          ! Transport particle until its death
          history: do
            call transOp % transport(p, tally, buffer, buffer)
            if(p % isDead) exit history

            call collOp % collide(p, tally, buffer, buffer)
            if(p % isDead) exit history
          end do history

          ! If buffer is quite full, shift some particles to the commonBuffer
          if (associated(self % commonBuffer) .and. (buffer % popSize() > self % bufferShift)) then
            bufferExtra = buffer % popSize() - self % bufferShift
            do j = 1, bufferExtra
              call buffer % release(transferP)
              call self % commonBuffer % detainCritical(transferP)
            end do
          end if

          ! Clear out buffer
          if (.not. buffer % isEmpty()) then
            call buffer % release(p)
            cycle
          end if

          ! Clear out common queue
          ! Note the apparently redundant critical sections (one here in PP, one in the dungeon).
          ! This is to prevent the situation where two threads both enter the conditional and compete
          ! for the final particle in the dungeon. The first thread would pop the particle while the
          ! second would try to pop from an empty dungeon.
          if (associated(self % commonBuffer)) then
            !$omp critical
            if (.not. self % commonBuffer % isEmpty()) then
              call self % commonBuffer % releaseCritical(p)
            end if
            !$omp end critical
            if (p % isDead) exit bufferLoop
          else
            exit bufferLoop
          end if

        end do bufferLoop

      end do gen
      !$omp end parallel do

      ! Update RNG
      call self % pRNG % stride(self % pop)

      ! Send end of cycle report
      call tally % reportCycleEnd(self % thisCycle)

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(N_cycles,defReal) * elapsed_T / i
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      call printFishLineR(i)
      print *
      print *, 'Source batch: ', numToChar(i), ' of ', numToChar(N_cycles)
      print *, 'Pop:          ', numToChar(self % pop)
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      call tally % display()
    end do

  end subroutine cycles

  !!
  !! Print calculation results to file
  !!
  subroutine collectResults(self)
    class(fixedSourcePhysicsPackage), intent(inout) :: self
    type(outputFile)                                :: out
    character(nameLen)                              :: name

    call out % init(self % outputFormat, filename=self % outputFile)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Source_batches'
    call out % printValue(self % N_cycles,name)

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
    class(fixedSourcePhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)                :: dict
    class(dictionary),pointer                       :: tempDict
    integer(shortInt)                               :: seed_temp, commonBufferSize
    integer(longInt)                                :: seed
    character(10)                                   :: time
    character(8)                                    :: date
    character(:),allocatable                        :: string
    character(nameLen)                              :: nucData, energy, geomName
    type(outputFile)                                :: test_out
    type(visualiser)                                :: viz
    character(100), parameter :: Here ='init (fixedSourcePhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_cycles,'cycles')
    call dict % get( nucData, 'XSdata')
    call dict % get( energy, 'dataType')

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

    ! Parallel buffer size
    call dict % getOrDefault( self % bufferSize, 'buffer', 50)

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
    geomName = 'fixedSourceGeom'
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

    ! Read variance reduction option as a geometry field
    if (dict % isPresent('varianceReduction')) then
      ! Build and initialise
      tempDict => dict % getDictPtr('varianceReduction')
      call new_field(tempDict, nameWW)
    end if

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
    ! Note no need to oversize for fixed source calculation
    allocate(self % thisCycle)
    call self % thisCycle % init(self % pop)

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

    call self % printSettings()

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(fixedSourcePhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(fixedSourcePhysicsPackage), intent(in) :: self

    print *, repeat("<>",50)
    print *, "/\/\ FIXED SOURCE CALCULATION /\/\"
    print *, "Source batches:       ", numToChar(self % N_cycles)
    print *, "Population per batch: ", numToChar(self % pop)
    print *, "Initial RNG Seed:     ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings

end module fixedSourcePhysicsPackage_class
