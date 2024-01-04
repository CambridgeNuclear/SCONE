module implicitPhysicsPackage_class

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
  use particle_class,                 only : particle, P_PHOTON, P_MATERIAL
  use particleDungeon_class,          only : particleDungeon
  use source_inter,                   only : source
  use RNG_class,                      only : RNG

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_addGeom => addGeom, &
                                             gr_geomIdx  => geomIdx
  use discretiseGeom_class,           only : discretise

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat ,&
                                             mm_matName        => matName
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_activate    => activate ,&
                                             ndReg_display     => display, &
                                             ndReg_kill        => kill, &
                                             ndReg_get         => get ,&
                                             ndReg_getMatNames => getMatNames
  use nuclearDatabase_inter,          only : nuclearDatabase
  use mgIMCDatabase_inter,            only : mgIMCDatabase, mgIMCDatabase_CptrCast
  use IMCMaterial_inter,              only : IMCMaterial, IMCMaterial_CptrCast

  ! Operators
  use collisionOperator_class,        only : collisionOperator
  use transportOperator_inter,        only : transportOperator

  ! Tallies
  use tallyCodes
  use tallyAdmin_class,               only : tallyAdmin
  use tallyResult_class,              only : tallyResult
  use energyWeightClerk_class,        only : energyWeightClerkResult

  ! Factories
  use transportOperatorFactory_func,  only : new_transportOperator
  use sourceFactory_func,             only : new_source

  use simulationTime_class
  use energyGridRegistry_mod,         only : define_energyGrid

  implicit none

  private

  ! Calculation Type
  integer(shortInt), parameter, public :: IMC  = 1
  integer(shortInt), parameter, public :: ISMC = 2

  !!
  !! Physics Package for Implicit Monte Carlo calculations
  !!
  !! Settings:
  !!  method                  -> IMC or ISMC
  !!  pop                     -> For IMC, approx. number of new particles to generate per time step
  !!                          -> For ISMC, starting population of material particles
  !!  limit                   -> Size of particle dungeons
  !!  steps                   -> Number of time steps to simulate
  !!  timeStep                -> Time step to be used
  !!  printUpdates (OPTIONAL) -> Prints material update info for first N materials
  !!
  type, public,extends(physicsPackage) :: implicitPhysicsPackage
    private
    ! Building blocks
    class(mgIMCDatabase), pointer          :: nucData  => null()
    class(geometry), pointer               :: geom     => null()
    integer(shortInt)                      :: geomIdx = 0
    type(collisionOperator)                :: collOp
    class(transportOperator), allocatable  :: transOp
    class(RNG), pointer                    :: pRNG    => null()
    type(tallyAdmin),pointer               :: tally   => null()
    type(tallyAdmin),pointer               :: energyWeightAtch => null()

    ! Settings
    integer(shortInt)  :: N_steps
    integer(shortInt)  :: pop
    integer(shortInt)  :: limit
    integer(shortInt)  :: method
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: nMat
    integer(shortInt)  :: printUpdates

    ! Calculation components
    type(particleDungeon), pointer :: thisStep     => null()
    type(particleDungeon), pointer :: nextStep     => null()
    type(particleDungeon), pointer :: temp_dungeon => null()
    class(source), allocatable     :: inputSource
    class(source), allocatable     :: matSource

    ! Timer bins
    integer(shortInt)  :: timerMain
    real (defReal)     :: CPU_time_start
    real (defReal)     :: CPU_time_end

  contains
    procedure :: init
    procedure :: printSettings
    procedure :: steps
    procedure :: collectResults
    procedure :: run
    procedure :: kill

  end type implicitPhysicsPackage

contains

  subroutine run(self)
    class(implicitPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ IMPLICIT CALCULATION /\/\"

    call self % steps(self % tally, self % energyWeightAtch, self % N_steps)
    call self % collectResults()

    print *
    print *, "\/\/ END OF IMPLICIT CALCULATION \/\/"
    print *
  end subroutine

  !!
  !! Run steps for calculation
  !!
  !!
  !! Notes differences between IMC and ISMC regarding particle generation:
  !!
  !!  -> IMC generates particles from material emission as well as input source, ISMC only
  !!     generates from input source.
  !!
  !!  -> Particles for IMC are killed when absorbed, but for ISMC remain as material particles.
  !!     This allows IMC to keep dungeon pop under control far better than ISMC. For ISMC we
  !!     generate particles such that pop is at limit at end of calculation, with runtime per
  !!     time step increasing approximatelty linearly as the calculation progresses
  !!
  subroutine steps(self, tally, tallyAtch, N_steps)
    class(implicitPhysicsPackage), intent(inout)    :: self
    type(tallyAdmin), pointer,intent(inout)         :: tally
    type(tallyAdmin), pointer,intent(inout)         :: tallyAtch
    integer(shortInt), intent(in)                   :: N_steps
    integer(shortInt)                               :: i, j, N, nFromMat, num, nParticles
    type(particle), save                            :: p
    real(defReal)                                   :: sourceWeight, elapsed_T, end_T, T_toEnd
    class(IMCMaterial), pointer                     :: mat
    character(100),parameter :: Here ='steps (implicitPhysicsPackage_class.f90)'
    class(tallyResult), allocatable                 :: tallyRes
    type(collisionOperator), save                   :: collOp
    class(transportOperator), allocatable, save     :: transOp
    type(RNG), target, save                         :: pRNG
    !$omp threadprivate(p, collOp, transOp, pRNG)

    !$omp parallel
    p % geomIdx = self % geomIdx

    ! Create a collision + transport operator which can be made thread private
    collOp = self % collOp
    transOp = self % transOp

    !$omp end parallel

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Generate starting population of material particles for ISMC
    if (self % method == ISMC) then
      call self % matSource % append(self % thisStep, self % pop, self % pRNG)
      nFromMat = 0
    end if

    do i=1,N_steps

      ! Generate particles while staying below dungeon limit (see note in subroutine description)
      if (self % method == IMC) then

        ! Reduce number of particles to generate if close to limit
        N = self % pop
        if (N + self % thisStep % popSize() > self % limit) then
          ! Fleck and Cummings IMC Paper, eqn 4.11
          N = self % limit - self % thisStep % popSize() - self % nMat - 1
          N = max(1, N)
        end if

        ! Calculate proportion to be generated from input source
        if (allocated(self % inputSource)) then
          sourceWeight = self % inputSource % sourceWeight
          nFromMat = int(N * (1 - sourceWeight/(sourceWeight + self % nucData % getEmittedRad())))
          ! Generate from input source
          call self % inputSource % append(self % thisStep, N - nFromMat, self % pRNG)
        end if

        ! Add to dungeon particles emitted from material
        call self % matSource % append(self % thisStep, nFromMat, self % pRNG)

      ! ISMC particle generation
      else if (allocated(self % inputSource)) then

        ! Generate particles such that pop is almost at limit at calculation end
        N = (self % limit - self % thisStep % popSize()) / (N_steps-i+1)
        call self % inputSource % append(self % thisStep, N, self % pRNG)

      end if

      if(self % printSource == 1) then
        call self % thisStep % printToFile(trim(self % outputFile)//'_source'//numToChar(i))
      end if

      call tally % reportCycleStart(self % thisStep)

      nParticles = self % thisStep % popSize()

      !$omp parallel do schedule(dynamic)
      gen: do num = 1, nParticles

        ! Create RNG which can be thread private
        pRNG = self % pRNG
        p % pRNG => pRNG
        call p % pRNG % stride(num)

        ! Obtain paticle from dungeon
        call self % thisStep % release(p)
        call self % geom % placeCoord(p % coords)

        ! Check particle type
        if (p % getType() /= P_PHOTON_MG .and. p % getType() /= P_MATERIAL_MG) then
          call fatalError(Here, 'Particle is not of type P_PHOTON_MG or P_MATERIAL_MG')
        end if

        ! Assign maximum particle time
        p % timeMax = time % stepEnd

        ! Check for time errors
        if (p % time < time % stepStart .or. p % time >= time % stepEnd) then
          call fatalError(Here, 'Particle time not within time step')
        end if

        ! Save state
        call p % savePreHistory()

          ! Transport particle until its death
          history: do
            call transOp % transport(p, tally, self % thisStep, self % nextStep)
            if(p % isDead) exit history

            if(p % fate == AGED_FATE) then
                ! Store particle for use in next time step
                p % fate = 0
                call self % nextStep % detain(p)
                exit history
            end if

            call collOp % collide(p, tally, self % thisStep, self % nextStep)

            ! Cycle if particle history not yet completed
            if (self % method == ISMC .or. p % type == P_PHOTON) cycle history

            ! If P_MATERIAL and IMC, kill particle and exit
            p % isDead = .true.
            p % fate = ABS_FATE
            call tally % reportHist(p)
            exit history

          end do history

      end do gen
      !$omp end parallel do

      ! Update RNG
      call self % pRNG % stride(nParticles)

      ! Send end of time step report
      call tally % reportCycleEnd(self % thisStep)

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(N_steps,defReal) * elapsed_T / i
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      call printFishLineR(i)
      print *
      print *
      print *, 'Source batch: ', numToChar(i), ' of ', numToChar(N_steps)
      print *, 'Pop:          ', numToChar(self % nextStep % popSize())
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      call tally % display()

      ! Obtain energy deposition tally results
      call tallyAtch % getResult(tallyRes, 'energyWeightTally')

      ! Update material properties using tallied energy
      select type(tallyRes)
        class is(energyWeightClerkResult)
          call self % nucData % updateProperties(tallyRes % materialEnergy, self % printUpdates)
        class default
          call fatalError(Here, 'Tally result class should be energyWeightClerkResult')
      end select

      ! Reset tally for next time step
      if (i /= N_Steps) call tallyAtch % reset('energyWeightTally')

      ! Advance to next time step
      call nextStep()

      ! Swap dungeons to store photons remaining from previous time step
      self % temp_dungeon => self % nextStep
      self % nextStep     => self % thisStep
      self % thisStep     => self % temp_dungeon
      call self % nextStep % cleanPop()

    end do

    ! Output final mat temperatures
    open(unit = 10, file = 'temps.txt')
    do j = 1, self % nMat
      mat => IMCMaterial_CptrCast(self % nucData % getMaterial(j))
      write(10, '(8A)') mm_matName(j), numToChar(mat % getTemp())
    end do
    close(10)

    ! Output final radiation energies
    open(unit = 11, file = 'radEnergy.txt')
    select type(tallyRes)
      class is(energyWeightClerkResult)
        write(11, '(8A)') numToChar(tallyRes % radiationEnergy)
      class default
          call fatalError(Here, 'Tally result class should be energyWeightClerkResult')
    end select
    close(11)

  end subroutine steps

  !!
  !! Print calculation results to file
  !!
  subroutine collectResults(self)
    class(implicitPhysicsPackage), intent(inout) :: self
    type(outputFile)                             :: out
    character(nameLen)                           :: name

    call out % init(self % outputFormat)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Source_batches'
    call out % printValue(self % N_steps,name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Transport_time'
    call out % printValue(timerTime(self % timerMain),name)

    ! Print tally
    call self % tally % print(out)

    call out % writeToFile(self % outputFile)

  end subroutine collectResults


  !!
  !! Initialise from individual components and dictionaries for source and tally
  !!
  subroutine init(self, dict)
    class(implicitPhysicsPackage), intent(inout)  :: self
    class(dictionary), intent(inout)              :: dict
    class(dictionary), pointer                    :: tempDict
    type(dictionary)                              :: locDict1, locDict2, locDict3
    integer(shortInt)                             :: seed_temp
    integer(longInt)                              :: seed
    character(10)                                 :: time
    character(8)                                  :: date
    character(:),allocatable                      :: string
    character(nameLen)                            :: nucData, geomName
    type(outputFile)                              :: test_out
    integer(shortInt)                             :: i
    character(nameLen), dimension(:), allocatable :: mats
    real(defReal)                                 :: timeStep
    type(dictionary),target                       :: newGeom, newData
    character(nameLen)                            :: method, units
    character(100), parameter :: Here ='init (implicitPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Get method
    call dict % getOrDefault(method, 'method', 'IMC')
    select case(method)
      case ('IMC')
        self % method = IMC
      case('ISMC')
        self % method = ISMC
      case default
        call fatalError(Here, 'Unrecognised method')
    end select

    ! Read calculation settings
    call dict % get(self % pop,'pop')
    call dict % get(self % limit,'limit')
    call dict % get(self % N_steps,'steps')
    call dict % get(timeStep,'timeStep')
    call dict % getOrDefault(self % printUpdates, 'printUpdates', 0)
    nucData = 'mg'

    ! Set time step after changing units if necessary
    call dict % getOrDefault(units, 'units', 'ns')
    select case(units)
      case('s')
        ! No change needed
      case('ns')
        ! Convert time step from ns to s
        timeStep = timeStep * 1e-9_defReal
      case('marshak')
        ! Special case where a = c = 1
        timeStep = timeStep/lightSpeed
        print *, 'WARNING: For Marshak wave, still need to manually change radiationConstant to 1'
      case default
        call fatalError(Here, 'Unrecognised units')
    end select
    call setStep(timeStep)

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

    ! Read whether to print particle source each time step
    call dict % getOrDefault(self % printSource, 'printSource', 0)

    ! Initialise energy grid in multi-frequency case
    if (dict % isPresent('energyGrid')) then
      tempDict => dict % getDictPtr('energyGrid')
      call define_energyGrid(nucData, tempDict)
      print *, 'Energy grid defined: ', nucData
    end if

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr('nuclearData'))

    ! Build geometry
    geomName = 'IMCGeom'
    call gr_addGeom(geomName, dict % getDictPtr('geometry'))
    self % geomIdx =  gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Activate Nuclear Data *** All materials are active
    call ndReg_activate(P_PHOTON_MG, nucData, self % geom % activeMats())
    self % nucData => mgIMCDatabase_CptrCast(ndReg_get(P_PHOTON_MG))

    call newGeom % kill()
    call newData % kill()

    ! Initialise material source
    call locDict1 % init(2)
    call locDict1 % store('type', 'materialSource')
    ! Tell source if we are using IMC or ISMC
    call locDict1 % store('calcType', self % method)
    call new_source(self % matSource, locDict1, self % geom)
    call locDict1 % kill()

    ! Read external particle source definition
    if( dict % isPresent('source') ) then
      tempDict => dict % getDictPtr('source')
      call new_source(self % inputSource, tempDict, self % geom)
    end if

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

    ! Provide materials with calculation type
    call self % nucData % setCalcType(self % method)

    ! Store number of materials
    self % nMat = mm_nMat()
    self % printUpdates = min(self % printUpdates, self % nMat)

    ! Create array of material names
    allocate(mats(self % nMat))
    do i=1, self % nMat
      mats(i) = mm_matName(i)
    end do


    ! Initialise energy weight tally attachment
    call locDict1 % init(1)
    call locDict2 % init(2)
    call locDict3 % init(2)
    call locDict3 % store('type','materialMap')
    call locDict3 % store('materials', [mats])
    call locDict2 % store('type','energyWeightClerk')
    call locDict2 % store('map', locDict3)
    call locDict1 % store('energyWeightTally', locDict2)

    allocate(self % energyWeightAtch)
    call self % energyWeightAtch % init(locDict1)
    call self % tally % push(self % energyWeightAtch)

    ! Size particle dungeons
    allocate(self % thisStep)
    call self % thisStep % init(self % limit)
    allocate(self % nextStep)
    call self % nextStep % init(self % limit)

    call self % printSettings()

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(implicitPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(implicitPhysicsPackage), intent(in) :: self

    print *, repeat("<>",50)
    print *, "/\/\ IMC CALCULATION /\/\"
    print *, "Source batches:       ", numToChar(self % N_steps)
    print *, "Population per batch: ", numToChar(self % pop)
    print *, "Initial RNG Seed:     ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings

end module implicitPhysicsPackage_class
