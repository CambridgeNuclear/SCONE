module IMCPhysicsPackage_class

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
  use particle_class,                 only : particle, P_PHOTON
  use particleDungeon_class,          only : particleDungeon
  use source_inter,                   only : source
  use RNG_class,                      only : RNG

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_addGeom => addGeom, &
                                             gr_geomIdx  => geomIdx

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
  use IMCMaterial_inter,              only : IMCMaterial, IMCMaterial_CptrCast
  use mgIMCMaterial_inter,            only : mgIMCMaterial

  ! Operators
  use collisionOperator_class,        only : collisionOperator
  use transportOperator_inter,        only : transportOperator

  ! Tallies
  use tallyCodes
  use tallyAdmin_class,               only : tallyAdmin
  use tallyResult_class,              only : tallyResult
  use absorptionClerk_class,          only : absClerkResult

  ! Factories
  use transportOperatorFactory_func,  only : new_transportOperator
  use sourceFactory_func,             only : new_source

  implicit none

  private

  !!
  !! Physics Package for IMC calculations
  !!
  type, public,extends(physicsPackage) :: IMCPhysicsPackage
    private
    ! Building blocks
    class(nuclearDatabase), pointer        :: nucData => null()
    class(geometry), pointer               :: geom    => null()
    integer(shortInt)                      :: geomIdx = 0
    type(collisionOperator)                :: collOp
    class(transportOperator), allocatable  :: transOp
    class(RNG), pointer                    :: pRNG    => null()
    type(tallyAdmin),pointer               :: tally   => null()
    type(tallyAdmin),pointer               :: imcWeightAtch => null()

    ! Settings
    integer(shortInt)  :: N_steps
    integer(shortInt)  :: pop
    integer(shortInt)  :: limit
    real(defReal)      :: deltaT
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: particleType
    logical(defBool)   :: sourceGiven = .false.
    integer(shortInt)  :: nMat
    integer(shortInt)  :: printUpdates

    ! Calculation components
    type(particleDungeon), pointer :: thisStep     => null()
    type(particleDungeon), pointer :: nextStep     => null()
    type(particleDungeon), pointer :: temp_dungeon => null()
    class(source), allocatable     :: inputSource
    class(source), allocatable     :: IMCSource

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

  end type IMCPhysicsPackage

contains

  subroutine run(self)
    class(IMCPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ IMC CALCULATION /\/\"

    call self % steps(self % tally, self % imcWeightAtch, self % N_steps)
    call self % collectResults()

    print *
    print *, "\/\/ END OF IMC CALCULATION \/\/"
    print *
  end subroutine

  !!
  !! Run steps for calculation
  !!
  subroutine steps(self, tally, tallyAtch, N_steps)
    class(IMCPhysicsPackage), intent(inout)         :: self
    type(tallyAdmin), pointer,intent(inout)         :: tally
    type(tallyAdmin), pointer,intent(inout)         :: tallyAtch
    integer(shortInt), intent(in)                   :: N_steps
    integer(shortInt)                               :: i, j, N
    type(particle)                                  :: p
    real(defReal)                                   :: elapsed_T, end_T, T_toEnd
    real(defReal), dimension(:), allocatable        :: tallyEnergy
    class(IMCMaterial), pointer                     :: mat
    character(100),parameter :: Here ='steps (IMCPhysicsPackage_class.f90)'
    class(tallyResult), allocatable                 :: tallyRes

    ! Attach nuclear data and RNG to particle
    p % pRNG   => self % pRNG
    p % geomIdx = self % geomIdx

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    allocate(tallyEnergy(self % nMat))

    do i=1,N_steps

      ! Swap dungeons to store photons remaining from previous time step
      self % temp_dungeon => self % nextStep
      self % nextStep     => self % thisStep
      self % thisStep     => self % temp_dungeon
      call self % nextStep % cleanPop()

      ! Select number of particles to generate - for now this is an equal number from each zone
      N = self % pop
      if(N + self % thisStep % popSize() > self % limit) then
        ! Fleck and Cummings IMC Paper, eqn 4.11
        N = self % limit - self % thisStep % popSize() - self % nMat - 1
      end if
      N = int(N/self % nMat)
      if (N == 0) N = 1

      ! Add to particle dungeon
      do j=1, self % nMat
        mat => IMCMaterial_CptrCast(self % nucData % getMaterial(j))
        if (mat % getTemp() > 0) then
          call self % IMCSource % append(self % thisStep, N, p % pRNG, j)
        end if
      end do

      ! Generate from input source
      if( self % sourceGiven ) then
        call self % inputSource % append(self % thisStep, 0, p % pRNG)
      end if

      if(self % printSource == 1) then
        call self % thisStep % printToFile(trim(self % outputFile)//'_source'//numToChar(i))
      end if

      call tally % reportCycleStart(self % thisStep)

      ! Assign new maximum particle time
      p % timeMax = self % deltaT * i

      gen: do
        ! Obtain paticle from dungeon
        call self % thisStep % release(p)
        call self % geom % placeCoord(p % coords)

        ! Check particle type
        if (p % getType() /= P_PHOTON_MG) then
          call fatalError(Here, 'Particle is not of type P_PHOTON_MG')
        end if

        ! For newly sourced particles, sample time uniformly within time step
        if (p % time == ZERO) then
          p % time = (p % pRNG % get() + i-1) * self % deltaT
        end if

        ! Check for time errors
        if (p % time >= p % timeMax .or. p % time < self % deltaT*(i-1)) then
          call fatalError(Here, 'Particle time is not within timestep bounds')
        else if (p % time /= p % time) then
          call fatalError(Here, 'Particle time is NaN')
        end if

        ! Save state
        call p % savePreHistory()

          ! Transport particle until its death
          history: do
            call self % transOp % transport(p, tally, self % thisStep, self % nextStep)
            if(p % isDead) exit history

            if(p % fate == AGED_FATE) then
                ! Store particle for use in next time step
                p % fate = 0
                call self % nextStep % detain(p)
                exit history
            end if

            call self % collOp % collide(p, tally, self % thisStep, self % nextStep)

            if(p % isDead) exit history

          end do history

        ! When dungeon is empty, exit
        if (self % thisStep % isEmpty()) exit gen

      end do gen

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
      call tallyAtch % getResult(tallyRes, 'imcWeightTally')

      select type(tallyRes)
        class is(absClerkResult)
          do j = 1, self % nMat
            tallyEnergy(j) = tallyRes % clerkResults(j)
          end do
        class default
          call fatalError(Here, 'Tally result class should be absClerkResult')
      end select

      ! Update material properties
      do j = 1, self % nMat
        mat => IMCMaterial_CptrCast(self % nucData % getMaterial(j))
        if (j <= self % printUpdates) then
          print *
          print *, "Material update:  ", mm_matName(j)
          call mat % updateMat(tallyEnergy(j), .true.)
        else
          call mat % updateMat(tallyEnergy(j), .false.)
        end if
      end do
      print *

      ! Reset tally for next time step
      call tallyAtch % reset('imcWeightTally')

      print *, 'Completed: ', numToChar(i), ' of ', numToChar(N_steps)

    end do

  end subroutine steps

  !!
  !! Print calculation results to file
  !!
  subroutine collectResults(self)
    class(IMCPhysicsPackage), intent(inout)         :: self
    type(outputFile)                                :: out
    character(nameLen)                              :: name

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
    class(IMCPhysicsPackage), intent(inout)         :: self
    class(dictionary), intent(inout)                :: dict
    class(dictionary),pointer                       :: tempDict
    type(dictionary)                                :: locDict1, locDict2, locDict3, locDict4, locDict5
    integer(shortInt)                               :: seed_temp
    integer(longInt)                                :: seed
    character(10)                                   :: time
    character(8)                                    :: date
    character(:),allocatable                        :: string
    character(nameLen)                              :: nucData, geomName
    type(outputFile)                                :: test_out
    integer(shortInt)                               :: i
    class(IMCMaterial), pointer                     :: mat
    character(nameLen), dimension(:), allocatable   :: mats
    character(100), parameter :: Here ='init (IMCPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Read calculation settings
    call dict % get(self % pop,'pop')
    call dict % get(self % limit, 'limit')
    call dict % get(self % N_steps,'steps')
    call dict % get(self % deltaT,'timeStepSize')
    call dict % getOrDefault(self % printUpdates, 'printUpdates', 0)
    self % particleType = P_PHOTON_MG
    nucData = 'mg'

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be cought early)
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

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'IMCGeom'
    call gr_addGeom(geomName, tempDict)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Activate Nuclear Data *** All materials are active
    call ndReg_activate(self % particleType, nucData, self % geom % activeMats())
    self % nucData => ndReg_get(self % particleType)

    ! Read particle source definition
    if( dict % isPresent('source') ) then
      tempDict => dict % getDictPtr('source')
      call tempDict % store('deltaT', self % deltaT)
      call new_source(self % inputSource, tempDict, self % geom)
      self % sourceGiven = .true.
    end if

    ! Initialise IMC source
    call locDict1 % init(1)
    call locDict1 % store('type', 'imcSource')
    call new_source(self % IMCSource, locDict1, self % geom)

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

    ! Store number of materials
    self % nMat = mm_nMat()

    ! Create array of material names
    allocate(mats(self % nMat))
    do i=1, self % nMat
      mats(i) = mm_matName(i)
    end do

    ! Provide each material with time step
    do i=1, self % nMat
      mat => IMCMaterial_CptrCast(self % nucData % getMaterial(i))
      call mat % setTimeStep(self % deltaT)
    end do

    ! Initialise imcWeight tally attachment
    call locDict2 % init(1)
    call locDict3 % init(4)
    call locDict4 % init(2)
    call locDict5 % init(1)

    call locDict5 % store('type', 'weightResponse')
    call locDict4 % store('type','materialMap')
    call locDict4 % store('materials', [mats])
    call locDict3 % store('response', ['imcWeightResponse'])
    call locDict3 % store('imcWeightResponse', locDict5)
    call locDict3 % store('type','absorptionClerk')
    call locDict3 % store('map', locDict4)
    call locDict2 % store('imcWeightTally', locDict3)

    allocate(self % imcWeightAtch)
    call self % imcWeightAtch % init(locDict2)

    call self % tally % push(self % imcWeightAtch)

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
    class(IMCPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(IMCPhysicsPackage), intent(in) :: self

    print *, repeat("<>",50)
    print *, "/\/\ IMC CALCULATION /\/\"
    print *, "Source batches:       ", numToChar(self % N_steps)
    print *, "Population per batch: ", numToChar(self % pop)
    print *, "Initial RNG Seed:     ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings


end module IMCPhysicsPackage_class
