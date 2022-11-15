module ISMCPhysicsPackage_class

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
  use imcWeightClerk_class,           only : imcWeightResult

  ! Factories
  use transportOperatorFactory_func,  only : new_transportOperator
  use sourceFactory_func,             only : new_source

  implicit none

  private

  !!
  !! Physics Package for ISMC calculations
  !!
  type, public,extends(physicsPackage) :: ISMCPhysicsPackage
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
    integer(shortInt)  :: imcSourceN
    logical(defBool)   :: sourceGiven = .false.
    integer(shortInt)  :: nMat
    integer(shortint)  :: printUpdates

    ! Calculation components
    type(particleDungeon), pointer :: thisStep     => null()
    type(particleDungeon), pointer :: nextStep     => null()
    type(particleDungeon), pointer :: temp_dungeon => null()
    type(particleDungeon), allocatable :: matPhotons
      ! Note that other physics packages used pointers for these particleDungeons ( => null() )
      ! I found it easier to get 'allocatable' to work, unsure if this needs to be changed
    class(source), allocatable     :: inputSource
    class(source), allocatable     :: ISMCSource
    integer(shortInt), dimension(:,:), pointer :: emptyArray => null()

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

  end type ISMCPhysicsPackage

contains

  subroutine run(self)
    class(ISMCPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ ISMC CALCULATION /\/\"

    call self % steps(self % tally, self % imcWeightAtch, self % N_steps)
    call self % collectResults()

    print *
    print *, "\/\/ END OF ISMC CALCULATION \/\/"
    print *
  end subroutine

  !!
  !! Run steps for calculation
  !!
  subroutine steps(self, tally, tallyAtch, N_steps)
    class(ISMCPhysicsPackage), intent(inout)         :: self
    type(tallyAdmin), pointer,intent(inout)         :: tally
    type(tallyAdmin), pointer,intent(inout)         :: tallyAtch
    integer(shortInt), intent(in)                   :: N_steps
    integer(shortInt)                               :: i, j, matIdx
    integer(shortInt), dimension(:), allocatable    :: Nm, Np
    type(particle)                                  :: p
    real(defReal)                                   :: elapsed_T, end_T, T_toEnd
    real(defReal), dimension(:), allocatable        :: tallyEnergy
    class(IMCMaterial), pointer                     :: mat
    character(100),parameter :: Here ='steps (ISMCPhysicsPackage_class.f90)'
    class(tallyResult), allocatable                 :: tallyRes

    ! Attach nuclear data and RNG to particle
    p % pRNG   => self % pRNG
    p % geomIdx = self % geomIdx

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    allocate(tallyEnergy(self % nMat))

    ! Generate initial material photons
    call self % ISMCSource % generate(self % nextStep, self % pop, p % pRNG)

    open(unit = 10, file = 'temps.txt')
    open(unit = 11, file = 'pops.txt')

    allocate(Nm(self % nMat))
    allocate(Np(self % nMat))

    ! Build connections between materials
    call self % transOp % buildMajMap(p % pRNG, self % nucData)

    do i=1,N_steps

      write(10, '(8A)') numToChar(i)

      Nm = 0
      Np = 0

      ! Swap dungeons to store photons remaining from previous time step
      self % temp_dungeon => self % nextStep
      self % nextStep     => self % thisStep
      self % thisStep     => self % temp_dungeon
      call self % nextStep % cleanPop()

      ! Generate from input source
      if( self % sourceGiven ) then

        ! Limit number of particles in each zone
        call self % thisStep % reduceSize(self % limit, self % emptyArray)

        ! Generate new particles
        call self % inputSource % append(self % thisStep, self % pop, p % pRNG)

      end if

      !if(self % printSource == 1) then
      !  call self % thisStep % printToFile(trim(self % outputFile)//'_source'//numToChar(i))
      !end if

      call tally % reportCycleStart(self % thisStep)

      ! Update majorants for transport operator
      call self % transOp % updateMajorants(p % pRNG)

      ! Assign new maximum particle time
      p % timeMax = self % deltaT * i

      gen: do
        ! Obtain paticle from dungeon
        call self % thisStep % release(p)
        call self % geom % placeCoord(p % coords)

        ! Check particle type
        if (p % getType() /= P_PHOTON_MG .and. p % getType() /= P_MATERIAL_MG) then
          call fatalError(Here, 'Particle is not of type P_PHOTON_MG or P_MATERIAL_MG')
        end if

        ! For newly sourced particles, sample time uniformly within time step
        if (p % time /= ZERO .and. p % time /= self % deltaT*(i-1)) then
          call fatalError(Here, 'Particle released from dungeon has incorrect time')
        else if (p % time == ZERO) then
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
            if(p % fate == LEAK_FATE) exit history
            
            if(p % fate == AGED_FATE) then
                if(p % type == P_PHOTON) then
                  matIdx = p % matIdx()
                  Np(matIdx) = Np(matIdx) + 1
                else if( p % type == P_MATERIAL ) then
                  matIdx = p % matIdx()
                  Nm(matIdx) = Nm(matIdx) + 1
                else
                  call fatalError(Here, 'Incorrect type')
                end if
                ! Store particle for use in next time step
                p % fate = 0
                call self % nextStep % detain(p)
                exit history
            end if

            if (p % type == P_MATERIAL) then
              call fatalError(Here, 'Material particle should not undergo collision')
            end if

            call self % collOp % collide(p, tally, self % thisStep, self % nextStep)

            if(p % isDead) call fatalError(Here, 'Particle should not be dead, check that collision &
                                                  &operator is of type "ISMCMGstd"')

          end do history

        ! When dungeon is empty, exit
        if( self % thisStep % isEmpty() ) then
          exit gen
        end if

      end do gen

      ! Send end of cycle report
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
      call tallyAtch % getResult(tallyRes, 'imcWeight')

      select type(tallyRes)
        class is(imcWeightResult)
          do j = 1, self % nMat
            tallyEnergy(j) = tallyRes % imcWeight(j)
          end do
        class default
          call fatalError(Here, 'Invalid result has been returned')
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

      ! Reset tally for next cycle
      call tallyAtch % reset('imcWeight')

      print *, 'Completed: ', numToChar(i), ' of ', numToChar(N_steps)

      write(11, '(8A)') 'M ', numToChar(Nm)
      write(11, '(8A)') 'P ', numToChar(Np)

    end do

    close(10)
    close(11)

  end subroutine steps

  !!
  !! Print calculation results to file
  !!
  subroutine collectResults(self)
    class(ISMCPhysicsPackage), intent(inout)         :: self
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
    class(ISMCPhysicsPackage), intent(inout)        :: self
    class(dictionary), intent(inout)                :: dict
    class(dictionary),pointer                       :: tempDict
    type(dictionary)                                :: locDict1, locDict2, locDict3, locDict4
    integer(shortInt)                               :: seed_temp
    integer(longInt)                                :: seed
    character(10)                                   :: time
    character(8)                                    :: date
    character(:),allocatable                        :: string
    character(nameLen)                              :: nucData, geomName
    type(outputFile)                                :: test_out
    integer(shortInt)                               :: i
    character(nameLen), dimension(:), allocatable   :: mats
    class(IMCMaterial), pointer                     :: mat
    character(100), parameter :: Here ='init (ISMCPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % getOrDefault( self % limit, 'limit', self % pop)
    call dict % get( self % N_steps,'steps')
    call dict % get( self % deltaT,'timeStepSize')
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

    ! Read whether to print particle source per cycle
    call dict % getOrDefault(self % printSource, 'printSource', 0)

    ! Build Nuclear Data 
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'ISMCGeom'
    call gr_addGeom(geomName, tempDict)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Activate Nuclear Data *** All materials are active
    call ndReg_activate(self % particleType, nucData, self % geom % activeMats())
    self % nucData => ndReg_get(self % particleType)

    ! Read particle source definition
    if( dict % isPresent('source') ) then
      tempDict => dict % getDictPtr('source')
      call new_source(self % inputSource, tempDict, self % geom)
      self % sourceGiven = .true.
    end if

    ! Initialise ISMC source
    call locDict1 % init(2)
    call locDict1 % store('type', 'ismcSource')
    call locDict1 % store('N', self % pop)
    call new_source(self % ISMCSource, locDict1, self % geom)

    ! Build collision operator
    tempDict => dict % getDictPtr('collisionOperator')
    call self % collOp % init(tempDict)

    ! Build transport operator
    tempDict => dict % getDictPtr('transportOperator')
    call tempDict % store('deltaT', self % deltaT)
    call new_transportOperator(self % transOp, tempDict, self % geom)

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

    ! Set calculation type for material objects and provide time step
    do i=1, self % nMat 
      mat => IMCMaterial_CptrCast(self % nucData % getMaterial(i))
      call mat % setType(ISMC)
      call mat % setTimeStep(self % deltaT)
    end do

    ! Initialise imcWeight tally attachment
    call locDict2 % init(1)
    call locDict3 % init(2)
    call locDict4 % init(2)

    call locDict4 % store('type','materialMap')
    call locDict4 % store('materials', [mats])
    call locDict3 % store('type','imcWeightClerk')
    call locDict3 % store('map', locDict4)
    call locDict2 % store('imcWeight', locDict3)

    allocate(self % imcWeightAtch)
    call self % imcWeightAtch % init(locDict2)

    call self % tally % push(self % imcWeightAtch)

    ! Size particle dungeon
    allocate(self % thisStep)
    call self % thisStep % init(self % limit * self % nMat *2)
    allocate(self % nextStep)
    call self % nextStep % init(self % limit * self % nMat *2)

    allocate(self % emptyArray(3, self % limit * self % nMat * 2))

    call self % printSettings()

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(ISMCPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(ISMCPhysicsPackage), intent(in) :: self

    print *, repeat("<>",50)
    print *, "/\/\ ISMC CALCULATION /\/\"
    print *, "Source batches:       ", numToChar(self % N_steps)
    print *, "Population per batch: ", numToChar(self % pop)
    print *, "Initial RNG Seed:     ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings


end module ISMCPhysicsPackage_class
