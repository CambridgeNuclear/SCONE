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
  use materialMenu_mod,               only : mm_nMat           => nMat
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

    ! Settings
    integer(shortInt)  :: N_cycles
    real(defReal)      :: timeStepSize
    integer(shortInt)  :: pop
    character(pathLen) :: outputFile
    character(nameLen) :: outputFormat
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: particleType

    ! Calculation components
    type(particleDungeon), pointer :: thisCycle       => null()
    type(particleDungeon), pointer :: nextCycle       => null()
    type(particleDungeon), pointer :: tempDungeon     => null()
    class(source), allocatable     :: IMCSource
    !integer(shortInt)              :: nTimeStep

    ! Timer bins
    integer(shortInt)  :: timerMain
    real (defReal)     :: CPU_time_start
    real (defReal)     :: CPU_time_end

  contains
    procedure :: init
    procedure :: printSettings
    procedure :: cycles
    procedure :: collectResults
    procedure :: run
    procedure :: kill

  end type IMCPhysicsPackage

contains

  subroutine run(self)
    class(IMCPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ IMC CALCULATION /\/\"

    call self % cycles(self % tally, self % N_cycles)
    call self % collectResults()

    print *
    print *, "\/\/ END OF IMC CALCULATION \/\/"
    print *
  end subroutine

  !!
  !!
  !!
  subroutine cycles(self, tally, N_cycles)
    class(IMCPhysicsPackage), intent(inout)         :: self
    type(tallyAdmin), pointer,intent(inout)         :: tally
    integer(shortInt), intent(in)                   :: N_cycles
    integer(shortInt)                               :: i, N
    type(particle)                                  :: p
    real(defReal)                                   :: elapsed_T, end_T, T_toEnd
    class(IMCMaterial), pointer                     :: mat
    character(100),parameter :: Here ='cycles (IMCPhysicsPackage_class.f90)'

    N = self % pop

    ! Attach nuclear data and RNG to particle
    p % pRNG   => self % pRNG
    p % k_eff = ONE
    p % geomIdx = self % geomIdx

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    do i=1,N_cycles

      ! Store photons remaining from previous cycle
      self % tempDungeon = self % nextCycle
      call self % nextCycle % cleanPop()

      ! Send start of cycle report
      call self % IMCSource % generate(self % thisCycle, N, p % pRNG)
      if(self % printSource == 1) then
        call self % thisCycle % printToFile(trim(self % outputFile)//'_source'//numToChar(i))
      end if

      call tally % reportCycleStart(self % thisCycle)

      self % thisCycle % endOfStepTime   = i * self % timeStepSize
      self % tempDungeon % endOfStepTime = i * self % timeStepSize 

      gen: do
        ! Obtain paticle from dungeon
        call self % thisCycle % release(p)
        call self % geom % placeCoord(p % coords)

        ! Save state
        call p % savePreHistory()

          ! Transport particle until its death
          history: do
            call self % transOp % transport(p, tally, self % thisCycle, self % nextCycle)
            if(p % isDead) exit history
            
            if(p % fate == TIME_FATE) then
                ! Store particle for use in next time step
                call self % nextCycle % detain(p)
                exit history
            end if

            call self % collOp % collide(p, tally, self % thisCycle, self % nextCycle)
            if(p % isDead) exit history

          end do history

        ! When both dungeons empty, exit
        if( self % thisCycle % isEmpty() .and. self % tempDungeon % isEmpty()) exit gen

        ! When source dungeon is emptied, switch to using particles remaining from previous cycle
        if( self % thisCycle % isEmpty()) then
          print *, "THIS CYCLE EMPTIED OF SOURCE >>> SWITCHING TO CENSUS PARTICLES"
          self % thisCycle = self % tempDungeon
          call self % tempDungeon % cleanPop()
        end if

      end do gen

      ! Send end of cycle report
      call tally % reportCycleEnd(self % thisCycle)

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(N_cycles,defReal) * elapsed_T / i
      T_toEnd = max(ZERO, end_T - elapsed_T)

      mat => IMCMaterial_CptrCast(self % nucData % getMaterial(1))

      !call mat % updateTemp()
 
      ! Display progress
      call printFishLineR(i)
      print *
      call mat % updateTemp()
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
    class(IMCPhysicsPackage), intent(inout)         :: self
    type(outputFile)                                :: out
    character(nameLen)                              :: name

    call out % init(self % outputFormat)

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

    call out % writeToFile(self % outputFile)

  end subroutine collectResults


  !!
  !! Initialise from individual components and dictionaries for source and tally
  !!
  subroutine init(self, dict)
    class(IMCPhysicsPackage), intent(inout)         :: self
    class(dictionary), intent(inout)                :: dict
    class(dictionary),pointer                       :: tempDict
    integer(shortInt)                               :: seed_temp
    integer(longInt)                                :: seed
    character(10)                                   :: time
    character(8)                                    :: date
    character(:),allocatable                        :: string
    character(nameLen)                              :: nucData, energy, geomName
    type(outputFile)                                :: test_out
    integer(shortInt)                               :: i
    character(100), parameter :: Here ='init (IMCPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_cycles,'cycles')
    call dict % get( self % timeStepSize,'timeStepSize')
    call dict % get( nucData, 'XSdata')
    call dict % get( energy, 'dataType')

    ! Process type of data
    select case(energy)
      case('mg')
        self % particleType = P_PHOTON_MG
      !case('ce')
      !  self % particleType = P_PHOTON_CE
      case default
        call fatalError(Here,"dataType must be 'mg' or 'ce'.")
    end select

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
    geomName = 'IMCGeom'
    call gr_addGeom(geomName, tempDict)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Activate Nuclear Data *** All materials are active
    call ndReg_activate(self % particleType, nucData, self % geom % activeMats())
    self % nucData => ndReg_get(self % particleType)

    ! Read particle source definition
    tempDict => dict % getDictPtr('source')
    call new_source(self % IMCSource, tempDict, self % geom)

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
    allocate(self % thisCycle)
    call self % thisCycle % init(3 * self % pop)
    allocate(self % nextCycle)
    call self % nextCycle % init(3 * self % pop)
    allocate(self % tempDungeon)
    call self % tempDungeon % init(3 * self % pop)

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
    print *, "Source batches:       ", numToChar(self % N_cycles)
    print *, "Population per batch: ", numToChar(self % pop)
    print *, "Initial RNG Seed:     ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings

  !!
  !! Return time at end of current time step
  !!
  !function endOfStepTime(self) result(time)
  !  implicit none
  !  class(IMCPhysicsPackage), intent(in) :: self
  !  real(defReal)                        :: time
  !
  !  time = self % timeStepSize * self % nTimeStep
  !end function endOfStepTime

end module IMCPhysicsPackage_class
