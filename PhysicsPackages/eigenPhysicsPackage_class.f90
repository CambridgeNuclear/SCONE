module eigenPhysicsPackage_class

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
  use RNG_class,                      only : RNG

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use cellGeometry_inter,             only : cellGeometry

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
  use keffAnalogClerk_class,          only : keffResult

  ! Factories
  use geometryFactory_func,           only : new_cellGeometry_ptr
  use transportOperatorFactory_func,  only : new_transportOperator

  ! Visualisation
  use visualiser_class,               only : visualiser

  implicit none
  private

  !!
  !! Physics Package for eigenvalue calculations
  !!
  type, public,extends(physicsPackage) :: eigenPhysicsPackage
    private
    ! Building blocks
    class(nuclearDatabase), pointer        :: nucData       => null()
    class(cellGeometry), pointer           :: geom          => null()
    type(collisionOperator)                :: collOp
    class(transportOperator), allocatable  :: transOp
    class(source), allocatable             :: initSource
    class(RNG), pointer                    :: pRNG          => null()
    type(tallyAdmin),pointer               :: inactiveTally => null()
    type(tallyAdmin),pointer               :: activeTally   => null()
    type(tallyAdmin),pointer               :: inactiveAtch  => null()
    type(tallyAdmin),pointer               :: activeAtch    => null()


    ! Settings
    integer(shortInt)  :: N_inactive
    integer(shortInt)  :: N_active
    integer(shortInt)  :: pop
    character(pathLen) :: outputFile
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: particleType

    ! Calculation components
    type(particleDungeon), pointer :: thisCycle    => null()
    type(particleDungeon), pointer :: nextCycle    => null()
    type(particleDungeon), pointer :: temp_dungeon => null()

    ! Timer bins
    integer(shortInt) :: timerMain
    real (defReal)    :: CPU_time_start
    real (defReal)    :: CPU_time_end

  contains
    procedure :: init
    procedure :: printSettings
    procedure :: cycles
    procedure :: generateInitialState
    procedure :: collectResults
    procedure :: run
    procedure :: kill

  end type eigenPhysicsPackage

contains

  subroutine run(self)
    class(eigenPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ EIGENVALUE CALCULATION /\/\"

    call self % generateInitialState()
    call self % cycles(self % inactiveTally, self % inactiveAtch, self % N_inactive)
    call self % cycles(self % activeTally, self % activeAtch, self % N_active)
    call self % collectResults()

    print *
    print *, "\/\/ END OF EIGENVALUE CALCULATION \/\/"
    print *
  end subroutine

  !!
  !!
  !!
  subroutine cycles(self, tally, tallyAtch, N_cycles)
    class(eigenPhysicsPackage), intent(inout) :: self
    type(tallyAdmin), pointer,intent(inout)   :: tally
    type(tallyAdmin), pointer,intent(inout)   :: tallyAtch
    integer(shortInt), intent(in)             :: N_cycles
    integer(shortInt)                         :: i, Nstart, Nend
    class(tallyResult),allocatable            :: res
    type(particle)                            :: neutron
    real(defReal)                             :: k_old, k_new
    real(defReal)                             :: elapsed_T, end_T, T_toEnd
    character(100),parameter :: Here ='cycles (eigenPhysicsPackage_class.f90)'


    ! Attach nuclear data and RNG to neutron
    neutron % pRNG   => self % pRNG

    ! Set initiial k-eff
    k_new = ONE

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    do i=1,N_cycles

      ! Send start of cycle report
      Nstart = self % thisCycle % popSize()
      call tally % reportCycleStart(self % thisCycle)

      gen: do
        ! Obtain paticle from current cycle dungeon
        call self % thisCycle % release(neutron)
        call self % geom % placeCoord(neutron % coords)

        ! Set k-eff for normalisation in the particle
        neutron % k_eff = k_new

        ! Save state
        call neutron % savePreHistory()

          ! Transport particle untill its death
          history: do
            call self % transOp % transport(neutron, tally, self % thisCycle, self % nextCycle)
            if(neutron % isDead) exit history

            call self % collOp % collide(neutron, tally ,self % thisCycle, self % nextCycle)
            if(neutron % isDead) exit history
          end do history

        if( self % thisCycle % isEmpty()) exit gen
      end do gen

      ! Send end of cycle report
      Nend = self % nextCycle % popSize()
      call tally % reportCycleEnd(self % nextCycle)

      ! Normalise population
      call self % nextCycle % normSize(self % pop, neutron % pRNG)

      if(self % printSource == 1) then
        call self % nextCycle % printToFile(trim(self % outputFile)//'_source'//numToChar(i))
      end if

      ! Flip cycle dungeons
      self % temp_dungeon => self % nextCycle
      self % nextCycle    => self % thisCycle
      self % thisCycle    => self % temp_dungeon

      ! Obtain estimate of k_eff
      call tallyAtch % getResult(res,'keff')

      select type(res)
        class is(keffResult)
          k_new = res % keff(1)

        class default
          call fatalError(Here, 'Invalid result has been returned')

      end select

      ! Load new k-eff estimate into next cycle dungeon
      k_old = self % nextCycle % k_eff
      self % nextCycle % k_eff = k_new

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(N_cycles,defReal) * elapsed_T / i
      T_toEnd = max(ZERO, end_T - elapsed_T)


      ! Display progress
      call printFishLineR(i)
      print *
      print *, 'Cycle: ', numToChar(i), ' of ', numToChar(N_cycles)
      print *, 'Pop: ', numToChar(Nstart) , ' -> ', numToChar(Nend)
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      call tally % display()
    end do
  end subroutine cycles

  !!
  !!
  !!
  subroutine generateInitialState(self)
    class(eigenPhysicsPackage), intent(inout) :: self
    character(100), parameter :: Here =' generateInitialState( eigenPhysicsPackage_class.f90)'

    ! Allocate and initialise particle Dungeons
    allocate(self % thisCycle)
    allocate(self % nextCycle)

    call self % thisCycle % init(3 * self % pop)
    call self % nextCycle % init(3 * self % pop)

    ! Generate initial surce
    print *, "GENERATING INITIAL FISSION SOURCE"
    call self % initSource % generate(self % thisCycle, self % pop, self % pRNG)
    print *, "DONE!"

  end subroutine generateInitialState

  !!
  !! Print calculation results to file
  !!
  subroutine collectResults(self)
    class(eigenPhysicsPackage), intent(inout) :: self
    type(outputFile)                          :: out
    character(pathLen)                        :: path
    character(nameLen)                        :: name

    name = 'asciiMATLAB'
    call out % init(name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Inactive_Cycles'
    call out % printValue(self % N_inactive,name)

    name = 'Active_Cycles'
    call out % printValue(self % N_active,name)

    ! Print Inactive tally
    call self % inactiveTally % print(out)

    ! Print Active attachment
    call self % activeAtch % print(out)

    call self % activeTally % print(out)

    path = trim(self % outputFile) // '.m'
    call out % writeToFile(path)

  end subroutine collectResults


  !!
  !! Initialise from individual components and dictionaries for inactive and active tally
  !!
  subroutine init(self, dict)
    class(eigenPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)          :: dict
    class(dictionary),pointer                 :: tempDict
    type(dictionary)                          :: locDict1, locDict2
    integer(shortInt)                         :: seed_temp
    integer(longInt)                          :: seed
    character(10)                             :: time
    character(8)                              :: date
    character(:),allocatable                  :: string
    character(nameLen)                        :: nucData, energy
    type(visualiser)                          :: viz
    class(geometry), pointer                  :: geom
    integer(shortInt)                         :: i
    character(100), parameter :: Here ='init (eigenPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_inactive,'inactive')
    call dict % get( self % N_active,'active')
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
    self % geom => new_cellGeometry_ptr(tempDict, ndReg_getMatNames())

    ! Activate Nuclear Data *** All materials are active
    call ndReg_activate(self % particleType, nucData, [(i, i=1, mm_nMat())])
    self % nucData => ndReg_get(self % particleType)

    ! Call visualisation
    if (dict % isPresent('viz')) then
      print *, "Initialising visualiser"
      tempDict => dict % getDictPtr('viz')
      geom => self % geom
      call viz % init(geom, tempDict)
      print *, "Constructing visualisation"
      call viz % makeViz()
      call viz % kill()
    endif

    ! Build collision operator
    tempDict => dict % getDictPtr('collisionOperator')
    call self % collOp % init(tempDict)

    ! Build transport operator
    tempDict => dict % getDictPtr('transportOperator')
    call new_transportOperator(self % transOp, tempDict, self % geom)

    ! Initialise active & inactive tally Admins
    tempDict => dict % getDictPtr('inactiveTally')
    allocate(self % inactiveTally)
    call self % inactiveTally % init(tempDict)

    tempDict => dict % getDictPtr('activeTally')
    allocate(self % activeTally)
    call self % activeTally % init(tempDict)

    ! Load Initial source
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
    call locDict1 % init(2)
    call locDict2 % init(2)

    call locDict2 % store('type','keffAnalogClerk')
    call locDict1 % store('keff', locDict2)
    call locDict1 % store('display',['keff'])

    allocate(self % inactiveAtch)
    call self % inactiveAtch % init(locDict1)

    call locDict2 % kill()
    call locDict1 % kill()

    ! Active tally attachment
    call locDict1 % init(2)
    call locDict2 % init(2)

    call locDict2 % store('type','keffImplicitClerk')
    call locDict1 % store('keff', locDict2)
    call locDict1 % store('display',['keff'])

    allocate(self % activeAtch)
    call self % activeAtch % init(locDict1)

    call locDict2 % kill()
    call locDict1 % kill()

    ! Attach attachments to result tallies
    call self % inactiveTally % push(self % inactiveAtch)
    call self % activeTally % push(self % activeAtch)


    call self % printSettings()

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(eigenPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(eigenPhysicsPackage), intent(in) :: self

    print *, repeat("<>",50)
    print *, "/\/\ EIGENVALUE CALCULATION WITH POWER ITERATION METHOD /\/\"
    print *, "Inactive Cycles:    ", numToChar(self % N_inactive)
    print *, "Active Cycles:      ", numToChar(self % N_active)
    print *, "Neutron Population: ", numToChar(self % pop)
    print *, "Initial RNG Seed:   ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings

end module eigenPhysicsPackage_class
