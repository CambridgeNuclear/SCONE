module dynamPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, printFishLineR, numToChar
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Particle classes and Random number generator
  use particle_class,                 only : particle, phaseCoord
  use particleDungeon_class,          only : particleDungeon
  use RNG_class,                      only : RNG

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry & Nuclear Data
  use cellGeometry_inter,             only : cellGeometry
  use nuclearData_inter,              only : nuclearData
  use transportNuclearData_inter,     only : transportNuclearData
  use perNuclideNuclearDataCE_inter,  only : perNuclideNuclearDataCE
  use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Operators
  use collisionOperatorBase_inter,    only : collisionOperatorBase
  use transportOperator_inter,        only : transportOperator

  ! Tallies
  use tallyCodes
  use tallyTimeAdmin_class,          only : tallyTimeAdmin
  use tallyInactiveAdmin_class,      only : tallyInactiveAdmin

  ! Factories
  use nuclearDataFactory_func,       only : new_nuclearData_ptr
  use geometryFactory_func,          only : new_cellGeometry_ptr
  use collisionOperatorFactory_func, only : new_collisionOperator_ptr
  use transportOperatorFactory_func, only : new_transportOperator_ptr

  implicit none
  private

  !!
  !! Physics Package for dynamic calculations
  !!
  type, public, extends(physicsPackage) :: dynamPhysicsPackage
    private
    ! Building blocks
    class(nuclearData), pointer            :: nucData         => null()
    class(transportNuclearData), pointer   :: transNucData    => null()
    class(cellGeometry), pointer           :: geom            => null()
    class(collisionOperatorBase), pointer  :: collOp          => null()
    class(transportOperator), pointer      :: dynamTransOp    => null()
    class(transportOperator), pointer      :: inactiveTransOp => null()
    class(RNG), pointer                    :: pRNG            => null()
    type(tallyTimeAdmin), pointer          :: timeTally       => null()
    type(tallyInactiveAdmin), pointer      :: inactiveTally   => null()

    ! Settings
    integer(shortInt)                        :: N_steps
    integer(shortInt)                        :: N_inactive
    integer(shortInt)                        :: pop
    real(defReal), dimension(:), allocatable :: stepLength
    real(defReal)                            :: initialStepLength
    logical(defBool)                         :: ageNeutrons = .false.
    logical(defBool)                         :: findFissionSource
    character(pathLen)                       :: outputFile

    ! Calculation components
    type(particleDungeon), pointer :: thisStep     => null()
    type(particleDungeon), pointer :: nextStep     => null()
    type(particleDungeon), pointer :: temp_dungeon => null()

  contains
    procedure :: init
    procedure :: printSettings
    procedure :: inactiveCycles       !! Perform inactive cycles to obtain fundamental mode
    procedure :: timeSteps            !! Perform time stepping
    procedure :: generateInitialState !! Generate initial neutrons for inactive cycles
    procedure :: generatePointSource  !! Generate isotropic point source of neutrons *** Replace in future
    procedure :: ageSourceNeutrons    !! Age initial neutrons to attempt to achieve more realistic distribution
    procedure :: collectResults       !! Print results to output file
    procedure :: run

  end type dynamPhysicsPackage

contains

  subroutine run(self)
    class(dynamPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ DYNAMIC CALCULATION /\/\" 

    ! Obtain an initial neutron source
    if (self % findFissionSource) then
      call self % generateInitialState()
      call self % inactiveCycles()
      if (self % ageNeutrons) call self % ageSourceNeutrons()
    ! Use an isotropic point source
    else
      call self % generatePointSource()
    end if
    call self % timeSteps()
    call self % collectResults()

  end subroutine

  !!
  !! Perform transport over time steps
  !!
  subroutine timeSteps(self)
    class(dynamPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    integer(shortInt)                         :: i, Nstart, Nend
    real(defReal)                             :: power_old, power_new, timeMax, genWeight
    integer(longInt) :: leak

    ! Attach nuclear data and RNG to neutron
    neutron % xsData => self % nucData
    neutron % pRNG   => self % pRNG

    ! Attach tally to operators
    self % collOP % tally => self % timeTally

    ! Load initial power
    power_old = self % timeTally % power()

    ! Initialise time max
    timeMax = ZERO
    leak = 0
    print *,'BEGINNING TIME STEPPING'

    do i=1,self % N_steps
      ! Send start of cycle report
      Nstart = self % thisStep % popSize()
      call self % timeTally % reportCycleStart(self % thisStep)

      ! Increment timeMax
      timeMax = timeMax + self % stepLength(i)

      step: do

        ! Obtain paticle from current cycle dungeon
        call self % thisStep % release(neutron)
        call self % geom % placeCoord(neutron % coords)
        neutron % timeMax = timeMax

          history: do
            call neutron % savePreHistory()
            call self % dynamTransOp % transport(neutron)

            ! Exit history - score leakage or place particle in next time point
            if(neutron % isDead) then
              call self % timeTally % reportHist(neutron)
              ! Revive neutron and place in next step if it became too old
              if (neutron % fate == aged_FATE) then
                neutron % isDead = .false.
                neutron % fate = 0
                call self % nextStep % detain(neutron)
              end if
              if (neutron % fate == leak_FATE) leak = leak + 1
              exit history
            end if

            call self % collOp % collide(neutron, self % thisStep, self % thisStep)

            ! Exit history and score absorbtion
            if(neutron % isDead) then
              neutron % fate = abs_FATE
              call self % timeTally % reportHist(neutron)
              exit history
            end if

          end do history

        if( self % thisStep % isEmpty()) exit step
      end do step

      ! Send end of cycle report
      Nend = self % nextStep % popSize()
      call self % timeTally % reportCycleEnd(self % nextStep)
      power_new = self % timeTally % power()

      ! Store initial weight and normalise population
      genWeight = self % nextStep % popWeight()
      call self % nextStep % normSize(self % pop, neutron % pRNG)

      ! Adjust weight for population growth or decay
      ! Only applied after the first step - first step is used for subsequent normalisation
      !call self % nextStep % normWeight(genWeight * power_new / power_old)
      call self % nextStep % normWeight(genWeight)

      ! Flip cycle dungeons
      self % temp_dungeon => self % nextStep
      self % nextStep     => self % thisStep
      self % thisStep     => self % temp_dungeon

      power_old = power_new

      ! Display progress
      call printFishLineR(i)
      print *
      print '(A,ES15.2,A)', 'Time: ', timeMax,'s'
      print *, 'Step: ', i, ' of ', self % N_steps,' Pop: ', Nstart, ' -> ',Nend
      call self % timeTally % display()
      call self % timeTally % incrementStep()

    end do
  end subroutine timeSteps

  !!
  !! Takes neutrons from inactive cycles and ages them by a chosen time step
  !! Used to obtain a more representative neutron distribution rather than using a fission source
  !!
  subroutine ageSourceNeutrons(self)
    class(dynamPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    integer(shortInt)                         :: i
    real(defReal)                             :: timeMax, genWeight

    ! Attach nuclear data and RNG to neutron
    neutron % xsData => self % nucData
    neutron % pRNG   => self % pRNG

    ! Attach tally to operators
    self % collOP % tally => self % timeTally

    print *,'AGING NEUTRONS WITH AN INITIAL TIME STEP'

    ! Send start of cycle report
    call self % timeTally % reportCycleStart(self % thisStep)

    ! Store initial weight and normalise population
    genWeight = self % thisStep % popWeight()

    ! Increment timeMax
    timeMax = self % initialStepLength

    step: do

      ! Obtain paticle from current cycle dungeon
      call self % thisStep % release(neutron)
      call self % geom % placeCoord(neutron % coords)
      neutron % timeMax = self % initialStepLength

        history: do
          call neutron % savePreHistory()
          call self % dynamTransOp % transport(neutron)

          ! Exit history - score leakage or place particle in next time point
          if(neutron % isDead) then
            call self % timeTally % reportHist(neutron)
            ! Revive neutron and place in next step if it became too old
            if (neutron % fate == aged_FATE) then
              neutron % isDead = .false.
              neutron % fate = 0
              neutron % time = ZERO ! Reset neutron's age
              call self % nextStep % detain(neutron)
            end if
            exit history
          end if

          call self % collOp % collide(neutron, self % thisStep, self % thisStep)

          ! Exit history and score absorbtion
          if(neutron % isDead) then
            exit history
          end if

        end do history

      if( self % thisStep % isEmpty()) exit step
    end do step

    ! Store initial weight and normalise population
    !genWeight = self % nextStep % popWeight()
    call self % nextStep % normSize(self % pop, neutron % pRNG)

    ! Adjust weight for population growth or decay
    call self % nextStep % normWeight(genWeight)

    ! Flip cycle dungeons
    self % temp_dungeon => self % nextStep
    self % nextStep     => self % thisStep
    self % thisStep     => self % temp_dungeon


    ! Display progress
    call printFishLineR(i)
    print *
    print *, 'Initial time step complete'

  end subroutine ageSourceNeutrons

  !!
  !! Perform inactive cycles
  !!
  subroutine inactiveCycles(self)
    class(dynamPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    integer(shortInt)                         :: i, Nstart, Nend
    real(defReal) :: k_old, k_new

    ! Attach nuclear data and RNG to neutron
    neutron % xsData => self % nucData
    neutron % pRNG   => self % pRNG

    ! Attach tally to operators
    self % collOP % tally => self % inactiveTally

    do i=1,self % N_inactive
      ! Send start of cycle report
      Nstart = self % thisStep % popSize()
      call self % inactiveTally % reportCycleStart(self % thisStep)

      gen: do

        ! Obtain paticle from current cycle dungeon
        call self % thisStep % release(neutron)
        call self % geom % placeCoord(neutron % coords)

          history: do
            call neutron % savePreHistory()

            call self % inactiveTransOp % transport(neutron)

            ! Exit history and score leakage
            if(neutron % isDead) then
              neutron % fate = leak_FATE
              call self % inactiveTally %  reportHist(neutron)
              exit history
            end if

            call self % collOp % collide(neutron, self % thisStep, self % nextStep)

            ! Exit history and score absorbtion
            if(neutron % isDead) then
              neutron % fate = abs_FATE
              call self % inactiveTally %  reportHist(neutron)
              exit history
            end if

          end do history

        if( self % thisStep % isEmpty()) exit gen
      end do gen

      ! Send end of cycle report
      Nend = self % nextStep % popSize()
      call self % inactiveTally % reportCycleEnd(self % nextStep)

      ! Normalise population
      call self % nextStep % normSize(self % pop, neutron % pRNG)

      ! Flip cycle dungeons
      self % temp_dungeon => self % nextStep
      self % nextStep     => self % thisStep
      self % thisStep     => self % temp_dungeon

      ! Load new k-eff estimate into next cycle dungeon
      k_old = self % nextStep % k_eff
      k_new = self % inactiveTally % keff()

      self % nextStep % k_eff = k_new

      ! Display progress
      call printFishLineR(i)
      print *
      print *, 'Cycle: ', i, ' of ', self % N_inactive,' Pop: ', Nstart, ' -> ',Nend
      call self % inactiveTally % display()
    end do

    ! Must return k_eff to 1, otherwise population will not shrink or grow
    self % thisStep % k_eff = ONE
    self % nextStep % k_eff = ONE

  end subroutine inactiveCycles

  !!
  !! Generate a source of neutrons with isotropic direction from the origin
  !! ***Will be replaced once fixed source input is implemented
  !!
  subroutine generatePointSource(self)
    class(dynamPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    integer(shortInt)                         :: i
    real(defReal)                             :: azim, polar

    allocate(self % thisStep)
    allocate(self % nextStep)

    call self % thisStep % init(5*self % pop)
    call self % nextStep % init(5*self % pop)

    print *, "GENERATING ISOTROPIC NEUTRON SOURCE AT ORIGIN"
    do i=1,self % pop
      neutron % E      = ONE
      !neutron % G      = 1
      call neutron % teleport([ZERO, ZERO, ZERO])
      ! Sample uniform angular distribution
      polar = TWO * self % pRNG % get() - ONE
      azim  = TWO * self % pRNG % get() - ONE
      call neutron % point([azim*polar, sqrt(ONE - azim*azim)*polar, sqrt(ONE - polar*polar)])
      neutron % w      = ONE
      neutron % isDead = .false.
      neutron % isMG   = .false.
      neutron % time   = ZERO
      call self % thisStep % detain(neutron)
    end do
    print *, "DONE!"

  end subroutine generatePointSource

  !!
  !! Generate initial neutron source for inactive cycles
  !!
  subroutine generateInitialState(self)
    class(dynamPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    integer(shortInt)                         :: i, matIdx, dummy
    real(defReal),dimension(6)                :: bounds
    real(defReal),dimension(3)                :: top, bottom
    real(defReal),dimension(3)                :: r
    real(defReal),dimension(3)                :: rand
    character(100), parameter :: Here =' generateInitialState( dynamPhysicsPackage_class.f90)'

    ! Allocate and initialise particle Dungeons
    allocate(self % thisStep)
    allocate(self % nextStep)

    call self % thisStep % init(5*self % pop)
    call self % nextStep % init(5*self % pop)

    ! Obtain bounds of the geometry
    bounds = self % geom % bounds()

    bottom = bounds([1,3,5])
    top    = bounds([2,4,6])

    ! Initialise iterator and attach RNG to neutron
    i = 0
    neutron % pRNG => self % pRNG

    print *, "GENERATING INITIAL FISSION SOURCE"
    ! Loop over requested population
    do while (i <= self % pop)
      ! Sample position
      rand(1) = neutron % pRNG % get()
      rand(2) = neutron % pRNG % get()
      rand(3) = neutron % pRNG % get()

      r = (top - bottom) * rand + bottom

      ! Find material under postision
      call self % geom % whatIsAt(r,matIdx,dummy)

      ! Resample position if material is not fissile
      if( .not.self % transNucData % isFissileMat(matIdx)) cycle

      ! Put material in neutron
      neutron % coords % matIdx = matIdx

      ! Generate and store fission site
      call self % transNucData % initFissionSite(neutron,r)
      call self % thisStep % detain(neutron)

      ! Update iterator
      i = i +1
    end do
    print *, "DONE!"

  end subroutine generateInitialState

  subroutine collectResults(self)
    class(dynamPhysicsPackage), intent(in) :: self
    type(outputFile)                       :: out
    character(pathLen)                     :: path
    character(nameLen)                     :: name

    name = 'asciiMATLAB'
    call out % init(name)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)
!
!    name = 'pop'
!    call out % printValue(self % pop,name)
!
!    name = 'Inactive_Cycles'
!    call out % printValue(self % N_inactive,name)
!
!    name = 'Active_Cycles'
!    call out % printValue(self % N_active,name)

    call self % timeTally % print(out)

    path = trim(self % outputFile) // '.m'
    call out % writeToFile(path)

  end subroutine collectResults


  !!
  !! Initialise from individual components and dictionaries for inactive and active tally
  !!
  subroutine init(self, dict)
    class(dynamPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)          :: dict
    type(dictionary)                          :: tempDict
    integer(shortInt)                         :: seed_temp
    integer(longInt)                          :: seed
    character(10)                             :: time
    character(8)                              :: date
    character(:),allocatable                  :: string
    class(nuclearData),pointer                :: nucData_ptr
    character(100), parameter :: Here ='init (dynamPhysicsPackage_class.f90)'

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_steps, 'nsteps')
    call dict % get( self % stepLength, 'dt')
    call dict % getOrDefault( self % N_inactive, 'inactive', 0)
    if (self % N_inactive > 0) self % findFissionSource = .true.
    call dict % getOrDefault( self % initialStepLength, 'age', -ONE)
    if (self % initialStepLength < 0) self % ageNeutrons = .false.

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

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

    ! Build nuclear data
    call dict % get(tempDict,'materials')
    nucData_ptr => new_nuclearData_ptr(tempDict)
    self % nucData => nucData_ptr

    ! Attach transport nuclear data
    select type(nucData_ptr)
      class is(transportNuclearData)
        self % transNucData => nucData_ptr

      class default
        call fatalError(Here,'Nuclear data needs be of class: transportNuclearData')

    end select

    ! Build geometry
    call dict % get(tempDict,'geometry')
    self % geom => new_cellGeometry_ptr(tempDict, self % nucData)

    ! Build collision operator
    call dict % get(tempDict,'collisionOperator')
    self % collOp => new_collisionOperator_ptr(self % nucData, tempDict)

    ! Build transport operator
    call dict % get(tempDict,'dynamicTransportOperator')
    self % dynamTransOp => new_transportOperator_ptr(self % nucData, self % geom, tempDict)

    ! Build inactive cycle transport operator if desired
    if (self % findFissionSource) then
      call dict % get(tempDict,'inactiveTransportOperator')
      self % inactiveTransOp => new_transportOperator_ptr(self % nucData, self % geom, tempDict)
    end if

    ! Initialise time tally Admin
    call dict % get(tempDict,'timeTally')
    allocate(self % timeTally)
    call self % timeTally % init(tempDict)

    ! Initialise inactive tally Admin
    call dict % get(tempDict,'inactiveTally')
    allocate(self % inactiveTally)
    call self % inactiveTally % init(tempDict)

    call self % printSettings()

  end subroutine init

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(dynamPhysicsPackage), intent(in) :: self

    print *, repeat("<>",50)
    if (self % findFissionSource) then
      print *, "/\/\ DYNAMIC CALCULATION WITH INACTIVE CYCLES /\/\"
      print *, "Number of inactive cycles: ", numToChar(self % N_inactive)
      if (self % ageNeutrons) then
        print *, "Initial neutron aging step length: ", numToChar(self % initialStepLength)
      end if
    else
      print *, "/\/\ DYNAMIC CALCULATION USING FIXED SOURCE /\/\"
    end if
    print *, "Number of time steps: ", numToChar(self % N_steps)
    print *, "Neutron Population: ", numToChar(self % pop)
    print *, "Initial RNG Seed:   ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings
    
end module dynamPhysicsPackage_class
