module dynamPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, printFishLineR, numToChar
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary

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
    class(nuclearData), pointer            :: nucData       => null()
    class(transportNuclearData), pointer   :: transNucData  => null()
    class(cellGeometry), pointer           :: geom          => null()
    class(collisionOperatorBase), pointer  :: collOp        => null()
    class(transportOperator), pointer      :: transOp       => null()
    class(RNG), pointer                    :: pRNG          => null()
    type(tallyTimeAdmin), pointer          :: timeTally     => null()

    ! Settings
    integer(shortInt)                        :: N_steps
    integer(shortInt)                        :: pop
    real(defReal), dimension(:), allocatable :: stepLength

    ! Calculation components
    type(particleDungeon), pointer :: thisStep     => null()
    type(particleDungeon), pointer :: nextStep     => null()
    type(particleDungeon), pointer :: temp_dungeon => null()

  contains
    procedure :: init
    procedure :: printSettings
    procedure :: timeSteps
    procedure :: generateInitialState
    procedure :: run

  end type dynamPhysicsPackage

contains

  subroutine run(self)
    class(dynamPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ DYNAMIC CALCULATION /\/\" 

    call self % generateInitialState()
    call self % timeSteps()

  end subroutine

  !!
  !! Perform transport over time steps
  !!
  subroutine timeSteps(self)
    class(dynamPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    type(phaseCoord)                          :: preState
    integer(shortInt)                         :: i, Nstart, Nend
    real(defReal)                             :: power_old, power_new, timeMax, genWeight

    ! Attach nuclear data and RNG to neutron
    neutron % xsData => self % nucData
    neutron % pRNG   => self % pRNG

    ! Attach tally to operators
    self % collOP % tally => self % timeTally

    ! Load initial power
    power_old = self % timeTally % power()
    print *,'POWER OLD = ',power_old

    ! Initialise time max
    timeMax = ZERO

    print *,'BEGINNING TIME STEPPING'

    do i=1,self % N_steps
      ! Send start of cycle report
      Nstart = self % thisStep % popSize()
      genWeight = self % thisStep % popWeight()
      call self % timeTally % reportCycleStart(self % thisStep)

      ! Increment timeMax
      timeMax = timeMax + self % stepLength(i)

      step: do

        ! Obtain paticle from current cycle dungeon
        call self % thisStep % release(neutron)
        call self % geom % placeCoord(neutron % coords)
        neutron % timeMax = timeMax

          history: do
            preState = neutron
            call self % transOp % transport(neutron)

            ! Exit history - score leakage or place particle in next time point
            if(neutron % isDead) then
              call self % timeTally % reportHist(preState,neutron,neutron%fate)
              ! Revive neutron and place in next step if it became too old
              if (neutron % fate == aged_FATE) then
                neutron % isDead = .false.
                neutron % fate = 0
                call self % nextStep % detain(neutron)
              end if
              exit history
            end if

            call self % collOp % collide(neutron, self % thisStep, self % thisStep)

            ! Exit history and score absorbtion
            if(neutron % isDead) then
              call self % timeTally % reportHist(preState,neutron,abs_FATE)
              exit history
            end if

          end do history

        if( self % thisStep % isEmpty()) exit step
      end do step

      ! Send end of cycle report
      Nend = self % nextStep % popSize()
      ! Must obtain power before ending cycle! Otherwise, zero power will be returned
      power_new = self % timeTally % power()
      call self % timeTally % reportCycleEnd(self % nextStep)

      ! Normalise population
      call self % nextStep % normSize(self % pop, neutron % pRNG)

      ! Adjust weight for population growth or decay
      ! Only applied after the first step - first step is used for subsequent normalisation
      if (i > 1) then
        call self % nextStep % normWeight(genWeight * power_new / power_old)
      end if
      print *,'POPWEIGHT = ',genWeight

      ! Flip cycle dungeons
      self % temp_dungeon => self % nextStep
      self % nextStep     => self % thisStep
      self % thisStep     => self % temp_dungeon

      power_old = power_new

      ! Display progress
      call printFishLineR(i)
      print *
      print *, 'Time: ', timeMax
      print *, 'Step: ', i, ' of ', self % N_steps,' Pop: ', Nstart, ' -> ',Nend
      call self % timeTally % display()

    end do
  end subroutine timeSteps

  subroutine generateInitialState(self)
    class(dynamPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    integer(shortInt)                         :: i
    real(defReal)                             :: azim, polar

    allocate(self % thisStep)
    allocate(self % nextStep)

    call self % thisStep % init(3*self % pop)
    call self % nextStep % init(3*self % pop)

    print *, "GENERATING ISOTROPIC NEUTRON SOURCE"
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
    print *, "NEUTRONS GENERATED"

  end subroutine generateInitialState


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
    call dict % get(tempDict,'transportOperator')
    self % transOp => new_transportOperator_ptr(self % nucData, self % geom, tempDict)

    ! Initialise time tally Admins
    call dict % get(tempDict,'timeTally')
    allocate(self % timeTally)
    call self % timeTally % init(tempDict)

    call self % printSettings()

  end subroutine init

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettings(self)
    class(dynamPhysicsPackage), intent(in) :: self

    print *, repeat("<>",50)
    print *, "/\/\ DYNAMIC CALCULATION USING FIXED SOURCE /\/\" 
    print *, "Number of time steps: ", numToChar(self % N_steps)
    print *, "Neutron Population: ", numToChar(self % pop)
    print *, "Initial RNG Seed:   ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettings
    
end module dynamPhysicsPackage_class
