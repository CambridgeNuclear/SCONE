module dynamPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,             only : fatalError
  use dictionary_class,              only : dictionary

  ! Particle classes and Random number generator
  use particle_class,                only : particle, phaseCoord
  use particleDungeon_class,         only : particleDungeon
  use RNG_class,                     only : RNG

  ! Geometry & Nuclear Data
  use geometry_class,                only : geometry
  use nuclearData_inter,             only : nuclearData

  ! Operators
  use collisionOperatorBase_inter,   only : collisionOperatorBase
  use transportOperator_inter,       only : transportOperator

  ! Tallies
  use tallyActiveAdmin_class,        only : tallyActiveAdmin

  ! Factories
  use nuclearDataFactory_func,       only : new_nuclearData_ptr
  use collisionOperatorFactory_func, only : new_collisionOperator_ptr
  use transportOperatorFactory_func, only : new_transportOperator_ptr

  implicit none
  private

  !!
  !! Physics Package for dynamic calculations
  !!
  type, public :: dynamPhysicsPackage
    ! private ** DEBUG
    ! Building blocks
    class(nuclearData), pointer            :: nucData       => null()
    class(geometry), pointer               :: geom          => null()
    class(collisionOperatorBase), pointer  :: collOp        => null()
    class(transportOperator), pointer      :: transOp       => null()
    class(RNG), pointer                    :: pRNG          => null()
    type(tallyActiveAdmin),pointer         :: timeTally     => null()

    ! Settings
    integer(shortInt)  :: N_steps    = 500
    integer(shortInt)  :: pop        = 5000
    real(defReal)      :: timeMax    = ZERO

    ! Calculation components
    type(particleDungeon), pointer :: thisStep     => null()
    type(particleDungeon), pointer :: nextStep     => null()
    type(particleDungeon), pointer :: temp_dungeon => null()

  contains
    procedure :: init
    procedure :: timeSteps
    procedure :: generateInitialState
    procedure :: run

  end type dynamPhysicsPackage

contains

  subroutine run(self)
    class(dynamPhysicsPackage), intent(inout) :: self

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
    real(defReal) :: power_old, power_new

    ! Attach nuclear data and RNG to neutron
    neutron % xsData => self % nucData
    neutron % pRNG   => self % pRNG

    ! Attach tally to operators
    self % collOP % tally => self % timeTally

    ! Load initial power
    power_old = self % timeTally % power()

    do i=1,self % N_steps
      ! Send start of cycle report
      Nstart = self % thisStep % popSize()
      call self % timeTally % reportCycleStart(self % thisStep)

      ! Increment timeMax
      self % timeMax = self % timeMax + self % timeTally % stepLengths(i)

      step: do

        ! Obtain paticle from current cycle dungeon
        call self % thisStep % release(neutron)
        call self % geom % placeParticle(neutron)
        neutron % timeMax = self % timeMax  ! Need to define where this comes from!!!

          history: do
            preState = neutron
            call self % transOp % transport(neutron)

            ! Exit history - score leakage or place particle in next time point
            if(neutron % isDead) then
              call self % timeTally % reportHist(preState,neutron,neutron%fate)
              ! Revive neutron and place in next step if it became too old
              if (neutron % fate == aged_FATE) then
                neutron % isDead = .false.
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
      call self % timeTally % reportCycleEnd(self % nextStep)
      power_new = self % timeTally % power()

      ! Normalise population
      call self % nextStep % normSize(self % pop, neutron % pRNG)

      ! Adjust weight for population growth or decay
      ! Only applied after the first step - first step is used for subsequent normalisation
      if (i > 1) then
        call self % nextStep % normWeight(thisStep % popWgt * power_new / power_old)
      end if

      ! Flip cycle dungeons
      self % temp_dungeon => self % nextStep
      self % nextStep     => self % thisStep
      self % thisStep     => self % temp_dungeon

      power_old = power_new

      ! Display progress
      print *, 'Step: ', i, ' of ', self % N_steps,' Pop: ', Nstart, ' -> ',Nend
      call self % timeTally % display()
    end do
  end subroutine timeSteps

  subroutine generateInitialState(self)
    class(dynamPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    integer(shortInt)                         :: i

    allocate(self % thisStep)
    allocate(self % nextStep)

    call self % thisStep % init(3*self % pop)
    call self % nextStep % init(3*self % pop)

     do i=1,self % pop
      neutron % E      = ONE
      !neutron % G      = 1
      call neutron % teleport([ZERO, ZERO, ZERO])
      call neutron % point([1.0_8, 0.0_8, 0.0_8]) ! Should distribute uniformly in angle
      neutron % w      = ONE
      neutron % isDead = .false.
      neutron % isMG   = .false.
      neutron % time   = ZERO
      call self % thisStep % detain(neutron)
    end do

  end subroutine generateInitialState


  !!
  !! Initialise from individual components and dictionaries for inactive and active tally
  !!
  subroutine init(self, matDict, geomDict, collDict ,transDict, inactiveDict, activeDict )
    class(eigenPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(in)             :: matDict
    class(dictionary), intent(inout)          :: geomDict  !*** Maybe change to intent(in) in geom?
    class(dictionary), intent(in)             :: collDict
    class(dictionary), intent(in)             :: transDict
    class(dictionary), intent(in)             :: timeDict

    ! Build nuclear data
    self % nucData => new_nuclearData_ptr(matDict)

    ! Build geometry *** Will be replaced by a factory at some point
    allocate(self % geom)
    call self % geom % init( geomDict, self % nucData)

    ! Build collision operator
    self % collOp => new_collisionOperator_ptr(self % nucData, collDict)

    ! Build transport operator *** Maybe put dictionary at the end in geometry as well?
    self % transOp => new_transportOperator_ptr(self % nucData, self % geom, transDict)

    allocate(self % activeTally)
    call self % timeTally % init(timeDict)

    ! Initialise RNG
    allocate(self % pRNG)
    call self % pRNG % init(768568_8)
    !call self % pRNG % init(67858567567_8)
  end subroutine init
    
end module dynamPhysicsPackage_class
