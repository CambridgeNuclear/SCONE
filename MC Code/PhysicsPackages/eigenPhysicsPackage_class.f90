module eigenPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,             only : fatalError
  use dictionary_class,              only : dictionary

  ! Particle classes and Random number generator
  use particle_class,                only : particle, phaseCoord
  use particleDungeon_class,         only : particleDungeon
  use RNG_class,                     only : RNG

  ! Geometry & Nuclear Data
  use cellGeometry_inter,            only : cellGeometry
  use basicCellCSG_class,            only : basicCellCSG !** Provisional
  use nuclearData_inter,             only : nuclearData

  ! Operators
  use collisionOperatorBase_inter,   only : collisionOperatorBase
  use transportOperator_inter,       only : transportOperator

  ! Tallies
  use tallyInactiveAdmin_class,      only : tallyInactiveAdmin
  use tallyActiveAdmin_class,        only : tallyActiveAdmin

  ! Factories
  use nuclearDataFactory_func,       only : new_nuclearData_ptr
  use collisionOperatorFactory_func, only : new_collisionOperator_ptr
  use transportOperatorFactory_func, only : new_transportOperator_ptr

  implicit none
  private

  !!
  !! Physics Package for eigenvalue calculations
  !!
  type, public :: eigenPhysicsPackage
    ! private ** DEBUG
    ! Building blocks
    class(nuclearData), pointer            :: nucData       => null()
    class(cellGeometry), pointer           :: geom          => null()
    class(collisionOperatorBase), pointer  :: collOp        => null()
    class(transportOperator), pointer      :: transOp       => null()
    class(RNG), pointer                    :: pRNG          => null()
    type(tallyInactiveAdmin),pointer       :: inactiveTally => null()
    type(tallyActiveAdmin),pointer         :: activeTally   => null()

    ! Settings
    integer(shortInt)  :: N_inactive = 300
    integer(shortInt)  :: N_active   = 500
    integer(shortInt)  :: pop        = 5000

    ! Calculation components
    type(particleDungeon), pointer :: thisCycle    => null()
    type(particleDungeon), pointer :: nextCycle    => null()
    type(particleDungeon), pointer :: temp_dungeon => null()

  contains
    procedure :: init
    procedure :: inactiveCycles
    procedure :: activeCycles
    procedure :: generateInitialState
    procedure :: run

  end type eigenPhysicsPackage

contains

  subroutine run(self)
    class(eigenPhysicsPackage), intent(inout) :: self

    call self % generateInitialState()
    call self % inactiveCycles()
    call self % activeCycles()

  end subroutine


  !!
  !! Perform inactive cycles
  !!
  subroutine inactiveCycles(self)
    class(eigenPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    type(phaseCoord)                          :: preState
    integer(shortInt)                         :: i, Nstart, Nend
    real(defReal) :: k_old, k_new

    ! Attach nuclear data and RNG to neutron
    neutron % xsData => self % nucData
    neutron % pRNG   => self % pRNG

    ! Attach tally to operators
    self % collOP % tally => self % inactiveTally

    do i=1,self % N_inactive
      ! Send start of cycle report
      Nstart = self % thisCycle % popSize()
      call self % inactiveTally % reportCycleStart(self % thisCycle)

      gen: do

        ! Obtain paticle from current cycle dungeon
        call self % thisCycle % release(neutron)
        call self % geom % placeCoord(neutron % coords)

          history: do
            preState = neutron

            call self % transOp % transport(neutron)

            ! Exit history and score leakage
            if(neutron % isDead) then
              call self % inactiveTally %  reportHist(preState,neutron,5001)
              exit history
            end if

            call self % collOp % collide(neutron, self % thisCycle, self % nextCycle)

            ! Exit history and score absorbtion
            if(neutron % isDead) then
              call self % inactiveTally %  reportHist(preState,neutron,5000)
              exit history
            end if

          end do history

        if( self % thisCycle % isEmpty()) exit gen
      end do gen

      ! Send end of cycle report
      Nend = self % nextCycle % popSize()
      call self % inactiveTally % reportCycleEnd(self % nextCycle)

      ! Normalise population
      call self % nextCycle % normSize(self % pop, neutron % pRNG)

      ! Flip cycle dungeons
      self % temp_dungeon => self % nextCycle
      self % nextCycle    => self % thisCycle
      self % thisCycle    => self % temp_dungeon

      ! Load new k-eff estimate into next cycle dungeon
      k_old = self % nextCycle % k_eff
      k_new = self % inactiveTally % keff()

      self % nextCycle % k_eff = k_new

      ! Display progress
      print *, 'Cycle: ', i, ' of ', self % N_inactive,' Pop: ', Nstart, ' -> ',Nend
      call self % inactiveTally % display()
    end do

  end subroutine inactiveCycles

  !!
  !! Perform active cycles
  !!
  subroutine activeCycles(self)
    class(eigenPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    type(phaseCoord)                          :: preState
    integer(shortInt)                         :: i, Nstart, Nend
    real(defReal) :: k_old, k_new

    ! Attach nuclear data and RNG to neutron
    neutron % xsData => self % nucData
    neutron % pRNG   => self % pRNG

    ! Attach tally to operators
    self % collOP % tally => self % activeTally

    do i=1,self % N_active
      ! Send start of cycle report
      Nstart = self % thisCycle % popSize()
      call self % activeTally % reportCycleStart(self % thisCycle)

      gen: do

        ! Obtain paticle from current cycle dungeon
        call self % thisCycle % release(neutron)
        call self % geom % placeCoord(neutron % coords)

          history: do
            preState = neutron
            call self % transOp % transport(neutron)

            ! Exit history and score leakage
            if(neutron % isDead) then
              call self % activeTally %  reportHist(preState,neutron,5001)
              exit history
            end if


            call self % collOp % collide(neutron, self % thisCycle, self % nextCycle)

            ! Exit history and score absorbtion
            if(neutron % isDead) then
              call self % activeTally %  reportHist(preState,neutron,5000)
              exit history
            end if

          end do history

        if( self % thisCycle % isEmpty()) exit gen
      end do gen

      ! Send end of cycle report
      Nend = self % nextCycle % popSize()
      call self % activeTally % reportCycleEnd(self % nextCycle)

      ! Normalise population
      call self % nextCycle % normSize(self % pop, neutron % pRNG)

      ! Flip cycle dungeons
      self % temp_dungeon => self % nextCycle
      self % nextCycle    => self % thisCycle
      self % thisCycle    => self % temp_dungeon

      ! Load new k-eff estimate into next cycle dungeon
      k_old = self % nextCycle % k_eff
      k_new = self % activeTally % keff()

      self % nextCycle % k_eff = k_new

      ! Display progress
      print *, 'Cycle: ', i, ' of ', self % N_active,' Pop: ', Nstart, ' -> ',Nend
      call self % activeTally % display()
    end do
  end subroutine activeCycles

  subroutine generateInitialState(self)
    class(eigenPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    integer(shortInt)                         :: i

    allocate(self % thisCycle)
    allocate(self % nextCycle)

    call self % thisCycle % init(3*self % pop)
    call self % nextCycle % init(3*self % pop)

     do i=1,self % pop
      neutron % E      = 0.5
      !neutron % G      = 1
      call neutron % teleport([0.1_8, 0.1_8, 0.1_8])
      call neutron % point([1.0_8, 0.0_8, 0.0_8])
      neutron % w      = 1.0
      neutron % isDead = .false.
      neutron % isMG   = .false.
      call self % thisCycle % detain(neutron)
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
    class(dictionary), intent(in)             :: inactiveDict
    class(dictionary), intent(in)             :: activeDict

    ! Build nuclear data
    self % nucData => new_nuclearData_ptr(matDict)

    ! Build geometry *** Will be replaced by a factory SOON !!!
    allocate(basicCellCSG :: self % geom)
    associate ( g => self % geom)
      select type (g)
        type is ( basicCellCSG)
          call g % init( geomDict, self % nucData)
      end select
    end associate

    ! Build collision operator
    self % collOp => new_collisionOperator_ptr(self % nucData, collDict)

    ! Build transport operator *** Maybe put dictionary at the end in geometry as well?
    self % transOp => new_transportOperator_ptr(self % nucData, self % geom, transDict)

    ! Initialise active & inactive tally Admins
    allocate(self % inactiveTally)
    call self % inactiveTally % init(inactiveDict)

    allocate(self % activeTally)
    call self % activeTally % init(activeDict)

    ! Initialise RNG
    allocate(self % pRNG)
    call self % pRNG % init(768568_8)
    !call self % pRNG % init(67858567567_8)


  end subroutine init


    
end module eigenPhysicsPackage_class
