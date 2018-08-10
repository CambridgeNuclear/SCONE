module eigenPhysicsPackage_class

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
  use tallyInactiveAdmin_class,       only : tallyInactiveAdmin
  use tallyActiveAdmin_class,         only : tallyActiveAdmin

  ! Factories
  use nuclearDataFactory_func,        only : new_nuclearData_ptr
  use geometryFactory_func,           only : new_cellGeometry_ptr
  use collisionOperatorFactory_func,  only : new_collisionOperator_ptr
  use transportOperatorFactory_func,  only : new_transportOperator_ptr

  implicit none
  private

  !!
  !! Physics Package for eigenvalue calculations
  !!
  type, public,extends(physicsPackage) :: eigenPhysicsPackage
     private
    ! Building blocks
    class(nuclearData), pointer            :: nucData       => null()
    class(transportNuclearData), pointer   :: transNucData  => null()
    class(cellGeometry), pointer           :: geom          => null()
    class(collisionOperatorBase), pointer  :: collOp        => null()
    class(transportOperator), pointer      :: transOp       => null()
    class(RNG), pointer                    :: pRNG          => null()
    class(tallyInactiveAdmin),pointer      :: inactiveTally => null()
    class(tallyActiveAdmin),pointer        :: activeTally   => null()

    ! Settings
    integer(shortInt)  :: N_inactive
    integer(shortInt)  :: N_active
    integer(shortInt)  :: pop
    character(pathLen) :: outputFile

    ! Calculation components
    type(particleDungeon), pointer :: thisCycle    => null()
    type(particleDungeon), pointer :: nextCycle    => null()
    type(particleDungeon), pointer :: temp_dungeon => null()

  contains
    procedure :: init
    procedure :: printSettings
    procedure :: inactiveCycles
    procedure :: activeCycles
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
    call self % inactiveCycles()
    call self % activeCycles()
    call self % collectResults()

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
      call printFishLineR(i)
      print *
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
      call printFishLineR(i)
      print *
      print *, 'Cycle: ', i, ' of ', self % N_active,' Pop: ', Nstart, ' -> ',Nend
      call self % activeTally % display()

    end do
  end subroutine activeCycles

  subroutine generateInitialState(self)
    class(eigenPhysicsPackage), intent(inout) :: self
    type(particle)                            :: neutron
    integer(shortInt)                         :: i, matIdx, dummy
    real(defReal),dimension(6)                :: bounds
    real(defReal),dimension(3)                :: top, bottom
    real(defReal),dimension(3)                :: r
    real(defReal),dimension(3)                :: rand
    character(100), parameter :: Here =' generateInitialState( eigenPhysicsPackage_class.f90)'

    ! Allocate and initialise particle Dungeons
    allocate(self % thisCycle)
    allocate(self % nextCycle)

    call self % thisCycle % init(3*self % pop)
    call self % nextCycle % init(3*self % pop)

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
      call self % thisCycle % detain(neutron)

      ! Update iterator
      i = i +1
    end do
    print *, "DONE!"

  end subroutine generateInitialState


  !!
  !! Print calculation results to file
  !!
  subroutine collectResults(self)
    class(eigenPhysicsPackage), intent(in) :: self
    type(outputFile)                       :: out
    character(pathLen)                     :: path
    character(nameLen)                     :: name

    name = 'asciiMATLAB'
    call out % init(name)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Inactive_Cycles'
    call out % printValue(self % N_inactive,name)

    name = 'Active_Cycles'
    call out % printValue(self % N_active,name)

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
    type(dictionary)                          :: tempDict
    integer(shortInt)                         :: seed_temp
    integer(longInt)                          :: seed
    character(10)                             :: time
    character(8)                              :: date
    character(:),allocatable                  :: string
    class(nuclearData),pointer                :: nucData_ptr
    character(100), parameter :: Here ='init (eigenPhysicsPackage_class.f90)'

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_inactive,'inactive')
    call dict % get( self % N_active,'active')

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')
    !self % outputFile = string

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

    ! Initialise active & inactive tally Admins
    call dict % get(tempDict,'inactiveTally')
    allocate(self % inactiveTally)
    call self % inactiveTally % init(tempDict)

    call dict % get(tempDict,'activeTally')
    allocate(self % activeTally)
    call self % activeTally % init(tempDict)


   ! call self % pRNG % init(768568_8)
    !call self % pRNG % init(6585886547_8)
    !call self % pRNG % init(67858567567_8)
    !call self % pRNG % init(5764746_8)
    !call self % pRNG % init(1346575219654672_8)

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
