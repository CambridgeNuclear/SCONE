module eigenPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, printFishLineR, numToChar
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Particle classes and Random number generator
  use particle_class,                 only : particle, phaseCoord, P_NEUTRON
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
  use collisionOperator_class,        only : collisionOperator
  use transportOperator_inter,        only : transportOperator

  ! Tallies
  use tallyCodes
  use tallyAdmin_class,               only : tallyAdmin
  use tallyResult_class,              only : tallyResult
  use keffAnalogClerk_class,          only : keffResult

  ! Factories
  use nuclearDataRegistry_mod,        only : build_NuclearData, getHandlePtr
  use geometryFactory_func,           only : new_cellGeometry_ptr
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
    type(collisionOperator)                :: collOp
    class(transportOperator), pointer      :: transOp       => null()
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

    ! Calculation components
    type(particleDungeon), pointer :: thisCycle    => null()
    type(particleDungeon), pointer :: nextCycle    => null()
    type(particleDungeon), pointer :: temp_dungeon => null()

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
    character(100),parameter :: Here ='cycles (eigenPhysicsPackage_class.f90)'

    ! Attach nuclear data and RNG to neutron
    neutron % xsData => self % nucData
    neutron % pRNG   => self % pRNG

    ! Attach tally to operators
    self % transOP % tally => tally

    ! Set initiial k-eff
    k_new = ONE

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

          history: do
            call neutron % savePreHistory()

            call self % transOp % transport(neutron)

            ! Exit history and score leakage
            if(neutron % isDead) then
              neutron % fate = leak_FATE
              call tally %  reportHist(neutron)
              exit history
            end if

            call self % collOp % collide(neutron, tally ,self % thisCycle, self % nextCycle)

            ! Exit history and score absorbtion
            if(neutron % isDead) then
              neutron % fate = abs_FATE
              call tally %  reportHist(neutron)
              exit history
            end if

          end do history

        if( self % thisCycle % isEmpty()) exit gen
      end do gen

      ! Send end of cycle report
      Nend = self % nextCycle % popSize()
      call tally % reportCycleEnd(self % nextCycle)

      ! Normalise population
      call self % nextCycle % normSize(self % pop, neutron % pRNG)

      if(self % printSource == 1) then
        call self % nextCycle % printSourceToFile(trim(self % outputFile)//'_source'//numToChar(i))
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

      ! Display progress
      call printFishLineR(i)
      print *
      print *, 'Cycle: ', i, ' of ', N_cycles,' Pop: ', Nstart, ' -> ',Nend
      call tally % display()
    end do
  end subroutine cycles

  !!
  !!
  !!
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

      ! Assign type of the particle
      neutron % type = P_NEUTRON

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
    class(dictionary),pointer                 :: tempDict
    type(dictionary)                          :: locDict1, locDict2
    integer(shortInt)                         :: seed_temp
    integer(longInt)                          :: seed
    character(10)                             :: time
    character(8)                              :: date
    character(:),allocatable                  :: string
    character(nameLen)                        :: nucData
    class(nuclearData),pointer                :: nucData_ptr
    character(100), parameter :: Here ='init (eigenPhysicsPackage_class.f90)'

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_inactive,'inactive')
    call dict % get( self % N_active,'active')
    call dict % get( nucData, 'XSdata')

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

    ! Read whether to print particle source per cycle
    call dict % getOrDefault(self % printSource, 'printSource', 0)

    ! Build nuclear data
    tempDict => dict % getDictPtr('nuclearData')
    call build_nuclearData(tempDict)
    nucData_ptr => getHandlePtr(nucData)
    self % nucData => nucData_ptr

    ! Attach transport nuclear data
    select type(nucData_ptr)
      class is(transportNuclearData)
        self % transNucData => nucData_ptr

      class default
        call fatalError(Here,'Nuclear data needs be of class: transportNuclearData')

    end select

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    self % geom => new_cellGeometry_ptr(tempDict, self % nucData)

    ! Build collision operator
    tempDict => dict % getDictPtr('collisionOperator')
    call self % collOp % init(tempDict)

    ! Build transport operator
    tempDict => dict % getDictPtr('transportOperator')
    self % transOp => new_transportOperator_ptr(self % nucData, self % geom, tempDict)

    ! Initialise active & inactive tally Admins
    tempDict => dict % getDictPtr('inactiveTally')
    allocate(self % inactiveTally)
    call self % inactiveTally % init(tempDict)

    tempDict => dict % getDictPtr('activeTally')
    allocate(self % activeTally)
    call self % activeTally % init(tempDict)

    ! Initialise active and inactive tally attachments

    ! Inactive tally attachment
    call locDict1 % init(2)
    call locDict2 % init(2)

    call locDict2 % store('type','keffAnalogClerk')
    call locDict2 % store('display','yes')
    call locDict1 % store('keff', locDict2)

    allocate(self % inactiveAtch)
    call self % inactiveAtch % init(locDict1)

    call locDict2 % kill()
    call locDict1 % kill()

    ! Active tally attachment
    call locDict1 % init(2)
    call locDict2 % init(2)

    call locDict2 % store('type','keffImplicitClerk')
    call locDict2 % store('display','yes')
    call locDict1 % store('keff', locDict2)

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
