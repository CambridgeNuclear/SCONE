module pointVolPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,    only : fatalError, numToChar, rotateVector, printFishLineR
  use hashFunctions_func,   only : FNV_1
  use dictionary_class,     only : dictionary
  use rng_class,            only : RNG
  use physicsPackage_inter, only : physicsPackage

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Geometry
  use coord_class,                    only : coordList
  use geometry_inter,                 only : geometry, distCache
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx, &
                                             gr_kill    => kill
  use geometryFactory_func,           only : new_geometry

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat, mm_matName => matName
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_getMatNames => getMatNames, &
                                             ndReg_kill        => kill

  implicit none
  private

  ! Parameters
  integer(shortInt), parameter :: CSUM = 1, CSUM2 = 2

  !!
  !! Physics package to perform point-sampling-based volume calculation
  !!
  !! Calculates relative volume of different materials in the problem by performing
  !! point-sampling in the geometry. The volume is normalised that the total domain
  !! volume is 1.0.
  !!
  !! Points are sampled uniformly within the geometry bounding box.
  !!
  !! This Physics Package exists to serve as a geometry debugging and benchmarking tool.
  !! Maybe in future it will be useful for burn-up?
  !!
  !! Sample Input Dictionary:
  !!   PP {
  !!     type pointVolPhysicsPackage;
  !!     pop 2000;      // Number of points per cycle
  !!     cycles 100;    // Number of cycles
  !!     #seed 86868;#  // Optional RNG seed
  !!     geometry {<Geometry Definition>}
  !!     nuclearData {<Nuclear data definition. Requires material names only>}
  !!   }
  !!
  !! Private Members
  !!   geom      -> Pointer to the geometry
  !!   geomIdx   -> Index of the geometry in geometry Registry
  !!   rand      -> Random number generator
  !!   timerMain -> Index of the timer defined to measure calculation time
  !!   N_cycles  -> Number of cycles
  !!   res       -> Array to accumulate total number of points in each material. Contains
  !!     cumulative sum over cycles (CSUM) and cumulative sum of squares over
  !!     cycles (CSUM2)
  !!   score     -> Array containing point scores on a given cycle.
  !!
  !! Interface:
  !!   physicsPackage interface
  !!
  type, public, extends(physicsPackage) :: pointVolPhysicsPackage
    private
    ! Components
    class(geometry), pointer :: geom
    integer(shortInt)        :: geomIdx   = 0
    type(RNG)                :: rand
    integer(shortInt)        :: timerMain = 0

    ! Settings
    integer(shortInt)  :: pop       = 0
    integer(shortInt)  :: N_cycles  = 0

    ! Results space
    integer(longInt), dimension(:), allocatable :: score
    real(defReal), dimension(:,:), allocatable :: res

  contains
    ! Superclass procedures
    procedure :: init
    procedure :: run
    procedure :: kill

    ! Private procedures
    procedure, private :: cycles
    procedure, private :: printResults
    procedure, private :: printSettings
  end type pointVolPhysicsPackage

contains

  !!
  !! Initialise Physics Package from dictionary
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine init(self,dict)
    class(pointVolPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)             :: dict
    integer(shortInt)                            :: seed_temp
    integer(longInt)                             :: seed
    character(10)                                :: time
    character(8)                                 :: date
    character(:),allocatable                     :: string
    class(dictionary),pointer                    :: tempDict
    character(nameLen)                           :: geomName
    character(100), parameter :: Here = 'init (pointVolPhysicsPackage_class.f90)'

    ! Load settings
    call dict % get(self % pop, 'pop')
    call dict % get(self % N_cycles, 'cycles')

    ! Check settings
    if (self % pop < 1) call fatalError(Here,'Must sample > 0 points in the geometry.')
    if (self % N_cycles < 1) call fatalError(Here,'Must perform at least > 0 sampling cycles.')

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Initialise RNG
    if( dict % isPresent('seed')) then
      call dict % get(seed_temp,'seed')

    else
      ! Obtain time string and hash it to obtain random seed
      call date_and_time(date, time)
      string = date // time
      call FNV_1(string,seed_temp)

    end if
    seed = seed_temp
    call self % rand % init(seed)

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'pointVolGeom'
    call new_geometry(tempDict, geomName)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Allocate results space
    allocate(self % res(mm_nMat(), 2))
    self % res = ZERO
    allocate(self % score(mm_nMat()))
    self % score = 0

  end subroutine init

  !!
  !! Run calculation
  !!
  !! See physicsPackage_inter for details
  !!
  subroutine run(self)
    class(pointVolPhysicsPackage), intent(inout) :: self

    call self % printSettings()
    call self % cycles(self % rand)
    call self % printResults()

  end subroutine run

  !!
  !! Perform cycles of the stochastic volume calculation with point sampling.
  !!
  !! Randomly places the starting point based on uniform distribution.
  !!
  !! Args:
  !!   rand [inout] -> Initialised random number generator
  !!
  !! NOTE:
  !!   RNG needs to be given as an argument `class(RNG)` to prevent inlining. Compiler (gcc 8.3)
  !!   produced erroneous code without it. Same random number would be produced for different calls
  !!   of `get` function.
  !!
  subroutine cycles(self, rand)
    class(pointVolPhysicsPackage), intent(inout) :: self
    class(RNG), intent(inout)                    :: rand
    real(defReal), dimension(3)                  :: rand3, bottom, top, u
    real(defReal), dimension(3), save            :: r
    integer(shortInt)                            :: gen, point
    integer(shortInt), save                      :: matIdx, uniqueId, i
    type(RNG), save                              :: pRNG
    real(defReal)                                :: elapsed_T, end_T, T_toEnd, cycle_T
    character(100), parameter :: Here = 'cycles (pointVolPhysicsPackage_class.f90)'
    !$omp threadprivate(pRNG, r, matIdx, uniqueId, i)

    !$omp parallel
    pRNG = rand
    !$omp end parallel

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    ! Get lower an upper corner of bounding box
    associate (aabb => self % geom % bounds())
      bottom = aabb(1:3)
      top    = aabb(4:6)
    end associate
        
    u = [ONE, ZERO, ZERO]

    ! Perform clculation
    do gen = 1, self % N_cycles
      !$omp parallel do
      do point = 1, self % pop

        ! Set seed
        call pRNG % stride( (gen-1) * self % pop + point )

        ! Find starting point that is inside the geometry
        i = 0

        rejection : do
          rand3(1) = pRNG % get()
          rand3(2) = pRNG % get()
          rand3(3) = pRNG % get()
          r = bottom + (top - bottom) * rand3

          ! Exit if point is inside the geometry
          call self % geom % whatIsAt(matIdx, uniqueId, r, u)
          if (matIdx /= OUTSIDE_MAT) exit rejection

          i = i + 1
          if (i > 1000) then
            call fatalError(Here, 'Infinite loop when searching for a point in the geometry.')
          end if
        end do rejection
    
        ! Found something    
        if (matIdx /= VOID_MAT) then
          !$omp atomic
          self % score(matIdx) = self % score(matIdx) + 1
        end if

      end do
      !$omp end parallel do

      ! Calculate times
      call timerStop(self % timerMain)
      cycle_T = timerTime(self % timerMain) - elapsed_T
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(self % N_cycles, defReal) * elapsed_T / gen
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      call printFishLineR(gen)
      print *
      print *, 'Cycle: ', numToChar(gen), ' of ', numToChar(self % N_cycles)
      print *, 'Pop: ', numToChar(self % pop)
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))

      ! Process scores
      ! Is this susceptible to round-off?
      self % res(:, CSUM)  = self % res(:, CSUM) + real(self % score(:), defReal) / self % pop
      self % res(:, CSUM2) = self % res(:, CSUM2) + (real(self % score(:), defReal) / self % pop)**2
      self % score = 0

    end do

  end subroutine cycles

  !!
  !! Output calculation results to the console
  !!
  !! Convert cumulative sums to mean and absolute standard deviation and
  !! print them to the console.
  !!
  !! Args:
  !!   None
  !!
  subroutine printResults(self)
    class(pointVolPhysicsPackage), intent(in) :: self
    real(defReal)                             :: mean, SD, var
    integer(shortInt)                         :: i

    print *
    print *, "RELATIVE VOLUME FOR MATERIALS: "
    do i = 1, mm_nMat()
      mean = self % res(i, CSUM) / self % N_cycles
      var = self % res(i, CSUM2) / self % N_cycles - mean**2
      SD = ONE/(self % N_cycles - 1) * sqrt(var)
      print '(A, A, A, ES12.5, A, ES12.5)', " Material: ", mm_matName(i), " Vol", mean, " +/-", SD
    end do

  end subroutine printResults

  !!
  !! Print settings of the point-sampling volume calculation
  !!
  !! Args:
  !!   None
  !!
  subroutine printSettings(self)
    class(pointVolPhysicsPackage), intent(in) :: self

    print *, repeat("<>", MAX_COL/2)
    print *, "/\/\ POINT-SAMPLING RELATIVE VOLUME CALCULATION /\/\"
    print *, "Total Cycles:    ", numToChar(self % N_cycles)
    print *, "Points per cycle: ", numToChar(self % pop)
    print *, "Initial RNG Seed:   ", numToChar(self % rand % getSeed())
    print *
    print *, repeat("<>", MAX_COL/2)

  end subroutine printSettings

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(pointVolPhysicsPackage), intent(inout) :: self

    ! Clean Nuclear Data & Geometries
    call gr_kill()
    call ndreg_kill()

    ! Clean contents
    self % geom    => null()
    self % geomIdx = 0
    !call self % rand % kill()
    self % timerMain = 0

    self % pop      = 0
    self % N_cycles = 0

    if (allocated(self % res)) deallocate(self % res)
    if (allocated(self % score)) deallocate(self % score)

  end subroutine kill

end module pointVolPhysicsPackage_class
