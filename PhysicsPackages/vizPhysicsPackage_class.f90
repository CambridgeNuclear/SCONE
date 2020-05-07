module vizPhysicsPackage_class

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
  use neutronMaterial_inter,          only : neutronMaterial, neutronMaterial_CptrCast
  use ceNeutronMaterial_class,        only : ceNeutronMaterial
  use mgNeutronMaterial_inter,        only : mgNeutronMaterial
  use fissionCE_class,                only : fissionCE, fissionCE_TptrCast
  use fissionMG_class,                only : fissionMG, fissionMG_TptrCast

  ! Factories
  use geometryFactory_func,           only : new_cellGeometry_ptr

  ! Visualisation
  use visualiser_class,               only : visualiser

  implicit none
  private

  !!
  !! Physics Package for eigenvalue calculations
  !!
  type, public,extends(physicsPackage) :: vizPhysicsPackage
    private
    ! Building blocks
    class(nuclearDatabase), pointer        :: nucData       => null()
    class(cellGeometry), pointer           :: geom          => null()
    class(RNG), pointer                    :: pRNG          => null()

    ! Settings
    character(pathLen) :: outputFile
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: particleType

    ! Timer bins
    integer(shortInt) :: timerMain

  contains
    procedure :: init
    procedure :: run
    procedure :: kill

  end type vizPhysicsPackage

contains

  !!
  !! Not implemented - visualisation occurs during init for now
  !!
  subroutine run(self)
    class(vizPhysicsPackage), intent(inout) :: self

  end subroutine

  !!
  !! Initialise from individual components and dictionaries
  !!
  subroutine init(self, dict)
    class(vizPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)        :: dict
    class(dictionary),pointer               :: tempDict
    integer(shortInt)                       :: seed_temp
    integer(longInt)                        :: seed
    character(10)                           :: time
    character(8)                            :: date
    character(:),allocatable                :: string
    character(nameLen)                      :: nucData, energy
    type(visualiser)                        :: viz
    class(geometry), pointer                :: geom
    integer(shortInt)                       :: i
    character(100), parameter :: Here ='init (vizPhysicsPackage_class.f90)'

    ! Read calculation settings
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
      print *, "CONSTRUCTING VISUALISATION"
      tempDict => dict % getDictPtr('viz')
      geom => self % geom
      call viz % init(geom, tempDict)
    else
      call fatalError(here,'Must provide viz dict for plotting')
    endif

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(vizPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

end module vizPhysicsPackage_class
