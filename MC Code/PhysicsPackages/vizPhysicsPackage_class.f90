module vizPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, printFishLineR, numToChar
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Particle classes and Random number generator
  use RNG_class,                      only : RNG

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry & Nuclear Data
  use geometry_inter,                 only : geometry
  use cellGeometry_inter,             only : cellGeometry
  use nuclearData_inter,              only : nuclearData
  use transportNuclearData_inter,     only : transportNuclearData
  use perNuclideNuclearDataCE_inter,  only : perNuclideNuclearDataCE
  use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Factories
  use nuclearDataRegistry_mod,        only : build_NuclearData, getHandlePtr
  use geometryFactory_func,           only : new_cellGeometry_ptr

  ! Visualisation
  use visualiser_class,               only : visualiser

  implicit none
  private

  !!
  !! Physics Package purely for visualisation
  !!
  type, public,extends(physicsPackage) :: vizPhysicsPackage
     private
    ! Building blocks
    class(nuclearData), pointer            :: nucData       => null()
    class(transportNuclearData), pointer   :: transNucData  => null()
    class(cellGeometry), pointer           :: geom          => null()
    class(RNG), pointer                    :: pRNG          => null()

    ! Settings
    character(pathLen) :: outputFile

  contains
    procedure :: init
    procedure :: run
    procedure :: kill

  end type vizPhysicsPackage

contains

  !!
  !! Empty procedure that is not called for vizPhysicsPackage
  !! 
  subroutine run(self)
    class(vizPhysicsPackage), intent(inout) :: self

  end subroutine run

  !!
  !! Initialise from individual components and dictionaries for inactive and active tally
  !!
  subroutine init(self, dict)
    class(vizPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)        :: dict
    class(dictionary),pointer               :: tempDict
    type(dictionary)                        :: locDict1, locDict2
    integer(shortInt)                       :: seed_temp
    integer(longInt)                        :: seed
    character(10)                           :: time
    character(8)                            :: date
    character(:),allocatable                :: string
    character(nameLen)                      :: nucData
    class(nuclearData),pointer              :: nucData_ptr
    type(visualiser)                        :: viz
    class(geometry), pointer                :: geom
    character(100), parameter :: Here ='init (vizPhysicsPackage_class.f90)'

    print *,"VISUALISATION ONLY"

    ! Read calculation settings
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

    ! Call visualisation
    if (dict % isPresent('viz')) then
      tempDict => dict % getDictPtr('viz')
      geom => self % geom
      call viz % init(geom, tempDict)
      print *,"Visualisation complete"
    else 
      print *, "No visualisation specified" 
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
