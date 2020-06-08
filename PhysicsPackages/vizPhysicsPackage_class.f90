module vizPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError 
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use cellGeometry_inter,             only : cellGeometry

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_getMatNames => getMatNames
  use nuclearDatabase_inter,          only : nuclearDatabase

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
    class(cellGeometry), pointer :: geom => null()
    type(visualiser)             :: viz

    ! Settings
    character(pathLen) :: outputFile

    ! Timer bins
    integer(shortInt) :: timerMain

  contains
    procedure :: init
    procedure :: run
    procedure :: kill

  end type vizPhysicsPackage

contains

  !!
  !! Calls visualiser to generate visualisation
  !!
  subroutine run(self)
    class(vizPhysicsPackage), intent(inout) :: self

    print *, "Constructing visualisation"
    call self % viz % makeViz()
    call self % viz % kill()

  end subroutine

  !!
  !! Initialise from individual components and dictionaries
  !!
  subroutine init(self, dict)
    class(vizPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)        :: dict
    class(dictionary),pointer               :: tempDict
    class(geometry), pointer                :: geom
    character(100), parameter :: Here ='init (vizPhysicsPackage_class.f90)'

    ! Read outputfile path
    call dict % getOrDefault(self % outputFile,'outputFile','./output')

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    self % geom => new_cellGeometry_ptr(tempDict, ndReg_getMatNames())

    ! Call visualisation
    if (dict % isPresent('viz')) then
      print *, "Initialising visualiser"
      tempDict => dict % getDictPtr('viz')
      geom => self % geom
      call self % viz % init(geom, tempDict)
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
