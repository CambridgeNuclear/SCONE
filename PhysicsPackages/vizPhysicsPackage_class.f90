module vizPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError
  use dictionary_class,               only : dictionary

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx
  use geometryFactory_func,           only : new_geometry

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_getMatNames => getMatNames
  use nuclearDatabase_inter,          only : nuclearDatabase

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
    class(geometry), pointer :: geom => null()
    integer(shortInt)        :: geomIdx = 0
    type(visualiser)         :: viz

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
    character(nameLen)                      :: geomName
    character(100), parameter :: Here ='init (vizPhysicsPackage_class.f90)'

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'visualGeom'
    call new_geometry(tempDict, geomName)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

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
