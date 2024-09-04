module collisionOperator_class

  use numPrecision
  use genericProcedures,     only : fatalError, numToChar
  use dictionary_class,      only : dictionary
  use particle_class,        only : particle, P_NEUTRON, P_PHOTON, printType
  use particleDungeon_class, only : particleDungeon

  ! Tally interfaces
  use tallyAdmin_class,      only : tallyAdmin

  ! Collision Processors
  use collisionProcessor_inter,       only : collisionProcessor
  use collisionProcessorFactory_func, only : new_collisionProcessor

  implicit none
  private

  !!
  !! Local parameters
  !!
  integer(shortInt), parameter :: MAX_P_ID = max(P_NEUTRON, P_PHOTON)
  integer(shortInt), parameter :: P_MG = 1, P_CE = 2
  integer(shortInt), parameter :: UNDEF_PHYSICS = -1
  integer(shortInt), parameter :: MAX_PHYSICS = 3


  !!
  !! Local helper type to store polymorphic collisionProcessors in an array
  !!
  type, private :: collProc
    class(collisionProcessor), allocatable :: proc
  end type collProc


  !!
  !! Scalar collision operator
  !!  -> Maps particles of diffrent types to approperiate implementation
  !!     of collision physics (collisionProcessor).
  !!  -> Uses lookup table for to quickly map diffrent combination of physical(neutron, photon ...)
  !!     and processing type (CE, MG) to approperiate physics.
  !!  -> Can store up to 3 diffrent physics types
  !!  -> Gives fatal error if particle type is not recognised or not-supported
  !!
  !! Sample dictionary input( provisional will change):
  !!  collOpName {
  !!    #neutronCE {<collisonProcessorDefinition>} #
  !!    #neutronMG {<collisonProcessorDefinition>} #
  !!  }
  !!
  type, public :: collisionOperator
    private
    type(collProc), dimension(MAX_PHYSICS)    :: physicsTable
    integer(shortInt), dimension(2, MAX_P_ID) :: lookupTable = UNDEF_PHYSICS

  contains
    ! Build procedures
    procedure :: init
    procedure :: kill

    ! Use procedures
    procedure :: collide


  end type collisionOperator

contains

  !!
  !! Initialise collision operator
  !!
  subroutine init(self, dict)
    class(collisionOperator), intent(inout) :: self
    class(dictionary), intent(in)           :: dict

    if(dict % isPresent('neutronCE')) then
      call new_collisionProcessor(self % physicsTable(1) % proc, dict % getDictPtr('neutronCE'))
      self % lookupTable(P_CE, P_NEUTRON) = 1
    end if

    if(dict % isPresent('neutronMG')) then
      call new_collisionProcessor(self % physicsTable(2) % proc, dict % getDictPtr('neutronMG'))
      self % lookupTable(P_MG, P_NEUTRON) = 2
    end if

  end subroutine init

  !!
  !! Clear collision operator. Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(collisionOperator), intent(inout) :: self
    integer(shortInt)                       :: i

    ! Set default to lookup table
    self % lookupTable = UNDEF_PHYSICS

    ! Deallocate collision processors
    do i=1,size(self % physicsTable)
      if(allocated(self % physicsTable(i) % proc)) then
        deallocate( self % physicsTable(i) % proc)
      end if
    end do

  end subroutine kill

  !!
  !! Determine type of the particle and call approperiate collisionProcessor
  !!
  subroutine collide(self, p, tally, thisCycle, nextCycle)
    class(collisionOperator), intent(inout) :: self
    class(particle), intent(inout)           :: p
    type(tallyAdmin), intent(inout)          :: tally
    class(particleDungeon),intent(inout)     :: thisCycle
    class(particleDungeon),intent(inout)     :: nextCycle
    integer(shortInt)                        :: idx, procType
    character(100), parameter :: Here = 'collide (collisionOperator_class.f90)'

    ! Select processing index with ternary expression
    if(p % isMG) then
      procType = P_MG
    else
      procType = P_CE
    end if

    ! Varify that type is valid
    if( p % type <= 0 .or. p % type > MAX_P_ID) then
      call fatalError(Here, 'Type of the particle is invalid: ' // numToChar(p % type))
    end if

    ! Get index
    idx = self % lookupTable(procType, p % type)

    ! Verify index
    if(idx == UNDEF_PHYSICS) then
      call fatalError(Here,'Physics is not defined for particle of type : '// p % typeToChar())
    end if

    ! Call physics
    call self % physicsTable(idx) % proc % collide(p, tally, thisCycle, nextCycle)

  end subroutine collide

end module collisionOperator_class
