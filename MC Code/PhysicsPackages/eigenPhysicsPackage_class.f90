module eigenPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,            only : fatalError
  use dictionary_class,             only : dictionary

  ! Particle classes and Random number generator
  use particle_class,               only : particle
  use particleDungeon_class,        only : particleDungeon
  use RNG_class,                    only : RNG

  ! Geometry & Nuclear Data
  use geometry_class,               only : geometry
  use nuclearData_inter,            only : nuclearData

  ! Operators
  use collisionOperatorBase_inter,  only : collisionOperatorBase
  use transportOperator_inter,      only : transportOperator

  ! Tallies
  use tallyInactiveAdmin_class,     only : tallyInacticeAdmin
  use tallyActiveAdmin_class,       only : tallyActiveAdmin

  implicit none
  private

  !!
  !! Physics Package for eigenvalue calculations
  !!
  type, public :: eigenPhysicsPackage
    private
    ! Building blocks
    class(nuclearData), pointer            :: nucData  => null()
    class(geometry), pointer               :: geometry => null()
    class(collisionOperatorBase), pointer  :: collOp   => null()
    class(transportOperator), pointer      :: transOp  => null()
    class(RNG), pointer                    :: pRNG     => null()
    type(tallyInactiveAdmin)               :: inactiveTally
    type(tallyActiveAdmin)                 :: activeAdmin

  contains
    procedure :: init

  end type eigenPhysicsPackage

contains

  !!
  !! Initialise from individual components and dictionaries for inactive and active tally
  !! *** This is far from final implementation. Ultimatly only single dictionery will be required
  !! *** Assumes all settings
  !!
  subroutine init(self,nucData,geometry,collOp,transOp,pRNG)
    class(eigenPhysicsPackage), intent(inout) :: self
  end subroutine init


    
end module eigenPhysicsPackage_class
