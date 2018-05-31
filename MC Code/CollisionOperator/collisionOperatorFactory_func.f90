module collisionOperatorFactory_func

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  ! Nuclear Data an geometry interface
  use nuclearData_inter, only : nuclearData

  ! Collision Operators
  use collisionOperatorBase_inter,    only : collisionOperatorBase
  use perMaterialCollisionOpMG_class, only : perMaterialCollisionOpMG
  use perNuclideCollisionOpCE_class,  only : perNuclideCollisionOpCE

  ! Tally Interfaces
  use tallyAdminBase_class,           only : tallyAdminBase

  implicit none

  ! *** ADD NAME OF A NEW TRANSPORT OPERATOR HERE ***!
  ! List that contains all accaptable types of transport operators
  ! It is printed if type was unrecognised
  ! NOTE:
  ! For now  it is necessary to adjust trailing blanks so all enteries have the same length
  character(nameLen),dimension(*),parameter :: AVALIBLE_collisionOps = [ 'perMaterialCollisionOpMG', &
                                                                         'perNuclideCollisionOpCE ']

contains

  function new_collisionOperator_ptr(nucData,dict) result(new)
    class(nuclearData),pointer,intent(in)    :: nucData
    class(dictionary),intent(in)             :: dict
    class(collisionOperatorBase),pointer     :: new
    character(nameLen)                       :: type
    character(100),parameter :: Here = 'new_collisionOperator_ptr (collisionOperatorFactory_func.f90)'

    ! Obtain string that specifies type to be built
    call dict % get(type,'type')

    ! Allocate approperiate subclass of nuclear data
    ! *** ADD CASE STATEMENT FOR A NEW COLLISION OPERATOR BELOW ***!
    select case(type)
      case('perMaterialCollisionOpMG')
        allocate(perMaterialCollisionOpMG :: new)
        call new % init(nucData, dict)

      case('perNuclideCollisionOpCE')
        allocate(perNuclideCollisionOpCE :: new)
        call new % init(nucData, dict)

      case default
        print *, AVALIBLE_collisionOps
        call fatalError(Here, 'Unrecognised type of collisionOperator: ' // trim(type))

    end select

  end function new_collisionOperator_ptr

end module collisionOperatorFactory_func
