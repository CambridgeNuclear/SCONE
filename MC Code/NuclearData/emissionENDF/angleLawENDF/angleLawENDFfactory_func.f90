module angleLawENDFfactory_func

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError
  use aceCard_class,     only : aceCard

  ! Diffrent angle laws
  use angleLawENDF_inter,     only : angleLawENDF
  use isotropicAngle_class,   only : isotropicAngle
  use noAngle_class,          only : noAngle
  use tabularAngle_class,     only : tabularAngle

  implicit none
  private

  ! Define public interface
  public :: new_angleLawENDF
  public :: new_angleLawENDF_ptr

contains

  !!
  !! Returns an allocatable angleLawENDF from aceCard and MT number
  !! aceCard can be in any poistion. Its position changes at output.
  !!
  function new_angleLawENDF(ACE,MT) result(new)
    type(aceCard), intent(inout)                     :: ACE
    integer(shortInt), intent(in)                    :: MT
    class(angleLawENDF),allocatable                  :: new
    integer(shortInt)                                :: LOCB
    character(100),parameter :: Here='new_angleLawENDF (angleLawENDFfactory_func.f90)'

    ! Check if reaction is absorbtion
    if( MT >= N_disap) then ! Reaction is absorbtion MT > 100
      !call fatalError(Here,'Cannot build angle law for absorbtion reaction')
      allocate(new, source = noAngle() )
      return

    end if

    ! Obtain LOCB
    LOCB = ACE % LOCBforMT(MT)

    ! Return error if LOCB corresponds to LAW 44
    if(LOCB == LOCB_CORRELATED) then
      call fatalError(Here,'Cannot build angle law for correlated angle-energy')

    end if

    ! Return if LOCB corresponds to LOCB_ISOTROPIC
    if(LOCB == LOCB_ISOTROPIC) then
      allocate(new, source = isotropicAngle() )
      return
    end if

    ! Set read head of ACE to beggining of angle data
    call ACE % setToAngleMT(MT)

    allocate(new, source = tabularAngle(ACE))


  end function new_angleLawENDF

  !!
  !! Returns a pointer to allocated angle LawENDF from aceCard and MT number
  !!
  function new_angleLawENDF_ptr(ACE,MT) result(new)
    type(aceCard), intent(inout)    :: ACE
    integer(shortInt), intent(in)   :: MT
    class(angleLawENDF),pointer     :: new

    ! Allocate pointer and copy data from local allocatable
    allocate( new, source = new_angleLawENDF(ACE,MT) )

  end function new_angleLawENDF_ptr

end module angleLawENDFfactory_func
