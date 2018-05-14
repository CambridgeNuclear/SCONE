module tallyResponseSlot_class

  use numPrecision
  use genericProcedures,    only : fatalError
  use particle_class,       only : particle, phaseCoord
  use tallyResponse_inter,  only : tallyResponse

  implicit none
  private

  !!
  !! Type to store a polymorphic tally response
  !!
  type, public,extends(tallyResponse) :: tallyResponseSlot
    private
    class(tallyResponse),allocatable :: response
  contains
    ! Interface procedures
    procedure :: getCollisionScore

    ! Build procedures
    generic           :: assignment(=) => assignResponse
    procedure,private :: assignResponse

  end type tallyResponseSlot

contains

  !!
  !! Given collision parameters return score function (f(x)) value for collision
  !! f needs to be multiplied by flux estimate before scoring
  !!
  function getCollisionScore(self,pre,post,MT,muL) result (f)
    class(tallyResponseSlot), intent(inout) :: self
    class(phaseCoord), intent(in)           :: pre
    class(particle), intent(in)             :: post
    integer(shortInt), intent(in)           :: MT
    real(defReal), intent(in)               :: muL
    real(defReal)                           :: f

    f = self % response % getCollisionScore(pre, post, MT, muL)

  end function getCollisionScore

  !!
  !! Load response into the slot
  !! Current content will be deallocated if present
  !! "response" argument variable will be deallocated in the subroutine
  !! (Allocation is moved from argument variable to the slot)
  !!
  subroutine assignResponse(LHS,RHS)
    class(tallyResponseSlot), intent(inout)        :: LHS
    class(tallyResponse),allocatable,intent(in)    :: RHS
    character(100),parameter        :: Here='assignResponse (tallyResponseSlot_class.f90)'

    ! Deallocate currently stored response if present
    if( allocated( LHS % response)) deallocate( LHS % response)

    if(.not.allocated(RHS)) then
      call fatalError(Here,'Trying to load unallocated tally response')
    end if


    ! Move allocation into slot
    !call move_alloc(response, self % response)

    allocate(LHS % response, source = RHS)

  end subroutine assignResponse

    
end module tallyResponseSlot_class
