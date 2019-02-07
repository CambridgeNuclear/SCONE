module angleLawENDFslot_class

  use numPrecision
  use angleLawENDF_inter,       only : angleLawENDF
  use RNG_class,                only : RNG
  use aceCard_class,            only : aceCard
  use angleLawENDFFactory_func, only : new_angleLawENDF


  implicit none
  private

  type, public,extends(angleLawENDF) :: angleLawENDFslot
    private
    class(angleLawENDF),allocatable :: slot
  contains
    ! Duplicate interfacte of the superclass
    procedure :: init
    procedure :: sample
    procedure :: probabilityOf

    ! Define assignment
    procedure :: moveAllocFrom

  end type angleLawENDFslot

contains

  !!
  !! Initialise from aceCard and MT number
  !! Use factory to allocate contents of the slot
  !!
  subroutine init(self, ACE, MT)
    class(angleLawENDFslot), intent(inout):: self
    class(aceCard), intent(inout)         :: ACE
    integer(shortInt), intent(in)         :: MT

    call new_angleLawENDF(self % slot, ACE, MT)

  end subroutine init

  !!
  !! Given collison energy and random number generator sample mu
  !!
  function sample(self,E,rand) result (mu)
    class(angleLawENDFslot), intent(in) :: self
    real(defReal), intent(in)           :: E
    class(RNG), intent(inout)           :: rand
    real(defReal)                       :: mu

    mu = self % slot % sample(E,rand)

  end function sample

  !!
  !! Return probability density of mu at collision energy E
  !!
  function probabilityOf(self,mu,E) result (prob)
    class(angleLawENDFslot), intent(in) :: self
    real(defReal), intent(in)           :: E, mu
    real(defReal)                       :: prob

    prob = self % slot % probabilityOf(mu,E)

  end function probabilityOf

  !!
  !! Move allocation from RHS to LHS slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(angleLawENDFslot), intent(inout) :: LHS
    type(angleLawENDFslot), intent(inout)  :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    call move_alloc(RHS % slot, LHS % slot)

  end subroutine moveAllocFrom

end module angleLawENDFslot_class
