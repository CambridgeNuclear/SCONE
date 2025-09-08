module tallyFilterSlot_class

  use numPrecision
  use dictionary_class,        only : dictionary
  use particle_class,          only : particleState
  use tallyFilter_inter,       only : tallyFilter
  use tallyFilterFactory_func, only : new_tallyFilter

  implicit none
  private

  !!
  !! Container for polymorphic instances of tallyFilter
  !! It is itself a tallyFilter
  !! Init functions uses tallyFilterFactory to build any type of tallyFilter as specified in
  !! the provided dictionary
  !!
  type, public,extends(tallyFilter) :: tallyFilterSlot
    private
    class(tallyFilter), allocatable :: slot
  contains
    ! Implementation of Superclass interface
    procedure :: init
    procedure :: isPass

    ! Class specific procedures
    procedure :: moveAllocFrom
    procedure :: kill

  end type tallyFilterSlot

contains

  !!
  !! Initialise tallyFilterSlot
  !! Builds an instance of a tallyFilter in tallyFilterFactory and stores it in a slot
  !!
  subroutine init(self, dict)
    class(tallyFilterSlot), intent(inout) :: self
    class(dictionary), intent(in)         :: dict

    call self % kill()
    call new_tallyFilter(self % slot, dict)

  end subroutine init

  !!
  !! Call filter inside the slot
  !!
  elemental function isPass(self, state) result(passed)
    class(tallyFilterSlot), intent(in) :: self
    class(particleState), intent(in)   :: state
    logical(defBool)                   :: passed

    passed = self % slot % isPass(state)

  end function isPass

  !!
  !! Move allocation from allocatable RHS into slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(tallyFilterSlot), intent(inout)          :: LHS
    class(tallyFilter), allocatable, intent(inout) :: RHS

    call LHS % kill()
    call move_alloc(RHS, LHS % slot)

  end subroutine moveAllocFrom

  !!
  !! Deallocate content
  !!
  subroutine kill(self)
    class(tallyFilterSlot), intent(inout) :: self

    if(allocated(self % slot)) deallocate(self % slot)

  end subroutine kill


end module tallyFilterSlot_class
