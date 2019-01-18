module tallyResponseSlot_class

  use numPrecision
  use dictionary_class,          only : dictionary
  use particle_class,            only : particle
  use tallyResponse_inter,       only : tallyResponse
  use tallyResponseFactory_func, only : new_tallyResponse

  implicit none
  private

  !!
  !! Container for polymorphic instances of tallyResponse
  !! It is itself a tallyResponse
  !! Init functions uses tallyResponseFactory to build any type of tallyResponse as specified in
  !! the provided dictionary
  !! If get is called on uninitialised slot result is undefined (probably SEG Error)
  !!
  type, public :: tallyResponseSlot
    private
    class(tallyResponse), allocatable :: slot
  contains
    ! Abstract interface implementation
    procedure :: init
    procedure :: get

    ! Class specific procedures
    procedure :: moveAllocFrom
    procedure :: kill

  end type tallyResponseSlot

contains

  !!
  !! Initialise tallyResponseSlot
  !! Builds an instance of a tallyResponse in tallyResponseFactory and stores it in a slot
  !!
  subroutine init(self, dict)
    class(tallyResponseSlot), intent(inout) :: self
    class(dictionary), intent(in)           :: dict

    call self % kill()
    call new_tallyResponse(self % slot, dict)

  end subroutine init

  !!
  !! Get value of the response from the content of the slot
  !! If slot is unallocated (uninitialised) result is undefined (probably SEG ERROR)
  !!
  elemental function get(self, p) result(value)
    class(tallyResponseSlot), intent(in) :: self
    class(particle), intent(in)          :: p
    real(defReal)                        :: value

    value = self % slot % get(p)

  end function get

  !!
  !! Move allocation from allocatable RHS into slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(tallyResponseSlot), intent(inout)          :: LHS
    class(tallyResponse), allocatable, intent(inout) :: RHS

    call LHS % kill()
    call move_alloc(RHS, LHS % slot)

  end subroutine moveAllocFrom

  !!
  !! Deallocate content
  !!
  subroutine kill(self)
    class(tallyResponseSlot), intent(inout) :: self

    if(allocated(self % slot)) deallocate(self % slot)

  end subroutine kill


end module tallyResponseSlot_class
