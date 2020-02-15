module tallyResponseSlot_class

  use numPrecision
  use dictionary_class,          only : dictionary
  use particle_class,            only : particle
  use tallyResponse_inter,       only : tallyResponse
  use tallyResponseFactory_func, only : new_tallyResponse

  ! Nuclear Data interface
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private

  !!
  !! Container for polymorphic instances of tallyResponse
  !!
  !! It is itself a tallyResponse
  !! Init functions uses tallyResponseFactory to build any type of tallyResponse as specified in
  !! the provided dictionary
  !! If get is called on uninitialised slot result is undefined (probably SEG Error)
  !!
  !! Private Members:
  !!   slot -> Allocatable tallyResponse content
  !!
  !! Interface:
  !!   tallyResponse Interface
  !!   moveAllocFrom -> Moves allocation from argument to slot
  !!
  type, public :: tallyResponseSlot
    private
    class(tallyResponse), allocatable :: slot
  contains
    ! Abstract interface implementation
    procedure :: init
    procedure :: get
    procedure :: kill

    ! Class specific procedures
    procedure :: moveAllocFrom

  end type tallyResponseSlot

contains

  !!
  !! Initialise tallyResponseSlot
  !!
  !! Builds an instance of a tallyResponse in tallyResponseFactory and stores it in a slot
  !!
  !! See tallyResponse_inter for details
  !!
  subroutine init(self, dict)
    class(tallyResponseSlot), intent(inout) :: self
    class(dictionary), intent(in)           :: dict

    call self % kill()
    call new_tallyResponse(self % slot, dict)

  end subroutine init

  !!
  !! Get value of the response from the content of the slot
  !!
  !! See tallyResponse_inter for details
  !!
  !! Errors:
  !!   If slot is unallocated (uninitialised) result is undefined (probably SEG ERROR)
  !!
  function get(self, p, xsData) result(value)
    class(tallyResponseSlot), intent(in)  :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    real(defReal)                         :: value

    value = self % slot % get(p, xsData)

  end function get

  !!
  !! Move allocation from allocatable RHS into slot
  !!
  !! RHS becomes unallocated on exit
  !!
  !! Args:
  !!   RHS [inout] -> Allocated tallyResponse to be loaded into slot
  !!
  !! Errors:
  !!   If RHS is unallocated has no effect
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(tallyResponseSlot), intent(inout)          :: LHS
    class(tallyResponse), allocatable, intent(inout) :: RHS

    call LHS % kill()
    call move_alloc(RHS, LHS % slot)

  end subroutine moveAllocFrom

  !!
  !! Return to Uninitialised State
  !!
  elemental subroutine kill(self)
    class(tallyResponseSlot), intent(inout) :: self

    if(allocated(self % slot)) deallocate(self % slot)

  end subroutine kill

end module tallyResponseSlot_class
