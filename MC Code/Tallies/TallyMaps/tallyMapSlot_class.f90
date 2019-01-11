module tallyMapSlot_class

  use numPrecision
  use particle_class,       only : particleState
  use outputFile_class,     only : outputFile
  use dictionary_class,     only : dictionary
  use tallyMap_inter,       only : tallyMap
  use tallyMapFactory_func, only : new_tallyMap

  implicit none
  private

  !!
  !!
  !!
  type, public,extends(tallyMap) :: tallyMapSlot
    private
    class(tallyMap),allocatable :: slot
  contains
    ! Superclass interface implementaction
    procedure  :: init        ! Initialise content from dictionary
    procedure  :: dimensions  ! Return number of dimensions
    procedure  :: bins        ! Return number of bins
    procedure  :: map         ! Map particle to a bin
    procedure  :: getAxisName ! Return character describing variable of devision
    procedure  :: print       ! Print values associated with bins to outputfile

    ! Class specific procedures
    generic   :: assignment(=) => copy
    procedure :: copy
    procedure :: moveAllocFrom
    procedure :: kill

  end type tallyMapSlot

contains

  !!
  !! Shortcut to factory.
  !! Builds object in a factory from a dictionary and stores it in a slot
  !!
  subroutine init(self, dict)
    class(tallyMapSlot), intent(inout) :: self
    class(dictionary), intent(in)      :: dict

    call self % kill()
    call new_tallyMap(self % slot, dict)

  end subroutine init

  !!
  !! Return total number of bins in this division
  !!
  elemental function bins(self, D) result(N)
    class(tallyMapSlot), intent(in) :: self
    integer(shortInt),intent(in)    :: D
    integer(shortInt)               :: N

    N = self % slot % bins(D)

  end function bins

  !!
  !! Return number of dimensions
  !!
  elemental function dimensions(self) result(D)
    class(tallyMapSlot), intent(in)    :: self
    integer(shortInt)                  :: D

    D = self % slot % dimensions()

  end function dimensions


  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  elemental function map(self,state) result(idx)
    class(tallyMapSlot), intent(in)  :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx

    idx = self % slot % map(state)

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  function getAxisName(self) result(name)
    class(tallyMapSlot), intent(in) :: self
    character(nameLen)              :: name

    name = self % slot % getAxisName()

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  subroutine print(self,out)
    class(tallyMapSlot), intent(in)  :: self
    class(outpuTFile), intent(inout) :: out

    call self % slot % print(out)

  end subroutine print

  !!
  !! Copy RHS into slot of LHS
  !! Be carefull about loading slots into slots
  !! It will work by function call chain may hurt performance
  !!
  subroutine copy(LHS,RHS)
    class(tallyMapSlot), intent(inout) :: LHS
    class(tallyMap),intent(in)         :: RHS

    if(allocated(LHS % slot)) deallocate(LHS % slot)

    allocate(LHS % slot, source = RHS )

  end subroutine copy

  !!
  !! Move allocation from allocatable RHS into slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(tallyMapSlot), intent(inout)          :: LHS
    class(tallyMap), allocatable, intent(inout) :: RHS

    if(allocated(LHS % slot)) deallocate(LHS % slot)
    call move_alloc(RHS, LHS % slot)

  end subroutine moveAllocFrom

  !!
  !! Deallocate content of the slot
  !!
  subroutine kill(self)
    class(tallyMapSlot), intent(inout) :: self

    if(allocated(self % slot)) deallocate(self % slot)

  end subroutine kill


end module tallyMapSlot_class
