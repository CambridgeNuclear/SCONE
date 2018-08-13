module tallyMapSlot_class

  use numPrecision
  use particle_class,   only : particle
  use outputFile_class, only : outputFile
  use tallyMap_inter,   only : tallyMap

  implicit none
  private



  type, public,extends(tallyMap) :: tallyMapSlot
    private
    class(tallyMap),allocatable :: slot
  contains
    ! Superclass interface implementaction
    procedure  :: bins        ! Return number of bins
    procedure  :: map         ! Map particle to a bin
    procedure  :: getAxisName ! Return character describing variable of devision
    procedure  :: print       ! Print values associated with bins to outputfile

    ! Class specific procedures
    generic   :: assignment(=) => copy
    procedure :: copy
    procedure :: moveAllocFrom
  end type tallyMapSlot

contains

  !!
  !! Return total number of bins in this division
  !!
  pure function bins(self) result(N)
    class(tallyMapSlot), intent(in) :: self
    integer(shortInt)               :: N

    N = self % slot % bins()

  end function bins

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  elemental function map(self,p) result(idx)
    class(tallyMapSlot), intent(in) :: self
    class(particle), intent(in)     :: p
    integer(shortInt)               :: idx

    idx = self % slot % map(p)

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  pure function getAxisName(self) result(name)
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


end module tallyMapSlot_class
