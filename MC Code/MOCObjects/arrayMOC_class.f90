!!
!! Class for flux and source arrays used by MOC
!! Array is sized to match the number of FSRs and energy groups
!! Can flatten, zero, scale or increment
!!
module arrayMOC_class

  use numPrecision
  use universalVariables
  use genericProcedures

  implicit none
  private

  type, public                                 :: arrayMOC
    integer(shortInt)                          :: nFSR       !! Number of FSRs
    integer(shortInt)                          :: nGroup     !! Number of energy groups
    real(defReal), dimension(:,:), allocatable :: value      !! Array of values by group and region
  contains
    procedure          :: init
    procedure          :: flatten
    procedure          :: zero
    procedure          :: increment
    procedure          :: scale
  end type arrayMOC

contains

  !!
  !! Constructor for array
  !!
  subroutine init(self, nFSR, nGroup)
    class(arrayMOC), intent(inout) :: self
    integer(shortInt), intent(in) :: nFSR
    integer(shortInt), intent(in) :: nGroup

    allocate(self % value(nGroup,nFSR))
    self % nFSR = nFSR
    self % nGroup = nGroup
    call self % zero()

  end subroutine init

  !!
  !! Flattens profile in a particular group to a specified value
  !!
  subroutine flatten(self, val, g)
    class(arrayMOC), intent(inout) :: self
    real(defReal), intent(in) :: val
    integer(shortInt), intent(in), optional :: g

    ! If a group is provided, flatten only in that group
    if (present(g)) then
      if ((g <= self % nGroup).AND.(g > 0)) then
        self % value(g,:) = val
      else
        call fatalError('flatten, arrayMOC','Invalid group index provided')
      end if
    ! Otherwise, flatten all groups to the same value
    else
      self % value(:,:) = val
    end if
  end subroutine flatten

  !!
  !! Zeros all array values
  !!
  subroutine zero(self)
    class(arrayMOC), intent(inout) :: self
    self % value(:,:) = ZERO
  end subroutine zero

  !!
  !! Increments an array value
  !!
  subroutine increment(self, inc, g, r)
    class(arrayMOC), intent(inout) :: self
    real(defReal), intent(in) :: inc
    integer(shortInt), intent(in) :: g
    integer(shortInt), intent(in) :: r

    self % value(g,r) = self % value(g,r) + inc
  end subroutine increment

  !!
  !! Scales values in array by multiplicative factor
  !! Optionally, only does so in a chosen energy group
  !!
  subroutine scale(self, fac, g)
    class(arrayMOC), intent(inout) :: self
    real(defReal), intent(in) :: fac
    integer(shortInt), intent(in), optional :: g

    ! If a group is provided, flatten only in that group
    if (present(g)) then
      if ((g <= self % nGroup).AND.(g > 0)) then
        self % value(g,:) = val
      else
        call fatalError('scale, arrayMOC','Invalid group index provided')
      end if
    ! Otherwise, flatten all groups to the same value
    else
      self % value(:,:) = fac * self % value(:,:)
    end if
  end subroutine scale

end module arrayMOC_class
