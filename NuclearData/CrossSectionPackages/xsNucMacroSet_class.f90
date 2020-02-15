module xsNucMacroSet_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  !!
  !! Stores and transfers macroscopic nuclide total XSs in a material.
  !!
  type, public :: xsNucMacroSet
    integer(shortInt), dimension(:), allocatable, private :: nucIndices
    real(defReal), dimension(:), allocatable              :: nucTotalXS
    real(defReal)                                         :: matTotalXS
  contains
    procedure :: invert
    procedure :: init
  end type xsNucMacroSet

  !!
  !! Pointer wrapper to allow safe access through the interface
  !!
  type, public :: xsNucMacroSet_ptr
    private
    type(xsNucMacroSet), pointer :: ptr => null()
  contains
    generic :: assignment(=) => shallowCopy, shallowCopy_pointer
    procedure :: shallowCopy
    procedure :: shallowCopy_pointer
    procedure :: invert => invert_ptr
  end type xsNucMacroSet_ptr

contains

  !!
  !! Allocate space and load nuclide indices
  !!
  subroutine init(self,nucIndices)
    class(xsNucMacroSet), intent(inout)            :: self
    integer(shortInt),dimension(:),intent(in)  :: nucIndices

    if(allocated(self % nucIndices)) deallocate (self % nucIndices)
    if(allocated(self % nucTotalXS)) deallocate (self % nucTotalXS)

    self % nucIndices = nucIndices
    allocate(self % nucTotalXS (size(nucIndices)))

  end subroutine init

  !!
  !! Invert probability distribution
  !!
  function invert(self,r) result(nucIdx)
    class(xsNucMacroSet), intent(in) :: self
    real(defReal), intent(in)        :: r
    integer(shortInt)                :: nucIdx
    real(defReal)                    :: r_scaled
    character(100),parameter         :: Here = 'invert (xsNucMacroSet_class.f90)'
    integer(shortInt)                :: i

    ! Scale r to total reaction cross-section
    r_scaled = r * self % matTotalXS

    ! Search nucTotalXSs (scaled PDF) with linear search
    do i = 1,size(self % nucTotalXS)

      ! Decrement scaled random number
      r_scaled = r_scaled - self % nucTotalXS(i)

      if( r_scaled <= 0) then
        nucIdx = self % nucIndices(i)
        return

      end if
    end do

    call fatalError(Here,'Provided number to invert cdf must be < 1 ')

    ! Avoid compiler warning (Serpent territory?...)
    NucIdx = 0

  end function invert

!!**************************************************************************************************
!! Pointer Wrapper Procedures
!!
!!**************************************************************************************************

  !!
  !! Assignment Pointer Wrapper to Pointer Wrapper
  !!
  subroutine shallowCopy(LHS,RHS)
    class(xsNucMacroSet_ptr), intent(inout) :: LHS
    type(xsNucMacroSet_ptr), intent(in)     :: RHS

    LHS % ptr => RHS % ptr
    
  end subroutine shallowCopy

  !!
  !! Assignment Pointer to Pointer Wrapper
  !!
  subroutine shallowCopy_pointer(LHS,RHS)
    class(xsNucMacroSet_ptr), intent(inout) :: LHS
    type(xsNucMacroSet),pointer, intent(in) :: RHS

    ! (Will kill for templates at this point... MAK)
    LHS % ptr => RHS

  end subroutine shallowCopy_pointer

  !!
  !! Access invert procedure
  !!
  function invert_ptr(self,r) result(nucIdx)
    class(xsNucMacroSet_ptr), intent(in) :: self
    real(defReal), intent(in)            :: r
    integer(shortInt)                    :: nucIdx

    nucIdx = self % ptr % invert(r)

  end function invert_ptr
end module xsNucMacroSet_class
