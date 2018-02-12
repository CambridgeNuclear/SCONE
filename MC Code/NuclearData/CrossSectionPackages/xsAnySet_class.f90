module xsAnySet_class

  use numPrecision
  use genericProcedures, only : fatalError, linFind, searchError

  implicit none
  private

  !! **** Unfinished ****!!
  !! **** Need to think more about the interface
  !! ***** UNFINISHED ***** OBSOLETE ****** DO NOT USE
  !!
  !! Set of reaction cross-sections and their coresponding MT numbers.
  !!
  type, public :: xsAnySet
    private
    real(defReal),dimension(:),allocatable    :: xs
    integer(defReal),dimension(:),pointer     :: MTmask
  contains
    procedure :: init
    procedure :: xsOf
   ! procedure :: subset

  end type xsAnySet

contains

  subroutine init(self,xs,MTmask)
    class(xsAnySet), intent(inout)                   :: self
    real(defReal),dimension(:), intent(in)           :: xs
    integer(defReal),dimension(:),pointer,intent(in) :: MTmask
    character(100),parameter                         :: Here = 'init (xsAnySet_class.f90)'

    ! Invalid input checks
    if(size(xs) /= size(MTmask)) call fatalError(Here,'xs and MTmask have diffrent size.')
    if(any(xs < 0.0 ))           call fatalError(Here,'-ve cross-sections are present!')
    if(.not.associated(MTmask))  call fatalError(Here,'MTmask is not associated!')

    ! Load data
    self % MTmask => MTmask
    self % xs     =  xs

  end subroutine init


  function xsOf(self,MT) result(xs)
    class(xsAnySet), intent(in)   :: self
    integer(shortInt), intent(in) :: MT
    real(defReal)                 :: xs
    integer(shortInt)             :: idx
    character(100),parameter      :: Here = 'xsOf (xsAnySet_class.f90)'

    idx = linFind(self % MT, MT)
    call searchError(idx,Here)

    xs = self % xs(idx)

  end function xsOf



end module xsAnySet_class
