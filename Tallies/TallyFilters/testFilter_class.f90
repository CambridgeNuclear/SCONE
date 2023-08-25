module testFilter_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use particle_class,    only : particleState
  use dictionary_class,  only : dictionary
  use tallyFilter_inter, only : tallyFilter

  implicit none
  private

  !!
  !! Very simple filter used for testing of other components only
  !!   Returns true if state % matIdx (mi) :  minIdx <= mi <= maxIdx
  !!
  type, public,extends(tallyFilter) :: testFilter
    private
    integer(shortInt) :: minIdx = 0
    integer(shortInt) :: maxIdx = 0
  contains
    procedure :: init
    procedure :: isPass

  end type testFilter

contains

  !!
  !! Initialise testFilter from dictionary
  !!
  subroutine init(self,dict)
    class(testFilter), intent(inout) :: self
    class(dictionary), intent(in)    :: dict
    character(100), parameter :: Here ='init (testFilter_class.f90)'

    call dict % get(self % minIdx,'minIdx')
    call dict % get(self % maxIdx,'maxIdx')

    ! Verify
    if( self % minIdx > self % maxIdx) call fatalError(Here,'minIdx > maxIdx')

  end subroutine init

  !!
  !! Returns true if energy value is between specified bounds
  !!
  elemental function isPass(self,state) result(passed)
    class(testFilter), intent(in)    :: self
    class(particleState), intent(in) :: state
    logical(defBool)                 :: passed

    passed = (self % minIdx <= state % matIdx) .and. (state % matIdx <= self % maxIdx)

  end function isPass

end module testFilter_class
