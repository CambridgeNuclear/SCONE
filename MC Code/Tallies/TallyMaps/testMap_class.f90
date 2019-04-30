module testMap_class

  use numPrecision
  use dictionary_class, only : dictionary
  use particle_class,   only : particleState
  use outputFile_class, only : outputFile
  use tallyMap1D_inter,   only : tallyMap1D, kill_super => kill

  implicit none
  private


  !!
  !! Very simple map used for testing of other components only
  !!   Given state it returns bin = matIdx or 0 if matIdx > maxIdx given in dictionary
  !!
  type, public,extends(tallyMap1D) :: testMap
    private
    integer(shortInt) :: maxIdx = 0
  contains
    procedure :: init          ! Initialise tallyMap from dictionary
    procedure :: bins          ! Return number of bins along dimension D
    procedure :: map           ! Map particle to a bin
    procedure :: getAxisName   ! Return character describing variable of devision
    procedure :: print         ! Print values associated with bins to outputfile
    procedure :: kill

  end type testMap

contains

  !!
  !! Initialise testMap from a dictionary
  !!
  subroutine init(self, dict)
    class(testMap), intent(inout)  :: self
    class(dictionary), intent(in)  :: dict

    call dict % get(self % maxIdx,'maxIdx')

  end subroutine init

  !!
  !! Return total number of bins in this division along Dimension D
  !! If D==0 return sum of bin numbers for all dimensions
  !!
  elemental function bins(self, D) result(N)
    class(testMap), intent(in)    :: self
    integer(shortInt), intent(in) :: D
    integer(shortInt)             :: N

    if (D == 1 .or. D == 0) then
      N = self % maxIdx
    else
      N = 0
    end if

  end function bins

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  elemental function map(self,state) result(idx)
    class(testMap), intent(in)       :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx

    if(state % matIdx < 0 .or. state % matIdx > self % maxIdx) then
      idx = 0
    else
      idx = state % matIdx
    end if

  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  function getAxisName(self) result(name)
    class(testMap), intent(in)  :: self
    character(nameLen)          :: name

    name='TestMap'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  subroutine print(self,out)
    class(testMap), intent(in)       :: self
    class(outputFile), intent(inout) :: out

    ! Do nothing for now

  end subroutine print

  !!
  !! Kill testMap
  !!
  elemental subroutine kill(self)
    class(testMap), intent(inout) :: self

    call kill_super(self)
    self % maxIdx = 0

  end subroutine kill
end module testMap_class
