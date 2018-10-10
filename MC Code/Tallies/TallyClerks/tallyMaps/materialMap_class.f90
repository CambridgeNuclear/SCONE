module materialMap_class

  use numPrecision
  use genericProcedures, only : fatalError
  use particle_class,    only : particle
  use outputFile_class,  only : outputFile
  use tallyMap_inter,    only : tallyMap

  implicit none
  private

  !!
  !! Map that divides based on the material in a particle
  !!
  type, public,extends(tallyMap) :: materialMap
    private

  contains
    ! Superclass interface implementaction
    procedure  :: bins        ! Return number of bins
    procedure  :: map         ! Map particle to a bin
    procedure  :: getAxisName ! Return character describing variable of devision
    procedure  :: print       ! Print values associated with bins to outputfile

  end type materialMap

contains

  !!
  !! Return total number of bins in this division
  !!
  pure function bins(self) result(N)
    class(materialMap), intent(in)  :: self
    integer(shortInt)               :: N



  end function bins

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  elemental function map(self,p) result(idx)
    class(materialMap), intent(in)     :: self
    class(particle), intent(in)     :: p
    integer(shortInt)               :: idx


  end function map

  !!
  !! Return string that describes variable used to divide event space
  !!
  pure function getAxisName(self) result(name)
    class(materialMap), intent(in)  :: self
    character(nameLen)              :: name

    name = 'Material'

  end function getAxisName

  !!
  !! Add information about division axis to the output file
  !!
  subroutine print(self,out)
    class(materialMap), intent(in)   :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

  end subroutine print


end module materialMap_class
