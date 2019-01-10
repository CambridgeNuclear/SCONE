module tallyMap1D_inter

  use numPrecision
  use particle_class,   only : particleState
  use dictionary_class, only : dictionary
  use outputFile_class, only : outputFile
  use tallyMap_inter,   only : tallyMap

  implicit none
  private

  !!
  !! Abstract interface to separate simple 1D maps from more complex types
  !!
  type, public,extends(tallyMap),abstract :: tallyMap1D
    private

  contains
    procedure(init),deferred        :: init        ! Initialise tallyMap from dictionary
    procedure(bins),deferred        :: bins        ! Return number of bins along dimension D
    procedure                       :: dimensions  ! Return number of dimensions
    procedure(map),deferred         :: map         ! Map particle to a bin
    procedure(getAxisName),deferred :: getAxisName ! Return character describing variable of devision
    procedure(print),deferred       :: print       ! Print values associated with bins to outputfile
  end type tallyMap1D

  abstract interface

    !!
    !! Initialise tallyMap1D from a dictionary
    !!
    subroutine init(self, dict)
      import :: tallyMap1D, &
                dictionary
      class(tallyMap1D), intent(inout) :: self
      class(dictionary), intent(in)  :: dict
    end subroutine init

    !!
    !! Return total number of bins in this division along Dimension D
    !!
    elemental function bins(self, D) result(N)
      import :: tallyMap1D, &
                shortInt
      class(tallyMap1D), intent(in)   :: self
      integer(shortInt), intent(in) :: D
      integer(shortInt)             :: N
    end function bins

    !!
    !! Map particle to a single bin. Return 0 for particle out of division
    !!
    elemental function map(self,state) result(idx)
      import :: tallyMap1D,      &
                particleState, &
                shortInt
      class(tallyMap1D), intent(in)      :: self
      class(particleState), intent(in) :: state
      integer(shortInt)                :: idx
    end function map

    !!
    !! Return string that describes variable used to divide event space
    !!
    function getAxisName(self) result(name)
      import :: tallyMap1D, &
                nameLen
      class(tallyMap1D), intent(in) :: self
      character(nameLen)          :: name
    end function getAxisName

    !!
    !! Add information about division axis to the output file
    !!
    subroutine print(self,out)
      import :: tallyMap1D, &
                outputFile
      class(tallyMap1D), intent(in)      :: self
      class(outpuTFile), intent(inout) :: out
    end subroutine print

  end interface

contains

  !!
  !! Return number of dimensions for 1D map (1)
  !!
  elemental function dimensions(self) result(D)
    class(tallyMap1D), intent(in) :: self
    integer(shortInt)             :: D

    D = 1

  end function dimensions
    
end module tallyMap1D_inter
