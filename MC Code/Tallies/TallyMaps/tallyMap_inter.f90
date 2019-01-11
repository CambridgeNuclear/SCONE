module tallyMap_inter

  use numPrecision
  use particle_class,   only : particleState
  use dictionary_class, only : dictionary
  use outputFile_class, only : outputFile

  implicit none
  private

  !!
  !! Abstract interface for maping event to approperiate bins
  !!   Given particle it returns an index of the bin
  !!   Can return index 0, which means that event should not be scored
  !!
  !!
  type, public,abstract :: tallyMap
    private

  contains
    procedure(init),deferred        :: init        ! Initialise tallyMap from dictionary
    procedure(bins),deferred        :: bins        ! Return number of bins along dimension D
    procedure(dimensions),deferred  :: dimensions  ! Return number of dimensions
    procedure(map),deferred         :: map         ! Map particle to a bin
    procedure(getAxisName),deferred :: getAxisName ! Return character describing variable of devision
    procedure(print),deferred       :: print       ! Print values associated with bins to outputfile
  end type tallyMap

  abstract interface

    !!
    !! Initialise tallyMap from a dictionary
    !!
    subroutine init(self, dict)
      import :: tallyMap, &
                dictionary
      class(tallyMap), intent(inout) :: self
      class(dictionary), intent(in)  :: dict
    end subroutine init

    !!
    !! Return total number of bins in this division along Dimension D
    !! If D==0 return sum of bin numbers for all dimensions
    !!
    elemental function bins(self, D) result(N)
      import :: tallyMap, &
                shortInt
      class(tallyMap), intent(in)   :: self
      integer(shortInt), intent(in) :: D
      integer(shortInt)             :: N
    end function bins

    !!
    !! Return number of dimensions
    !!
    elemental function dimensions(self) result(D)
      import :: tallyMap, &
                shortInt
      class(tallyMap), intent(in)    :: self
      integer(shortInt)              :: D
    end function dimensions

    !!
    !! Map particle to a single bin. Return 0 for particle out of division
    !!
    elemental function map(self,state) result(idx)
      import :: tallyMap,      &
                particleState, &
                shortInt
      class(tallyMap), intent(in)      :: self
      class(particleState), intent(in) :: state
      integer(shortInt)                :: idx
    end function map

    !!
    !! Return string that describes variable used to divide event space
    !!
    function getAxisName(self) result(name)
      import :: tallyMap, &
                nameLen
      class(tallyMap), intent(in) :: self
      character(nameLen)          :: name
    end function getAxisName

    !!
    !! Add information about division axis to the output file
    !!
    subroutine print(self,out)
      import :: tallyMap, &
                outputFile
      class(tallyMap), intent(in)      :: self
      class(outputFile), intent(inout) :: out
    end subroutine print

  end interface
end module tallyMap_inter
