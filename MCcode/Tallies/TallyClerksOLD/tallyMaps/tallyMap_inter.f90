module tallyMap_inter

  use numPrecision
  use particle_class,   only : particle
  use outputFile_class, only : outputFile

  implicit none
  private


  !!
  !! Abstract interface for maping event to approperiate bins
  !!   Given particle it returns an index of the bin
  !!
  !!
  type, public,abstract :: tallyMap
    private

  contains
    procedure(bins),deferred        :: bins        ! Return number of bins
    procedure(map),deferred         :: map         ! Map particle to a bin
    procedure(getAxisName),deferred :: getAxisName ! Return character describing variable of devision
    procedure(print),deferred       :: print       ! Print values associated with bins to outputfile
  end type tallyMap

  abstract interface
    !!
    !! Return total number of bins in this division
    !!
    pure function bins(self) result(N)
      import :: tallyMap, &
                shortInt
      class(tallyMap), intent(in) :: self
      integer(shortInt)           :: N
    end function bins

    !!
    !! Map particle to a single bin. Return 0 for particle out of division
    !!
    function map(self,p) result(idx)
      import :: tallyMap, &
                particle, &
                shortInt
      class(tallyMap), intent(in) :: self
      class(particle), intent(in) :: p
      integer(shortInt)           :: idx
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
      class(outpuTFile), intent(inout) :: out
    end subroutine print

  end interface
end module tallyMap_inter
