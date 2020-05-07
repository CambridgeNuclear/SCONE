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
  !! Interface:
  !!   init -> initialises map from a dictionary
  !!   bins -> returns number of bins in dimension D. For D=0 returns number of all
  !!       bins in a map.
  !!   dimensions -> returns number of dimensions of a map.
  !!   getAxisName -> returns character with the type of the map
  !!   print -> prints information about the map to current block of the output file
  !!   binArrayShape -> returns a shape vector of the map
  !!
  type, public,abstract :: tallyMap
    private

  contains
    procedure(init),deferred        :: init
    procedure(bins),deferred        :: bins
    procedure(dimensions),deferred  :: dimensions
    procedure(getAxisName),deferred :: getAxisName
    procedure                       :: binArrayShape
    procedure(map),deferred         :: map
    procedure(print),deferred       :: print
    procedure                       :: kill
  end type tallyMap

  ! Procedures extendable in subclasses
  public :: kill


  abstract interface

    !!
    !! Initialise tallyMap from a dictionary
    !!
    !! Args:
    !!   dict [in] -> dictionary
    !!
    !! Errors:
    !!   Returns fatalError for invalid input
    !!
    subroutine init(self, dict)
      import :: tallyMap, &
                dictionary
      class(tallyMap), intent(inout) :: self
      class(dictionary), intent(in)  :: dict
    end subroutine init

    !!
    !! Return total number of bins in this division along Dimension D
    !! If D==0 returns total number of bins
    !!
    !! Args:
    !!   D [in] -> integer dimension
    !!
    !! Result:
    !!   Number of bins along dimension D. Returns total number of bins for D==0.
    !!   Returns 0 for invalid D (e.g. -ve)
    !!
    elemental function bins(self, D) result(N)
      import :: tallyMap, &
                shortInt
      class(tallyMap), intent(in)   :: self
      integer(shortInt), intent(in) :: D
      integer(shortInt)             :: N
    end function bins

    !!
    !! Returns number of dimensions
    !!
    !! Args:
    !!  None
    !!
    !! Result:
    !!   Integer giving number of dimensions in the map
    !!
    elemental function dimensions(self) result(D)
      import :: tallyMap, &
                shortInt
      class(tallyMap), intent(in)    :: self
      integer(shortInt)              :: D
    end function dimensions

    !!
    !! Return string that describes variable used to divide event space
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Left-adjusted, nameLen long character with the name of the axis type (e.g. x-coord)
    !!
    function getAxisName(self) result(name)
      import :: tallyMap, &
                nameLen
      class(tallyMap), intent(in) :: self
      character(nameLen)          :: name
    end function getAxisName

    !!
    !! Map particle to a single bin. Return 0 for particle out of division
    !!
    !! Args:
    !!   state [in] -> particleState to map
    !!
    !! Result:
    !!   Integer specifying the bin index for given particle state.
    !!   Returns 0 if the particle is outside the mapable range
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
    !! Add information about division axis to the output file
    !!
    !! Args:
    !!   out [inout] -> initialised outputFile
    !!
    subroutine print(self, out)
      import :: tallyMap, &
                outputFile
      class(tallyMap), intent(in)      :: self
      class(outputFile), intent(inout) :: out
    end subroutine print

  end interface

contains

  !!
  !! Returns array of sizes in each dimension
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   ShortInt array of size equal to number of dimensions.
  !!   Number of bins in each dimensions is listed in the array
  !!
  !! NOTE:
  !!   The result of this function is equivalent to the result of 'shape' intrinsic procedure
  !!   as if the tally map was a multi-dimensional matrix.
  !!
  pure function binArrayShape(self) result(sh)
    class(tallyMap), intent(in)                      :: self
    integer(shortInt),dimension(self % dimensions()) :: sh
    integer(shortInt)                                :: i

    do i=1,self % dimensions()
      sh(i) = self % bins(i)
    end do

  end function binArrayShape

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(tallyMap), intent(inout) :: self

  end subroutine kill



end module tallyMap_inter
