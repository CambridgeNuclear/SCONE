module multiMap_class

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary
  use outputFile_class,  only : outputFile
  use particle_class,    only : particleState

  use tallyMap_inter,     only : tallyMap
  use tallyMapSlot_class, only : tallyMapSlot

  implicit none
  private

  !!
  !! Multi Dimensional Tally Map
  !!
  !! Combines arbitrary number of 1D maps into a single multidimensional map
  !!
  !! Private Members:
  !!
  !! Interface:
  !!   See tallyMap
  !!
  !! Sample Input Dictionary:
  !!
  !! myMap {
  !!   type multiMap;
  !!   maps (map1 map2);
  !!   map1 { <1D map definition> }
  !!   map2 { <1D map definition> }
  !!   map3 { <1D map definition> }
  !! }
  !!
  type, public, extends(tallyMap) :: multiMap
    private
    type(tallyMapSlot), dimension(:), allocatable :: maps

  contains
    procedure :: init
    procedure :: bins
    procedure :: dimensions
    procedure :: getAxisName
    procedure :: map
    procedure :: print
    procedure :: kill
  end type multiMap

contains

  !!
  !! Initialise tallyMap from a dictionary
  !!
  !! See tallyMap for specification.
  !!
  subroutine init(self, dict)
    class(multiMap), intent(inout)                :: self
    class(dictionary), intent(in)                 :: dict
    character(nameLen), dimension(:), allocatable :: mapNames
    integer(shortInt)                             :: i
    character(100), parameter :: Here = 'init (multiMap_class.f90)'

    ! Read order of maps requested
    call dict % get(mapNames, 'maps')

    ! Allocate space
    allocate(self % maps(size(mapNames)))

    ! Build maps
    do i=1,size(self % maps)
      call self % maps(i) % init( dict % getDictPtr(mapNames(i)))
    end do

    ! Verify that all maps are 1D
    if(any(self % maps % dimensions() /= 1)) then
      call fatalError(Here, 'Multidimensional maps are not allowed in multiMap!')
    end if

  end subroutine init

  !!
  !! Return total number of bins in this division along Dimension D
  !!
  !! See tallyMap for specification.
  !!
  elemental function bins(self, D) result(N)
    class(multiMap), intent(in)   :: self
    integer(shortInt), intent(in) :: D
    integer(shortInt)             :: N
    integer(shortInt)             :: i

    if(D == 0) then
      ! Perform multiplicative reduction over all dimensions
      N = 1
      do i = 1,size(self % maps)
        N = N * self % maps(i) % bins(1)
      end do

    else if(D <= size(self % maps)) then
      N = self % maps(D) % bins(1)

    else
      N = 0

    end if

  end function bins

  !!
  !! Return number of dimensions
  !!
  !! See tallyMap for specification.
  !!
  elemental function dimensions(self) result(D)
    class(multiMap), intent(in)    :: self
    integer(shortInt)              :: D

    D = size(self % maps)

  end function dimensions

  !!
  !! Return string that describes variable used to divide event space
  !!
  !! See tallyMap for specification
  !!
  function getAxisName(self) result(name)
    class(multiMap), intent(in) :: self
    character(nameLen)          :: name

    name ='multiMap'

  end function getAxisName

  !!
  !! Map particle to a single bin. Return 0 for particle out of division
  !!
  !! See tallyMap for specification.
  !!
  elemental function map(self,state) result(idx)
    class(multiMap), intent(in)      :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx

    idx = 1

  end function map

  !!
  !! Add information about division axis to the output file
  !!
  !! See tallyMap for specification.
  !!
  subroutine print(self,out)
    class(multiMap), intent(in)      :: self
    class(outputFile), intent(inout) :: out
    character(nameLen)               :: name
    integer(shortInt)                :: i

    ! Print map order information
    name = 'multiMapOrder'

    call out % startArray(name, [size(self % maps)])
    do i = 1,size(self % maps)
      call out % addValue( self % maps(i) % getAxisName() )

    end do
    call out % endArray()

    ! Print individual maps information
    do i = 1,size(self % maps)
      call self % maps(i) % print(out)

    end do

  end subroutine print

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(multiMap), intent(inout) :: self

    ! Kill individual maps

    ! Free space
    deallocate(self % maps)

  end subroutine kill

end module multiMap_class
