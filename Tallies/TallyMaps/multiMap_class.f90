module multiMap_class

  use numPrecision
  use universalVariables, only : INF
  use genericProcedures,  only : fatalError
  use dictionary_class,   only : dictionary
  use outputFile_class,   only : outputFile
  use particle_class,     only : particleState

  use tallyMap_inter,         only : tallyMap
  use tallyMap1D_inter,       only : tallyMap1D
  use tallyMap1DFactory_func, only : new_tallyMap1D

  implicit none
  private

  !!
  !! Helper type to store polymorphic instances of 1D tallyMaps
  !!
  type map1D
    class(tallyMap1D), allocatable :: slot
  end type map1D

  !!
  !! Multi Dimensional Tally Map
  !!
  !! Combines arbitrary number of 1D maps into a single multidimensional map
  !!
  !! Private Members:
  !!   maps -> tallyMaps. Variation in column-major order. Leftmost fastest
  !!   multi -> bin location multipliers. Stride of each change in a row.
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
    type(map1D), dimension(:), allocatable        :: maps
    integer(shortInt), dimension(:), allocatable  :: multi

  contains
    procedure :: init
    procedure :: bins
    procedure :: dimensions
    procedure :: getAxisName
    procedure :: map
    procedure :: distance
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
    integer(shortInt)                             :: i, mul
!    character(100), parameter :: Here = 'init (multiMap_class.f90)'

    ! Read order of maps requested
    call dict % get(mapNames, 'maps')

    ! Allocate space
    allocate(self % maps(size(mapNames)))
    allocate(self % multi(size(mapNames)))

    ! Build maps
    do i = 1, size(self % maps)
      call new_tallyMap1D(self % maps(i) % slot, dict % getDictPtr(mapNames(i)))
    end do

    ! Calculate multipliers
    mul = 1
    do i = 1,size(self % maps)
      self % multi(i) = mul
      mul = mul * self % maps(i) % slot % bins(1)
    end do

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
        N = N * self % maps(i) % slot % bins(1)
      end do

    else if( 0 < D .and. D <= size(self % maps) ) then
      N = self % maps(D) % slot % bins(1)

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
  elemental function map(self, state) result(idx)
    class(multiMap), intent(in)      :: self
    class(particleState), intent(in) :: state
    integer(shortInt)                :: idx, binIdx
    integer(shortInt)                :: i

    idx = 1
    do i=1,size(self % maps)
      binIdx = self % maps(i) % slot % map(state)

      ! Short-circuit evaluation if out-of division
      if(binIdx == 0) then
        idx = 0
        return
      end if

      idx = idx + (binIdx-1) * self % multi(i)

    end do

  end function map
  
  !!
  !! Compute particle's distance to boundary
  !!
  !! See tallyMap for specification
  !!
  elemental function distance(self, state) result(d)
    class(multiMap), intent(in)      :: self
    class(particleState), intent(in) :: state
    real(defReal)                    :: d
    integer(shortInt)                :: i
    real(defReal)                    :: d0

    d = INF
    
    do i = 1, size(self % maps)
      d0 = self % maps(i) % slot % distance(state)
      d = min(d, d0)
    end do

  end function distance

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
      call out % addValue( self % maps(i) % slot % getAxisName() )

    end do
    call out % endArray()

    ! Print individual maps information
    do i = 1,size(self % maps)
      call self % maps(i) % slot % print(out)

    end do

  end subroutine print

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(multiMap), intent(inout) :: self
    integer(shortInt)              :: i

    ! Kill maps
    if(allocated(self % maps)) then
      do i=1,size(self % maps)
        call self % maps(i) % slot % kill()
      end do
      deallocate(self % maps)
    end if

    ! Free rest of space space
    if(allocated(self % multi)) deallocate(self % multi)

  end subroutine kill

end module multiMap_class
