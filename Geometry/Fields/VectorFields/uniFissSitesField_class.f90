module uniFissSitesField_class

  use numPrecision
  use genericProcedures,     only : fatalError, numToChar
  use dictionary_class,      only : dictionary
  use particle_class,        only : particle, particleState
  use particleDungeon_class, only : particleDungeon
  use field_inter,           only : field
  use vectorField_inter,     only : vectorField

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: uniFissSitesField_TptrCast

  !!
  !! Uniform Fission Sites Field
  !!
  !! Returns a 2D vector with:
  !! - fraction of fissionable material volume occupied by the required cell
  !! - percentage of fission sites in the required cell
  !! NOTE: as things are now, results are correct only if all the spatial bins
  !!       include the same amount of fissionable material (are of equal size) !!!!!
  !!
  !! Sample Dictionary Input:
  !!   uniformFissionSites { map { <map definition> }}
  !!
  !! Public Members:
  !!   net ->  map that lays over the geometry. It should be a spatial map, an
  !!           energy map wouldn't make much sense!
  !!   N   ->  total number of map bins
  !!   sourceFraction -> array with the percentage of fission sites in each bin
  !!   buildSource    -> array used to 'tally' fission sites
  !!
  !! Interface:
  !!   vectorField interface
  !!
  type, public, extends(vectorField) :: uniFissSitesField
    class(tallyMap), allocatable :: net
    integer(shortInt)            :: N = 0
    real(defReal), dimension(:), allocatable     :: sourceFraction
    real(defReal), dimension(:), allocatable     :: buildSource
  contains
    ! Superclass interface
    procedure :: init
    procedure :: kill
    procedure :: at
    procedure :: storeFS
    procedure :: updateMap
  end type uniFissSitesField

contains

  !!
  !! Initialise from dictionary
  !!
  !! See field_inter for details
  !!
  subroutine init(self, dict)
    class(uniFissSitesField), intent(inout) :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt), parameter  :: ALL = 0
    character(100), parameter     :: Here = 'init (uniFissSitesField_class.f90)'

    ! Initialise overlay map
    call new_tallyMap(self % net, dict % getDictPtr('map'))
    self % N = self % net % bins(ALL)

    allocate(self % sourceFraction(self % N), self % buildSource(self % N),self % weights(self % N))
    self % sourceFraction = ONE/self % N
    self % buildSource = ZERO

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(uniFissSitesField), intent(inout) :: self

    call self % net % kill()
    deallocate(self % net)
    deallocate(self % sourceFraction)
    deallocate(self % buildSource)

    self % N = 0

  end subroutine kill

  !!
  !! Get value of the vector field given the phase-space location of a particle
  !!
  !! See vectorField_inter for details
  !!
  function at(self, p) result(val)
    class(uniFissSitesField), intent(in) :: self
    class(particle), intent(inout)       :: p
    real(defReal), dimension(3)          :: val
    type(particleState)                  :: state
    integer(shortInt)                    :: binIdx

    ! Get current particle state
    state = p

    ! Read map bin index
    binIdx = self % net % map(state)

    ! Return if invalid bin index
    if (binIdx == 0) then
      val = ONE
      return
    end if

    val(1) = ONE / self % N
    val(2) = self % sourceFraction(binIdx)
    val(3) = ZERO

  end function at

  !!
  !! Store the fission sites generated in a vector
  !!
  !! Args:
  !! state [in] -> particle state of the fission site
  !!
  subroutine storeFS(self, state)
    class(uniFissSitesField), intent(inout) :: self
    type(particleState), intent(in)         :: state
    integer(shortInt)                       :: idx

    idx = self % net % map(state)
    if (idx == 0) return
    ! Add fission sites where appropriate
    self % buildSource(idx) = self % buildSource(idx) + state % wgt

  end subroutine storeFS

  !!
  !! Calculates fission site probability distribution
  !! It accounts for possible ares of the map having zero events
  !!
  subroutine updateMap(self)
    class(uniFissSitesField), intent(inout) :: self
    integer(shortInt) :: i

    ! Eliminate zeros in the distribution
    do i = 1, self % N
      if (self % buildSource(i) == ZERO) self % buildSource(i) = ONE
    end do
    ! Normalise to calculate probability
    self % sourceFraction = self % buildSource / sum(self % buildSource)

    self % buildSource = ZERO

  end subroutine updateMap

  !!
  !! Cast field pointer to uniFissSitesField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of uniFissSitesField
  !!   Pointer to source if source is uniFissSitesField type
  !!
  pure function uniFissSitesField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(uniFissSitesField), pointer  :: ptr

    select type (source)
      type is (uniFissSitesField)
        ptr => source

      class default
        ptr => null()
    end select

  end function uniFissSitesField_TptrCast


end module uniFissSitesField_class
