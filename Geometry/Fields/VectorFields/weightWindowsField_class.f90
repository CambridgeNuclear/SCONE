module weightWindowsField_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use dictionary_class,  only : dictionary
  use dictParser_func,   only : fileToDict
  use particle_class,    only : particle, particleState
  use field_inter,       only : field
  use vectorField_inter, only : vectorField

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: weightWindowsField_TptrCast

  !!
  !! Weight Windows Vector Field
  !!
  !! Returns a vector with lower, upper and survival weight
  !!
  !! Sample Dictionary Input:
  !!   varianceReduction { type weightWindowsField; file /home/WWfile.txt ;}
  !!
  !! Public Members:
  !!   net ->  map that lays over the geometry. Can be any type of Tally Map; it can
  !!           also include undefined bins if weight windows values are given
  !!   N   ->  total number of map bins
  !!   lowerW ->  array with the weight windows lower weights;
  !!   upperW ->  array with the weight windows upper weights;
  !!   constSurvival ->  constant to be multiplied to lower weight for getting
  !!                     roulette survival weight
  !!
  !! Interface:
  !!   vectorField interface
  !!
  type, public, extends(vectorField) :: weightWindowsField
    class(tallyMap), allocatable :: net
    integer(shortInt)            :: N
    real(defReal), dimension(:), allocatable :: lowerW
    real(defReal), dimension(:), allocatable :: upperW
    real(defReal) :: constSurvival = ZERO
  contains
    ! Superclass interface
    procedure :: init
    procedure :: kill
    procedure :: at
  end type weightWindowsField

contains

  !!
  !! Initialise from dictionary
  !!
  !! See field_inter for details
  !!
  subroutine init(self, dict)
    class(weightWindowsField), intent(inout) :: self
    class(dictionary), intent(in)            :: dict
    type(dictionary)                         :: dict2
    character(pathLen)                       :: path
    integer(shortInt), parameter  :: ALL = 0
    character(100), parameter :: Here = 'init (weightWindowsField_class.f90)'

    call dict % get(path,'file')

    ! Load dictionary
    call fileToDict(dict2, path)

    ! Initialise overlay map
    call new_tallyMap(self % net, dict2 % getDictPtr('map'))
    self % N = self % net % bins(ALL)

    ! Read weight values from file
    call dict2 % get(self % constSurvival, 'constSurvival')
    call dict2 % get(self % lowerW, 'wLower')
    call dict2 % get(self % upperW, 'wUpper')

    ! Check data consistency
    if (size(self % lowerW) /= self % N) then
      call fatalError(Here, 'Size of weight window net does not match the map provided')
    end if

    if (self % constSurvival <= ONE) call fatalError(Here, 'Survival constant must be bigger than one')
    if (any(self % lowerW * self % constSurvival > self % upperW)) then
      call fatalError(Here, 'Russian Roulette survival weight must be smaller &
                           & than the weight window upper weight')
    end if


  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(weightWindowsField), intent(inout) :: self

    call self % net % kill()
    deallocate(self % net)

  end subroutine kill

  !!
  !! Get value of the vector field given the phase-space location of a particle
  !!
  !! See vectorField_inter for details
  !!
  function at(self, p) result(val)
    class(weightWindowsField), intent(in) :: self
    class(particle), intent(inout)        :: p
    real(defReal), dimension(3)           :: val
    type(particleState)                   :: state
    integer(shortInt)                     :: binIdx

    ! Get current particle state
    state = p

    ! Read map bin index
    binIdx = self % net % map(state)

    ! Return if invalid bin index
    if (binIdx == 0) then
      val = ZERO
      return
    end if

    val(1) = self % lowerW(binIdx)
    val(2) = self % upperW(binIdx)
    val(3) = val(1) * self % constSurvival

  end function at

  !!
  !! Cast field pointer to weightWindowsField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of weightWindowsField
  !!   Pointer to source if source is weightWindowsField type
  !!
  pure function weightWindowsField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(weightWindowsField), pointer :: ptr

    select type (source)
      type is (weightWindowsField)
        ptr => source

      class default
        ptr => null()
    end select

  end function weightWindowsField_TptrCast


end module weightWindowsField_class
