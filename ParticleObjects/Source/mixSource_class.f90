module mixSource_class

  use numPrecision
  use errors_mod,        only : fatalError
  use genericProcedures, only : numToChar, linearSearchFloor
  use dictionary_class,  only : dictionary
  use RNG_class,         only : RNG
  use particle_class,    only : particleState

  ! Source
  use source_inter,      only : source

  ! Geometry
  use geometry_inter,    only : geometry

  implicit none
  private

  !!
  !! Public Pointer Type Cast
  !! Needed for initialisation in sourceFactory_func
  !!
  public mixSource_TptrCast

  !!
  !! Helper type to store polymorphic instances of sources
  !!
  type source_slot
    class(source), allocatable :: sourceItem
  end type source_slot

  !!
  !! Class describing a mixture of sources
  !!
  !! It builds and samples from a discrete cdf which source the particle should be
  !! sampled from. The cdf is built from the sources' weights, provided by the user.
  !!
  !! Private members:
  !!   sources -> array of sources
  !!   cdf     -> cdf used to select a source
  !!
  !! Interface:
  !!   init              -> initialise mixture of sources
  !!   sampleParticle    -> sample which source to sample the particle from
  !!   kill              -> terminate source
  !!
  !! Sample Dictionary Input:
  !!   source {
  !!       type mixSource;
  !!       sources (src1 src2);
  !!       weights (0.5 0.5);
  !!       src1 { < source 1 definition> }
  !!       src1 { < source 2 definition> }
  !!      }
  !!
  type, public, extends(source) :: mixSource
    type(source_slot), dimension(:), allocatable :: sources
    real(defReal), dimension(:), allocatable     :: cdf
  contains
    procedure :: init
    procedure :: sampleParticle
    procedure :: kill
  end type mixSource

contains

  !!
  !! Initialise from dictionary
  !!
  !! See source_inter for details
  !!
  !! Errors:
  !!   - when the number of weights doesn't correspond to the number of sources
  !!   - if any of the weights is zero
  !!
  subroutine init(self, dict, geom)
    class(mixSource), intent(inout)          :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    integer(shortInt)                        :: i, N
    real(defReal), dimension(:), allocatable :: weights
    character(100), parameter :: Here = 'init (mixSource_class.f90)'

    ! Get number of sources
    ! Note that sources have been allocated and initialised in sourceFactory
    N = size(self % sources)

    ! Retrieve weights relative to each source and check size
    call dict % get(weights, 'weights')
    if (size(weights) /= N) then
      call fatalError(Here, 'The number of sources and weights is not the same')
    end if
    if (any(weights <= ZERO)) call fatalError(Here, 'Some of the weights are not positive')

    ! Normalise weights
    weights = weights / sum(weights)

    ! Build cdf
    allocate(self % cdf(N + 1))
    self % cdf(1) = ZERO

    do i = 1, N
      self % cdf(i + 1) = self % cdf(i) + weights(i)
    end do

  end subroutine init

  !!
  !! Sample particle's phase space co-ordinates
  !!
  !!
  function sampleParticle(self, rand) result(p)
    class(mixSource), intent(inout) :: self
    class(RNG), intent(inout)       :: rand
    type(particleState)             :: p
    integer(shortInt)               :: idx

    idx = linearSearchFloor(self % cdf, rand % get())

    ! Sample particle from the sampled source
    p = self % sources(idx) % sourceItem % sampleParticle(rand)

  end function sampleParticle

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(mixSource), intent(inout) :: self

    ! Kill local components
    if (allocated(self % sources)) deallocate(self % sources)
    if (allocated(self % cdf)) deallocate(self % cdf)

  end subroutine kill

  !!
  !! Cast source pointer to mixSource pointer
  !!
  !! Args:
  !!   src [in] -> source pointer of class source
  !!
  !! Result:
  !!   Null if src is not of mixSource class
  !!   Target points to src if src is mixSource class
  !!
  pure function mixSource_TptrCast(src) result(ptr)
    class(source), target, intent(in) :: src
    type(mixSource), pointer          :: ptr

    select type(src)
      class is(mixSource)
        ptr => src

      class default
        ptr => null()
    end select

  end function mixSource_TptrCast

end module mixSource_class
