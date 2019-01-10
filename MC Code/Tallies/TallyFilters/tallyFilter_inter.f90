module tallyFilter_inter

  use numPrecision
  use particle_class,   only : particleState
  use dictionary_class, only : dictionary

  implicit none
  private

  !!
  !! Abstract interface of tallyFilters
  !!
  !! All tallyFilters given a particle state return .true. or .false.
  !! This determines whether an event should be tallied or not
  !!
  type, public,abstract :: tallyFilter
    private
  contains
    procedure(init),deferred   :: init
    procedure(filter),deferred :: filter
  end type tallyFilter

  abstract interface
    !!
    !! Initialise filter from dictionary
    !!
    subroutine init(self, dict)
      import :: tallyFilter, &
                dictionary
      class(tallyFilter), intent(inout) :: self
      class(dictionary), intent(in)     :: dict
    end subroutine init

    !!
    !! Return .true. if state passes filter test
    !! Return .false. otherwise or if test is undefined
    !!
    elemental function filter(self, state) result(passed)
      import :: tallyFilter,  &
                particleState,&
                defBool
      class(tallyFilter), intent(in)   :: self
      class(particleState), intent(in) :: state
      logical(defBool)                 :: passed
    end function filter

  end interface

end module tallyFilter_inter
