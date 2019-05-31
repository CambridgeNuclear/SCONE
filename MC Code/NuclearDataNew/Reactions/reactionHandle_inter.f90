module reactionHandle_inter

  use numPrecision
  use dataDeck_inter, only : dataDeck

  implicit none
  private

  !!
  !! Abstract reaction handle
  !!
  !! Does nothing or has nothing by itself. It only exists to group all types or
  !! reaction interfaces into a single object fimaly
  !!
  type, public,abstract :: reactionHandle
    private
  contains
    procedure(init),deferred :: init
    procedure(kill),deferred :: kill
  end type reactionHandle

  abstract interface
    !!
    !! Initialise reaction object
    !!
    !! Args:
    !!   data [in] -> A dataDeck object with data required for build
    !!   MT [in]   -> MT number of reaction to be build
    !!
    !! Errors:
    !!   May return fatalError if type of dataDeck is not supported or
    !!   if build fails for any reason
    !!
    subroutine init(self, data, MT)
      import :: reactionHandle, dataDeck, shortInt
      class(reactionHandle), intent(inout) :: self
      class(dataDeck), intent(inout)       :: data
      integer(shortInt), intent(in)        :: MT
    end subroutine init

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: reactionHandle
      class(reactionHandle), intent(inout) :: self
    end subroutine kill

  end interface


end module reactionHandle_inter
