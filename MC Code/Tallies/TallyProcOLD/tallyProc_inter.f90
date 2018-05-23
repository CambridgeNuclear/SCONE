module tallyProc_inter

  use numPrecision
  use particle_class,     only : particle, phaseCoord
  use tallyCounter_class, only : tallyCounter

  implicit none
  private

  !!
  !! Abstract interface for a tally processor
  !! Given event information returns pointer to tallyCounters that apply
  !!
  type, public,abstract :: tallyProc
    private
  contains
    procedure(filterCollision), deferred :: filterCollision
    procedure(numOfBins), deferred       :: numOfBins
    generic                              :: assignBins  => assignBins_ptr, assignBins_loc
    procedure(assignBins_ptr),deferred   :: assignBins_ptr
    procedure(assignBins_loc),deferred   :: assignBins_loc
  end type tallyProc
    
  abstract interface

    !!
    !! Given information about the collision and response code
    !! Return pointer to bins that need to be scored
    !!
    function filterCollision(self,pre,post,MT,muL,resCode) result (bins)
      import :: tallyProc, &
                particle ,&
                phaseCoord ,&
                defReal, &
                shortInt, &
                tallyCounter
      class(tallyProc), intent(in)            :: self
      class(phaseCoord), intent(in)           :: pre
      class(particle), intent(in)             :: post
      integer(shortInt), intent(in)           :: MT
      real(defReal), intent(in)               :: muL
      integer(shortInt), intent(in)           :: resCode
      type(tallyCounter),dimension(:),pointer :: bins
    end function filterCollision

    !!
    !! Return number of bins tallyProc uses to store all data
    !!
    function numOfBins(self) result(N)
      import :: tallyProc, &
                shortInt
      class(tallyProc), intent(in) :: self
      integer(shortInt)            :: N
    end function numOfBins

    !!
    !! Given pointer to array of bins of size equal to tallyProc numOfBins
    !! Use the bins provided by a pointer for tallying
    !!
    subroutine assignBins_ptr(self,binMem)
      import :: tallyProc, &
                tallyCounter
      class(tallyProc),intent(inout)                     :: self
      type(tallyCounter),dimension(:),pointer,intent(in) :: binMem
    end subroutine assignBins_ptr

    !!
    !! Allocate space for tallying bins
    !!
    subroutine assignBins_loc(self)
      import :: tallyProc
      class(tallyProc), intent(inout) :: self
    end subroutine assignBins_loc

  end interface

end module tallyProc_inter
