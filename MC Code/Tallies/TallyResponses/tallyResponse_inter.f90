module tallyResponse_inter

  use numPrecision
  use dictionary_class, only : dictionary
  use particle_class,   only : particle

  implicit none
  private


  !!
  !! Abstract interface for all tallyResponses
  !! Very simple class, which given a particle returns a real number
  !! Real number is used to weight FLUX sample to score reaction rates etc.
  !! Returns only scalar to move all logic for dealing with multiple responces to tallyClerks.
  !! Thus tallyResponses should be quick and easy to write
  !!
  type, public,abstract :: tallyResponse
    private

  contains
    procedure(init), deferred :: init
    procedure(get), deferred  :: get

  end type tallyResponse

  abstract interface

    !!
    !! Initialise Response from dictionary
    !!
    subroutine init(self, dict)
      import :: tallyResponse, &
                dictionary
      class(tallyResponse), intent(inout) :: self
      class(dictionary), intent(in)       :: dict

    end subroutine init

    !!
    !! Get value of response for particle p
    !!
    elemental function get(self, p) result(val)
      import :: tallyResponse, &
                particle, &
                defReal
      class(tallyResponse), intent(in) :: self
      class(particle), intent(in)      :: p
      real(defReal)                    :: vale
    end function get
  end interface
    
end module tallyResponse_inter
