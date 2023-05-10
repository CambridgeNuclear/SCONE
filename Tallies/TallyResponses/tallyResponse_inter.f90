module tallyResponse_inter

  use numPrecision
  use dictionary_class,      only : dictionary
  use particle_class,        only : particle

  ! Nuclear Data interface
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private


  !!
  !! Abstract interface for all tallyResponses
  !!
  !! Very simple class, which given a particle returns a real number
  !! Real number is used to weight FLUX sample to score reaction rates etc.
  !! Returns only scalar to move all logic for dealing with multiple responces to tallyClerks.
  !! Thus tallyResponses should be quick and easy to write
  !!
  !! Interface:
  !!   init -> Initialise
  !!   get  -> Get value of the response
  !!   kill -> Return to uninitialised state
  !!
  type, public, abstract :: tallyResponse
    private

  contains
    procedure(init), deferred :: init
    procedure(get), deferred  :: get
    procedure(kill), deferred :: kill

  end type tallyResponse

  abstract interface

    !!
    !! Initialise Response from dictionary
    !!
    !! Args:
    !!   dict [in] -> DIctionary with the data
    !!
    !! Errors:
    !!   Depend on specific implementation.
    !!   fatalError if there is a mistake in definition
    !!
    subroutine init(self, dict)
      import :: tallyResponse, &
                dictionary
      class(tallyResponse), intent(inout) :: self
      class(dictionary), intent(in)       :: dict

    end subroutine init

    !!
    !! Get value of response
    !!
    !! Args:
    !!   p [in]         -> Particle to provide state
    !!   xsData [inout] -> Nuclear Database used by the particle
    !!
    !! Result:
    !!   Value of the response for particle p
    !!
    !! Errors:
    !!   Depend on specific implementation
    !!
    function get(self, p, xsData) result(val)
      import :: tallyResponse, &
                particle, &
                nuclearDatabase, &
                defReal
      class(tallyResponse), intent(in)      :: self
      class(particle), intent(in)           :: p
      class(nuclearDatabase), intent(inout) :: xsData
      real(defReal)                         :: val
    end function get

    !!
    !! Return to uninitialised state
    !!
    !! Args:
    !!   None
    !!
    !! Errors:
    !!   None
    !!
    elemental subroutine kill(self)
      import :: tallyResponse
      class(tallyResponse), intent(inout) :: self
    end subroutine kill

  end interface

end module tallyResponse_inter
