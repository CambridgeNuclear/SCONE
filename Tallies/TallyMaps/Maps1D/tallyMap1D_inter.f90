module tallyMap1D_inter

  use numPrecision
  use particle_class,   only : particleState
  use dictionary_class, only : dictionary
  use outputFile_class, only : outputFile
  use tallyMap_inter,   only : tallyMap, kill_super => kill

  implicit none
  private

  !!
  !! Abstract interface to separate simple 1D maps from more complex types
  !!
  !! Interface:
  !!   See tallyMap
  !!
  type, public,extends(tallyMap),abstract :: tallyMap1D
    private
  contains
    procedure                       :: dimensions
    procedure                       :: kill
  end type tallyMap1D

  ! Procedures extendable in subclasses
  public :: kill

contains

  !!
  !! Return number of dimensions for 1D map (1)
  !!
  !! See tallyMap for specification
  !!
  elemental function dimensions(self) result(D)
    class(tallyMap1D), intent(in) :: self
    integer(shortInt)             :: D

    D = 1

  end function dimensions
    
  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(tallyMap1D), intent(inout) :: self

    ! Call superclass procedure
    call kill_super(self)

  end subroutine kill

end module tallyMap1D_inter
