module tallyMap1D_inter

  use numPrecision
  use particle_class,   only : particleState
  use dictionary_class, only : dictionary
  use outputFile_class, only : outputFile
  use tallyMap_inter,   only : tallyMap

  implicit none
  private

  !!
  !! Abstract interface to separate simple 1D maps from more complex types
  !!
  type, public,extends(tallyMap),abstract :: tallyMap1D
    private
  contains
    procedure :: dimensions  ! Return number of dimensions
  end type tallyMap1D

contains

  !!
  !! Return number of dimensions for 1D map (1)
  !!
  elemental function dimensions(self) result(D)
    class(tallyMap1D), intent(in) :: self
    integer(shortInt)             :: D

    D = 1

  end function dimensions
    
end module tallyMap1D_inter
