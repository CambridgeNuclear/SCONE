module dataDeck_inter


  implicit none
  private

  !!
  !! Abstract Nuclear Data Library Access interface
  !!
  !! Name is chosen for reto-associations
  !!
  !! This abstract type creates a family of objects that can be used to build
  !! Nuclear Data objects for execution (e.g. reactions, nuclides, materials etc.)
  !! This allow to create an easy interface for building nuclear data object from
  !! multiple library formats (ACE, HDF5, SCONE Dictionary etc.)
  !!
  !! It is up to every nuclearData object to decide how it can be built!
  !!
  !! Interface:
  !!   myType -> Return character with type name. Usefull for error messages!
  !!
  type, public, abstract :: dataDeck

  contains
    procedure(myType),deferred :: myType
  end type dataDeck

  abstract interface
    !!
    !! Return String with type name
    !!
    !! Returns:
    !!   An allocatable character with type name without trailink blanks
    !!
    pure function myType(self) result(type)
      import :: dataDeck
      class(dataDeck), intent(in) :: self
      character(:),allocatable    :: type
    end function myType

  end interface


contains


    
end module dataDeck_inter
