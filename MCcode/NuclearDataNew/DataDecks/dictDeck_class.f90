module dictDeck_class

  use numPrecision
  use dictionary_class, only : dictionary
  use dataDeck_inter,   only : dataDeck

  implicit none
  private

  !!
  !! SCONE dictionary based data deck
  !!
  !! Just a simple wrapper around the pointer to the dictionary
  !!
  !! Public Members:
  !!   dict -> pointer to the dictionary
  !!
  type, public, extends(dataDeck) :: dictDeck
    class(dictionary),pointer :: dict => null()
  contains
    procedure :: myType
  end type dictDeck

contains

  !!
  !! Return String with type name
  !!
  !! Returns:
  !!   An allocatable character with type name without trailing blanks
  !!
  pure function myType(self) result(type)
    class(dictDeck), intent(in) :: self
    character(:),allocatable    :: type

    type = 'dictDeck'
  end function myType

end module dictDeck_class
