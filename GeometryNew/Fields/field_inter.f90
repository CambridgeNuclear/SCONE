module field_inter

  use numPrecision
  use dictionary_class, only : dictionary

  implicit none
  private

  !!
  !! Abstract class to combine scalar & vector fields into single family
  !!
  !! Does nothing by itself
  !!
  !! Interface:
  !!   init -> Initialise from dictionary
  !!   kill -> Return to uninitialised state
  !!
  type, public, abstract :: field
  contains
    procedure(init), deferred :: init
    procedure(kill), deferred :: kill
  end type field

  abstract interface

    !!
    !! Initialise field from dictionary
    !!
    !! Args:
    !!   dict [in] -> Dictionary with input data
    !!
    subroutine init(self, dict)
      import :: field, dictionary
      class(field), intent(inout)   :: self
      class(dictionary), intent(in) :: dict
    end subroutine init

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: field
      class(field), intent(inout) :: self
    end subroutine kill
  end interface

end module field_inter
