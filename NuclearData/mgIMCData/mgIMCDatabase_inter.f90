module mgIMCDatabase_inter

  ! Nuclear Data Interfaces & Objects
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: mgIMCDatabase_CptrCast

  !!
  !! An abstract class that groups all MG IMC Data objects
  !!
  !! It does nothing, It adds nothing,
  !! It just provides a common superclass for related classes
  !!
  type, public, abstract, extends(nuclearDatabase) :: mgIMCDatabase

  end type mgIMCDatabase

contains

  !!
  !! Cast nuclearDatabase pointer to mgIMCDatabase class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of mgIMCDatabase class
  !!   Target points to source if source is mgIMCDatabase class
  !!
  pure function mgIMCDatabase_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(mgIMCDatabase), pointer               :: ptr

    select type(source)
      class is(mgIMCDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function mgIMCDatabase_CptrCast


end module mgIMCDatabase_inter
