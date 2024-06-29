module mgNeutronDatabase_inter

  use numPrecision
  
  ! Nuclear Data Interfaces & Objects
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: mgNeutronDatabase_CptrCast

  !!
  !! An abstract class that groups all MG Neutron Data objects
  !!
  !! It does nothing, It adds nothing, (other than give the number of groups)
  !! It just provides a common superclass for related classes
  !! and holds the number of energy groups
  !!
  !! Public members: 
  !!   nG -> number of energy groups
  !!
  type, public, abstract, extends(nuclearDatabase) :: mgNeutronDatabase
    integer(shortInt)  :: nG = 0
  contains
  end type mgNeutronDatabase

  abstract interface

    !!
    !! Returns the number of energy groups used
    !!
    !! Args:
    !!   None
    !!
    !! Result:
    !!   Integer number of energy groups
    !!
    pure function nGroups(self) result(ng)
      import :: mgNeutronDatabase, shortInt
      class(mgNeutronDatabase), intent(in) :: self
      integer(shortInt)                    :: ng
    end function nGroups


  end interface

contains

  !!
  !! Cast nuclearDatabase pointer to mgNeutronDatabase class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of mgNeutronDatabase class
  !!   Target points to source if source is mgNeutronDatabase class
  !!
  pure function mgNeutronDatabase_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(mgNeutronDatabase), pointer           :: ptr

    select type(source)
      class is(mgNeutronDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function mgNeutronDatabase_CptrCast


end module mgNeutronDatabase_inter
