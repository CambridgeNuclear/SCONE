module mgIMCDatabase_inter

  use numPrecision

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

  contains
    procedure(getTotalEnergy), deferred   :: getTotalEnergy
    procedure(updateProperties), deferred :: updateProperties
    procedure(setTimeStep), deferred      :: setTimeStep

  end type mgIMCDatabase

  abstract interface

    !!
    !! Return total energy to be emitted during current time step
    !!
    function getTotalEnergy(self) result(energy)
      import :: mgIMCDatabase, defReal
      class(mgIMCDatabase), intent(in) :: self
      real(defReal)                    :: energy
    end function getTotalEnergy

    !!
    !! Update material properties based on energy absorbed during the time step
    !!
    subroutine updateProperties(self, tallyEnergy, printUpdates)
      import mgIMCDatabase, defReal, shortInt
      class(mgIMCDatabase), intent(inout)     :: self
      real(defReal), dimension(:), intent(in) :: tallyEnergy
      integer(shortInt), intent(in)           :: printUpdates
    end subroutine updateProperties

    !!
    !! Provide each material with time step to calculate initial fleck factor
    !!
    subroutine setTimeStep(self, deltaT)
      import mgIMCDatabase, defReal
      class(mgIMCDatabase), intent(inout) :: self
      real(defReal), intent(in)               :: deltaT
    end subroutine setTimeStep

  end interface

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
