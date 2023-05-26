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
    procedure(getEmittedRad), deferred     :: getEmittedRad
    procedure(getMaterialEnergy), deferred :: getMaterialEnergy
    procedure(updateProperties), deferred  :: updateProperties
    procedure(setTimeStep), deferred       :: setTimeStep
    procedure(setCalcType), deferred       :: setCalcType

  end type mgIMCDatabase

  abstract interface

    !!
    !! Return energy to be emitted during current time step
    !!
    !! Args:
    !!   matIdx [in] [optional] -> If provided, return the energy to be emitted from only matIdx
    !!                             Otherwise, return total energy to be emitted from all mats
    !!
    function getEmittedRad(self, matIdx) result(energy)
      import :: mgIMCDatabase, shortInt, defReal
      class(mgIMCDatabase), intent(in)        :: self
      integer(shortInt), intent(in), optional :: matIdx
      real(defReal)                           :: energy
    end function getEmittedRad

    !!
    !! Return material energy
    !!
    !! Args:
    !!   matIdx [in] [optional] -> If provided, return the energy of only matIdx
    !!                             Otherwise, return total energy of all mats
    !!
    function getMaterialEnergy(self, matIdx) result(energy)
      import :: mgIMCDatabase, shortInt, defReal
      class(mgIMCDatabase), intent(in)        :: self
      integer(shortInt), intent(in), optional :: matIdx
      real(defReal)                           :: energy
    end function getMaterialEnergy


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

    !!
    !! Tell each material if we are using IMC or ISMC
    !!
    subroutine setCalcType(self, type)
      import mgIMCDatabase, shortInt
      class(mgIMCDatabase), intent(inout) :: self
      integer(shortInt), intent(in)       :: type
    end subroutine setCalcType

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
