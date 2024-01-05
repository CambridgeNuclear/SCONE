module mgIMCDatabase_inter

  use numPrecision
  use RNG_class,        only : RNG

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
    procedure(getEmittedRad), deferred       :: getEmittedRad
    procedure(getMaterialEnergy), deferred   :: getMaterialEnergy
    procedure(updateProperties), deferred    :: updateProperties
    procedure(setCalcType), deferred         :: setCalcType
    procedure(sampleEnergyGroup), deferred   :: sampleEnergyGroup
    procedure(sampleTransformTime), deferred :: sampleTransformTime

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
      import :: mgIMCDatabase, defReal, shortInt
      class(mgIMCDatabase), intent(inout)     :: self
      real(defReal), dimension(:), intent(in) :: tallyEnergy
      integer(shortInt), intent(in)           :: printUpdates
    end subroutine updateProperties

    !!
    !! Tell each material if we are using IMC or ISMC
    !!
    subroutine setCalcType(self, type)
      import :: mgIMCDatabase, shortInt
      class(mgIMCDatabase), intent(inout) :: self
      integer(shortInt), intent(in)       :: type
    end subroutine setCalcType

    !!
    !! Sample energy group of a particle emitted from material matIdx
    !!
    !! Args:
    !!   matIdx [in] -> index of material to sample from
    !!   rand   [in] -> RNG
    !!
    !! Result:
    !!   G -> energy group of sampled particle
    !!
    function sampleEnergyGroup(self, matIdx, rand) result(G)
      import :: mgIMCDatabase, shortInt, RNG
      class(mgIMCDatabase), intent(inout) :: self
      integer(shortInt), intent(in)       :: matIdx
      class(RNG), intent(inout)           :: rand
      integer(shortInt)                   :: G
    end function sampleEnergyGroup

    !!
    !! Sample the time taken for a material particle to transform into a photon
    !! Used for ISMC only
    !!
    function sampleTransformTime(self, matIdx, rand) result(t)
      import :: mgIMCDatabase, shortInt, RNG, defReal
      class(mgIMCDatabase), intent(inout) :: self
      integer(shortInt), intent(in)       :: matIdx
      class(RNG), intent(inout)           :: rand
      real(defReal)                       :: t
    end function sampleTransformTime

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
