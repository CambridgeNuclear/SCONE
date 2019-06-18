module nuclearDatabase_inter

  use numPrecision
  use particle_class, only : particle
  use charMap_class,  only : charMap

  implicit none
  private

  !!
  !! Nuclear Database contains nuclear interaction data
  !!
  !! This is a top abstract class for all types of nuclear data.
  !! So nuclear data objects for CE or MG neutron, photons or other particles are
  !! subclasses of this type.
  !!
  !! Interface:
  !!   getTransMatXS -> returns transport Material XS given a particle
  !!   getTotalMatXS -> returns total Material XS fiven a particle
  !!   getMajorantXS -> returns majorant XS given particle and list of active materials
  !!   matNamesMap   -> returns pointer to map of material names to matIdx
  !!   kill          -> return to uninitialised state, clean memory
  !!
  type, public,abstract :: nuclearDatabase
  contains
    procedure(getTransMatXS), deferred :: getTransMatXS
    procedure(getTotalMatXS), deferred :: getTotalMatXS
    procedure(getMajorantXS), deferred :: getMajorantXS
    procedure(matNamesMap), deferred   :: matNamesMap
    procedure(kill), deferred          :: kill
  end type nuclearDatabase

  abstract interface

    !!
    !! Return value of Material Transport XS for a particle
    !!
    !! Reads all relevalnt state information from the particle (e.g. E or G)
    !! Usually is the same as material total XS, but it may be not always the case
    !!
    !! Args:
    !!   p [in]      -> particle at a given state
    !!   matIdx [in] -> Material index
    !!
    !! Result:
    !!   Value of a material transport XS [1/cm]
    !!
    !! Errors:
    !!   Undefined behaviour if the state of the particle is invalid e.g. -ve energy
    !!   Undefined behavior if matIdx does not correspond to a defined material
    !!
    function getTransMatXS(self, p, matIdx) result(xs)
      import :: nuclearDatabase, particle, shortInt, defReal
      class(nuclearDatabase), intent(inout) :: self
      class(particle), intent(in)           :: p
      integer(shortInt), intent(in)         :: matIdx
      real(defReal)                         :: xs
    end function getTransMatXS

    !!
    !! Return value of Material Total XS for a particle
    !!
    !! Reads all relevalnt state information from the particle (e.g. E or G)
    !!
    !! Args:
    !!   p [in]      -> particle at a given state
    !!   matIdx [in] -> Material index
    !!
    !! Result:
    !!   Value of a material total XS [1/cm]
    !!
    !! Errors:
    !!   Undefined behaviour if the state of the particle is invalid e.g. -ve energy
    !!   Undefined behavior if matIdx does not correspond to a defined material
    !!
    function getTotalMatXS(self, p, matIdx) result(xs)
      import :: nuclearDatabase, particle, shortInt, defReal
      class(nuclearDatabase), intent(inout) :: self
      class(particle), intent(in)           :: p
      integer(shortInt), intent(in)         :: matIdx
      real(defReal)                         :: xs
    end function getTotalMatXS

    !!
    !! Return value of Majorant XS for a particle
    !!
    !! Reads all relevalnt state information from the particle (e.g. E or G)
    !! Majorant XS is the largest of TRANSPORT XSs for ACTIVE materials
    !!
    !! Args:
    !!   p [in] -> particle at a given state
    !!
    !! Result:
    !!   Value of a majorant XS [1/cm]
    !!
    !! Errors:
    !!   Undefined behaviour if the state of the particle is invalid e.g. -ve energy
    !!
    function getMajorantXS(self, p) result(xs)
      import :: nuclearDatabase, particle, shortInt, defReal
      class(nuclearDatabase), intent(inout) :: self
      class(particle), intent(in)           :: p
      real(defReal)                         :: xs
    end function getMajorantXS

    !!
    !! Return pointer to material names map
    !!
    !! Result:
    !!   A pointer to map that translates material names -> matIdx
    !!
    !! Errors:
    !!   fatalError if the map is not initialised
    !!
    function matNamesMap(self) result(map)
      import :: nuclearDatabase, charMap
      class(nuclearDatabase), intent(in) :: self
      type(charMap), pointer             :: map
    end function matNamesMap

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: nuclearDatabase
      class(nuclearDatabase), intent(inout) :: self
    end subroutine kill
  end interface
    
end module nuclearDatabase_inter
