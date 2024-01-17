module nuclearDatabase_inter

  use numPrecision
  use dictionary_class, only : dictionary
  use particle_class,   only : particle
  use charMap_class,    only : charMap
  use RNG_class,        only : RNG

  ! Nuclear Data Handles
  use nuclideHandle_inter,  only : nuclideHandle
  use materialHandle_inter, only : materialHandle
  use reactionHandle_inter, only : reactionHandle

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
  !!   getMaterial   -> returns a pointer to a material handle for given matIdx
  !!   getNuclide    -> returns a pointer to a nuclide handle for given nucIdx
  !!   getReaction   -> returns a pointer to a reaction for given matidx or nucIdx and MT number
  !!   initMajorant  -> initialises the majorant cross section for delta tracking
  !!   kill          -> return to uninitialised state, clean memory
  !!
  type, public,abstract :: nuclearDatabase
  contains
    procedure(init), deferred          :: init
    procedure(activate), deferred      :: activate
    procedure(getTransMatXS), deferred :: getTransMatXS
    procedure(getTotalMatXS), deferred :: getTotalMatXS
    procedure(getMajorantXS), deferred :: getMajorantXS
    procedure(matNamesMap), deferred   :: matNamesMap
    procedure(getMaterial), deferred   :: getMaterial
    procedure(getNuclide), deferred    :: getNuclide
    procedure(getReaction), deferred   :: getReaction
    procedure(initMajorant), deferred  :: initMajorant
    procedure(kill), deferred          :: kill
  end type nuclearDatabase

  abstract interface
    !!
    !! Initialise Database from dictionary and pointer to self
    !!
    !! Args:
    !!   dict   [in] -> Dictionary with the settings
    !!   ptr    [in] -> Pointer to self (of class nuclearDatabase)
    !!   silent [in] -> Optional. If set to .true. disables console output
    !!
    !! Errors
    !!   FatalError is ptr is not assosiated with self
    !!
    subroutine init(self, dict, ptr, silent)
      import :: nuclearDatabase, dictionary, defBool
      class(nuclearDatabase), target, intent(inout) :: self
      class(dictionary), intent(in)                 :: dict
      class(nuclearDatabase), pointer, intent(in)   :: ptr
      logical(defBool), optional, intent(in)        :: silent
    end subroutine

    !!
    !! Activate this nuclearDatabase
    !!
    !! Will configure relevant cache(s) for the data
    !!
    !! Args:
    !!   activeMat [in] -> Array of matIdx of materials active in the simulation
    !!
    !! Errors:
    !!   fatalError if activeMat contains materials not defined in the instance
    !!
    subroutine activate(self, activeMat)
      import :: nuclearDatabase, shortInt
      class(nuclearDatabase), intent(inout)       :: self
      integer(shortInt), dimension(:), intent(in) :: activeMat
    end subroutine activate

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
    !! Return pointer to material in a database
    !!
    !! Allows to retrive an access to material data for all databases types
    !!
    !! NOTE: This function can be used to inquire about the presence of matIdx in the database!
    !!
    !! Args:
    !!   matIdx [in] -> material index of required material
    !!
    !! Result:
    !!   Pointer to a material of class materialHandle
    !!
    !! Errors:
    !!   Return null() pointer for invalid material index (not present in database)
    !!
    function getMaterial(self, matIdx) result(mat)
      import :: nuclearDatabase, shortInt, materialHandle
      class(nuclearDatabase), intent(in) :: self
      integer(shortInt), intent(in)      :: matIdx
      class(materialHandle), pointer     :: mat
    end function getMaterial

    !!
    !! Return pointer to nuclide in a database
    !!
    !! Allows to retrive an access to nuclide data for all databases types
    !! If database does not contain nuclides (e.g. MG data) just returns null() pointer
    !!
    !! NOTE: This function can be used to inquire about the presence of nucIdx in the database!
    !!
    !! Args:
    !!   nucIdx [in] -> nuclide index of required material
    !!
    !! Result:
    !!   Pointer to a nuclide of class nuclideHandle
    !!
    !! Errors:
    !!   Return null() pointer for invalid nuclide index (not present in database)
    !!
    function getNuclide(self, nucIdx) result(nuc)
      import :: nuclearDatabase, shortInt, nuclideHandle
      class(nuclearDatabase), intent(in) :: self
      integer(shortInt), intent(in)      :: nucIdx
      class(nuclideHandle), pointer      :: nuc
    end function getNuclide

    !!
    !! Return a pointer to a reaction
    !!
    !! Allows to retrive an access to reaction data for all databases types
    !! Reactions can be associated either with nuclides or materials. Thus, there is ambiguity
    !! whether material or nuclide should be asked to provide reaction data (using matIdx or nuIdx)
    !!
    !! This ambiguity is resolved by following convenction:
    !!   if MT < 0 then reaction is associated with material: idx -> matIdx
    !!   if MT > 0 then reaction is associated with nuclide: idx -> nucIdx
    !!
    !! NOTE: This function can be used to enquire about the presence of data. If the data is
    !!       not present null() pointer is always returned!
    !!
    !! Args:
    !!   MT [in]  -> MT number (identifier) of the requested reaction
    !!   idx [in] -> if MT < 0 matIdx; if MT > 0 nucIdx
    !!
    !! Result:
    !!   Pointer to a reaction data of reactionHandle class
    !!
    !! Error:
    !!   If MT and idx combination is either invalid or not present in the database then null()
    !!   pointer is returned.
    !!
    function getReaction(self, MT, idx) result(reac)
      import :: nuclearDatabase, shortInt, reactionHandle
      class(nuclearDatabase), intent(in) :: self
      integer(shortInt), intent(in)      :: MT
      integer(shortInt), intent(in)      :: idx
      class(reactionHandle),pointer      :: reac
    end function getReaction

    !!
    !! Subroutine that precomputes the majorant cross section to be used by DT
    !!
    !! NOTE: Assumes that the nuclear database has been initialised and activated
    !!
    subroutine initMajorant(self, rand)
      import :: nuclearDatabase, RNG
      class(nuclearDatabase), intent(inout) :: self
      class(RNG), intent(inout)             :: rand
    end subroutine initMajorant

    !!
    !! Return to uninitialised state
    !!
    elemental subroutine kill(self)
      import :: nuclearDatabase
      class(nuclearDatabase), intent(inout) :: self
    end subroutine kill
  end interface

end module nuclearDatabase_inter
