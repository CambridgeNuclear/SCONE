module aceNeutronDatabase_class

  use numPrecision
  use endfConstants
  use universalVariables
  use errors_mod,         only : fatalError
  use genericProcedures,  only : numToChar, removeDuplicatesSorted, binarySearch
  use dictionary_class,   only : dictionary
  use RNG_class,          only : RNG
  use charMap_class,      only : charMap
  use intMap_class,       only : intMap

  ! Nuclear Data Interfaces
  use nuclearDatabase_inter,        only : nuclearDatabase
  use materialHandle_inter,         only : materialHandle
  use nuclideHandle_inter,          only : nuclideHandle
  use reactionHandle_inter,         only : reactionHandle
  use ceNeutronDatabase_inter,      only : ceNeutronDatabase, ceNeutronDatabase_CptrCast
  use neutronXSPackages_class,      only : neutronMicroXSs
  use ceNeutronMaterial_class,      only : ceNeutronMaterial

  ! Material Menu
  use materialMenu_mod,             only : materialItem, nuclideInfo, mm_nMat => nMat, &
                                           mm_getMatPtr => getMatPtr, mm_nameMap => nameMap

  ! ACE CE Nuclear Data Objects
  use aceLibrary_mod,               only : new_neutronAce, new_moderACE, aceLib_load => load, aceLib_kill => kill
  use aceCard_class,                only : aceCard
  use aceSabCard_class,             only : aceSabCard
  use aceNeutronNuclide_class,      only : aceNeutronNuclide


  ! CE NEUTRON CACHE
  use ceNeutronCache_mod,           only : cache_nuclideCache => nuclideCache, &
                                           cache_materialCache => materialCache, &
                                           cache_majorantCache => majorantCache, &
                                           cache_zaidCache => zaidCache, &
                                           cache_init => init

  ! Scattering procedures
  use scatteringKernels_func,  only : relativeEnergy_constXS, dopplerCorrectionFactor

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: aceNeutronDatabase_TptrCast
  public :: aceNeutronDatabase_CptrCast

  !!
  !! A CE Neutron Database based on ACE file format
  !!
  !! It's possible to use probability tables in the unresolved resonance range if
  !! ures is included in the input file
  !!
  !! Sample input:
  !!   nuclearData {
  !!   handles {
  !!   ce { type aceNeutronDatabase; DBRC (92238 94242); ures < 1 or 0 >;
  !!        majorant < 1 or 0 >; aceLibrary <nuclear data path> ;} 
  !!        #avgDist 3.141;# }
  !!
  !! Public Members:
  !!   nuclides    -> array of aceNeutronNuclides with data
  !!   materials   -> array of ceNeutronMaterials with data
  !!   eBounds     -> array with bottom (1) and top (2) energy bound
  !!   majorant    -> unionised majorant cross section
  !!   eGridUnion  -> unionised energy grid
  !!   activeMat   -> array of materials present in the geometry
  !!   nucToZaid   -> map to link nuclide index to zaid index
  !!   hasUrr      -> ures probability tables flag, it's false by default
  !!   hasDBRC     -> DBRC flag, it's false by default
  !!   hasMajorant -> unionised majorant cross section flag
  !!
  !! Interface:
  !!   nuclearData Interface
  !!   ceNeutronDatabase Interface
  !!
  type, public, extends(ceNeutronDatabase) :: aceNeutronDatabase
    type(aceNeutronNuclide),dimension(:),pointer :: nuclides  => null()
    type(ceNeutronMaterial),dimension(:),pointer :: materials => null()
    real(defReal), dimension(:), allocatable     :: majorant
    real(defReal), dimension(:), allocatable     :: eGridUnion
    real(defReal), dimension(2)                  :: eBounds   = ZERO
    integer(shortInt),dimension(:),allocatable   :: activeMat

    ! Probability tables data
    integer(shortInt),dimension(:),allocatable   :: nucToZaid
    logical(defBool)                             :: hasUrr  = .false.
    logical(defBool)                             :: hasDBRC = .false.
    logical(defBool)                             :: hasMajorant = .false.

  contains

    ! nuclearDatabase Procedures
    procedure :: kill
    procedure :: matNamesMap
    procedure :: getMaterial
    procedure :: getNuclide
    procedure :: getReaction
    procedure :: init
    procedure :: activate

    ! ceNeutronDatabase Procedures
    procedure :: energyBounds
    procedure :: updateMajorantXS
    procedure :: updateTrackMatXS
    procedure :: updateTotalMatXS
    procedure :: updateMacroXSs
    procedure :: updateTotalNucXS
    procedure :: updateMicroXSs
    procedure :: updateTotalTempNucXS
    procedure :: getScattMicroMajXS

    ! class Procedures
    procedure :: initUrr
    procedure :: initDBRC
    procedure :: initMajorant
    procedure :: updateTotalTempMajXS
    procedure :: updateRelEnMacroXSs
    procedure :: makeNuclideName

  end type aceNeutronDatabase


contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(aceNeutronDatabase), intent(inout) :: self

    ! Clean
    if(associated(self % nuclides)) then
      call self % nuclides % kill()
      deallocate(self % nuclides)
    end if

    if(associated(self % materials)) then
      call self % materials % kill()
      deallocate(self % materials)
    end if

    self % eBounds = ZERO

    if(allocated(self % activeMat)) deallocate(self % activeMat)

  end subroutine kill

  !!
  !! Return pointer to material names map
  !!
  !! See nuclearData_inter for  more details
  !!
  function matNamesMap(self) result(map)
    class(aceNeutronDatabase), intent(in) :: self
    type(charMap), pointer                :: map

    map => mm_nameMap

  end function matNamesMap

  !!
  !! Return pointer to material in a database
  !!
  !! See nuclearData_inter for  more details
  !!
  function getMaterial(self, matIdx) result(mat)
    class(aceNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)         :: matIdx
    class(materialHandle), pointer        :: mat

    ! Check bounds and return
    if( 1 <= matIdx .and. matIdx <= size(self % materials)) then
      mat => self % materials(matIdx)
    else
      mat => null()
    end if

  end function getMaterial

  !!
  !! Return pointer to nuclide in a database
  !!
  !! See nuclearData_inter for  more details
  !!
  function getNuclide(self, nucIdx) result(nuc)
    class(aceNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)         :: nucIdx
    class(nuclideHandle), pointer         :: nuc

    ! Check bounds and return
    if( 1 <= nucIdx .and. nucIdx <= size(self % nuclides)) then
      nuc => self % nuclides(nucIdx)
    else
      nuc => null()
    end if

  end function getNuclide

  !!
  !! Return a pointer to a reaction
  !!
  !! See nuclearData_inter for  more details
  !!
  function getReaction(self, MT, idx) result(reac)
    class(aceNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)         :: MT
    integer(shortInt), intent(in)         :: idx
    class(reactionHandle), pointer        :: reac
    integer(shortInt)                     :: idxMT, sabIdx

    ! Catch case of invalid reaction
    !   MT < 0 -> material reaction
    !   MT = 0 -> does not exist
    !   MT = 1 -> N_total has no reaction object
    if (MT <= 1) then
      reac => null()
      return
    end if

    ! Detect invalid indices
    if (idx < 1 .or. idx > size(self % nuclides)) then
      reac => null()
      return
    end if

    ! Get nuclide reaction
    if (MT == N_N_elastic) then
      reac => self % nuclides(idx) % elasticScatter

    else if (MT == N_N_ThermEL) then
      sabIdx  = cache_nuclideCache(idx) % sabIdx
      reac => self % nuclides(idx) % thData(sabIdx) % elasticOut

    else if (MT == N_fission) then
      reac => self % nuclides(idx) % fission

    else if (MT == N_N_ThermINEL) then
      sabIdx  = cache_nuclideCache(idx) % sabIdx
      reac => self % nuclides(idx) % thData(sabIdx) % inelasticOut

    else
      ! Find index of MT reaction
      idxMT = self % nuclides(idx) % idxMT % getOrDefault(MT, 0)
      ! See if the MT is present or not
      if (idxMT == 0) then
        reac => null()
      else
        reac => self % nuclides(idx) % MTdata(idxMT) % kinematics
      end if

    end if

  end function getReaction

  !!
  !! Returns the elastic scattering majorant cross section for a nuclide
  !!
  !! See ceNeutronDatabase for more details
  !!
  function getScattMicroMajXS(self, E, kT, A, nucIdx) result(maj)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    real(defReal), intent(in)             :: kT
    real(defReal), intent(in)             :: A
    integer(shortInt), intent(in)         :: nucIdx
    real(defReal)                         :: maj
    real(defReal)                         :: eUpper, eLower, eMin, eMax
    real(defReal)                         :: alpha

    ! Find energy limits to define majorant calculation range
    alpha = 3.0_defReal * sqrt( kT / (E * A) )
    eUpper = E * (ONE + alpha) * (ONE + alpha)
    eLower = E * (ONE - alpha) * (ONE - alpha)

    ! Find system minimum and maximum energies
    call self % energyBounds(eMin, eMax)

    ! Avoid energy limits being outside system range
    if (eLower < eMin .or. ONE < alpha) eLower = eMin
    if (eUpper > eMax) eUpper = eMax

    ! Find largest elastic scattering xs in energy range given by E_lower and E_upper
    maj = self % nuclides(nucIdx) % getMajXS(eLower, eUpper, N_N_ELASTIC)

  end function getScattMicroMajXS

  !!
  !! Return energy bounds for data in the database
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine energyBounds(self, eMin, eMax)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(out)            :: eMin
    real(defReal), intent(out)            :: eMax

    eMin = self % eBounds(1)
    eMax = self % eBounds(2)

  end subroutine energyBounds

  !!
  !! Make sure that the majorant of ALL active materials is at energy E
  !! in ceNeutronCache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateMajorantXS(self, E, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    class(RNG), optional, intent(inout)   :: rand
    integer(shortInt)                     :: idx, i, matIdx
    real(defReal)                         :: f
    character(100), parameter :: Here = 'updateMajorantXS (aceNeutronDatabase_class.f90)'

    associate (maj => cache_majorantCache(1))
      maj % E  = E

      ! Get majorant via the precomputed unionised cross section
      if (self % hasMajorant) then
        idx = binarySearch(self % eGridUnion, E)

        if (idx <= 0) then
          call fatalError(Here,'Failed to find energy: '//numToChar(E)//&
                               ' in unionised majorant grid')

        end if

        associate (E_top => self % eGridUnion(idx + 1), E_low  => self % eGridUnion(idx))
          f = (E - E_low) / (E_top - E_low)
        end associate

        maj % xs = self % majorant(idx+1) * f + (ONE - f) * self % majorant(idx)

      else ! Compute majorant on the fly

        maj % xs = ZERO

        ! Loop over materials
        do i = 1, size(self % activeMat)
          matIdx = self % activeMat(i)

          ! Update if needed
          if (cache_materialCache(matIdx) % E_track /= E) then
            call self % updateTrackMatXS(E, matIdx, rand)
          end if

          maj % xs = max(maj % xs, cache_materialCache(matIdx) % trackXS)
        end do

      end if

    end associate

  end subroutine updateMajorantXS

  !!
  !! Make sure that trackXS of material with matIdx is at energy E = E_track
  !! in ceNeutronChache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateTrackMatXS(self, E, matIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    class(RNG), optional, intent(inout)   :: rand

    associate (matCache => cache_materialCache(matIdx), &
               mat      => self % materials(matIdx))

      ! Set new energy
      matCache % E_track = E

      if (mat % useTMS(E)) then
        ! The material tracking xs is the temperature majorant in the case of TMS
        call self % updateTotalTempMajXS(E, matIdx)

      else
        ! When TMS is not in use, the material tracking xs is equivalent to the total
        call self % updateTotalMatXS(E, matIdx, rand)
        matCache % trackXS = matCache % xss % total

      end if

    end associate

  end subroutine updateTrackMatXS

  !!
  !! Subroutine to update the temperature majorant in a given material at given temperature
  !!
  !! The function finds the upper and lower limits of the energy range the nuclide majorant
  !! is included into, then adds up the nuclide temperature majorant multiplied by the Doppler
  !! correction factor
  !!
  !! Args:
  !!   E [in]         -> Incident neutron energy for which temperature majorant is found
  !!   matIdx [in]    -> Index of material for which the material temperature majorant is found
  !!
  subroutine updateTotalTempMajXS(self, E, matIdx)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    integer(shortInt)                     :: nucIdx, i
    real(defReal)                         :: dens, corrFact, nucTempMaj

    associate (matCache => cache_materialCache(matIdx), &
               mat      => self % materials(matIdx))

      ! Clean current total XS
      matCache % trackXS = ZERO

      ! loop through all nuclides in material and find sum of majorants
      do i = 1, size(mat % nuclides)

        ! Get nuclide data
        nucIdx  = mat % nuclides(i)
        dens    = mat % dens(i)

        call self % updateTotalTempNucXS(E, mat % kT, nucIdx)

        ! Sum nuclide majorants to find material majorant
        corrFact   = cache_nuclideCache(nucIdx) % doppCorr
        nucTempMaj = cache_nuclideCache(nucIdx) % tempMajXS * corrFact
        matCache % trackXS = matCache % trackXS + dens * nucTempMaj

      end do

    end associate

  end subroutine updateTotalTempMajXS

  !!
  !! Make sure that totalXS of material with matIdx is at energy E
  !! in ceNeutronCache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateTotalMatXS(self, E, matIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    class(RNG), optional, intent(inout)   :: rand
    integer(shortInt)                     :: i, nucIdx
    real(defReal)                         :: dens

    associate (matCache => cache_materialCache(matIdx), &
               mat      => self % materials(matIdx))

      ! Set new energy and clean current total XS
      matCache % E_tot = E
      matCache % xss % total = ZERO

      if (mat % useTMS(E)) then
        ! When TMS is in use, the total xs is retrieved sampling the nuclides' relative
        ! energies given the temperature difference between material temperature and
        ! temperature of the nuclides' base cross sections
        call self % updateRelEnMacroXSs(E, matIdx, rand)

      else
        ! Construct total macro XS
        do i = 1, size(mat % nuclides)
          dens   = mat % dens(i)
          nucIdx = mat % nuclides(i)

          ! Update if needed
          if (cache_nuclideCache(nucIdx) % E_tot /= E) then
            call self % updateTotalNucXS(E, nucIdx, mat % kT, rand)
          end if

          ! Add microscopic XSs
          matCache % xss % total = matCache % xss % total + &
                                   dens * cache_nuclideCache(nucIdx) % xss % total
        end do

      end if

    end associate

  end subroutine updateTotalMatXS

  !!
  !! Make sure that the macroscopic XSs for the material with matIdx are set
  !! to energy E in ceNeutronCache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateMacroXSs(self, E, matIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    class(RNG), optional, intent(inout)   :: rand
    integer(shortInt)                     :: i, nucIdx
    real(defReal)                         :: dens

    associate(mat      => self % materials(matIdx), &
              matCache => cache_materialCache(matIdx))

      ! Clean current xss
      call matCache % xss % clean()

      if (mat % useTMS(E)) then
        ! When TMS is in use, the xss are retrieved sampling the nuclides' relative
        ! energies given the temperature difference between material temperature and
        ! temperature of the nuclides' base cross sections
        call self % updateRelEnMacroXSs(E, matIdx, rand)

      else

        ! Set new energy
        matCache % E_tot  = E
        matCache % E_tail = E

        ! Construct microscopic XSs
        do i = 1, size(mat % nuclides)
          dens   = mat % dens(i)
          nucIdx = mat % nuclides(i)

          ! Update if needed
          if (cache_nuclideCache(nucIdx) % E_tail /= E .or. cache_nuclideCache(nucIdx) % E_tot /= E) then
            call self % updateMicroXSs(E, nucIdx, mat % kT, rand)
          end if

          ! Add microscopic XSs
          call matCache % xss % add(cache_nuclideCache(nucIdx) % xss, dens)
        end do

      end if

    end associate

  end subroutine updateMacroXSs

  !!
  !! Subroutine to update the macroscopic cross sections in a given material
  !! at given temperature, sampling a relative energy per each nuclide and applying
  !! a Doppler correction factor
  !!
  !! Args:
  !!   E [in]         -> Incident neutron energy for which the relative energy xss are found
  !!   matIdx [in]    -> Index of material for which the relative energy xss are found
  !!
  subroutine updateRelEnMacroXSs(self, E, matIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: matIdx
    class(RNG), optional, intent(inout)   :: rand
    integer(shortInt)                     :: i, nucIdx
    real(defReal)                         :: dens, nuckT, A, deltakT, eRel, eMin, &
                                             eMax, doppCorr
    character(100), parameter :: Here = 'updateRelEnMacroXSs (aceNeutronDatabase_class.f90)'

    associate(mat      => self % materials(matIdx), &
              matCache => cache_materialCache(matIdx))

      ! Check if relative energy cross sections have been retrieved before
      if (E /= matCache % E_rel) then

        ! Clean current xss
        call matCache % xssRel % clean()
        matCache % E_rel = E

        ! Construct microscopic XSs
        do i = 1, size(mat % nuclides)

          dens   = mat % dens(i)
          nucIdx = mat % nuclides(i)
          nuckT  = self % nuclides(nucIdx) % getkT()
          A      = self % nuclides(nucIdx) % getMass()
          deltakT = mat % kT - nuckT

          eRel = relativeEnergy_constXS(E, A, deltakT, rand)

          ! Call through system minimum and maximum energies
          call self % energyBounds(eMin, eMax)

          ! avoid sampled relative energy from MB dist extending into energies outside system range
          if (eRel < eMin) eRel = eMin
          if (eMax < eRel) eRel = eMax

          associate(nucCache => cache_nuclideCache(nucIdx))

            ! Doppler correction factor for low energies
            doppCorr = dopplerCorrectionFactor(E, A, deltakT)

            ! Update if needed
            if (nucCache % E_tail /= eRel .or. nucCache % E_tot /= eRel) then
              call self % updateMicroXSs(eRel, nucIdx, mat % kT, rand)
            end if

            ! Add microscopic XSs
            call matCache % xssRel % add(nucCache % xss, dens * doppCorr)

          end associate

        end do

      end if

      ! Update cache, and ensure that the energy indicators are reset to avoid wrong look-ups
      matCache % xss = matCache % xssRel
      matCache % E_tot  = ZERO
      matCache % E_tail = ZERO

    end associate

  end subroutine updateRelEnMacroXSs

  !!
  !! Make sure that totalXS of nuclide with nucIdx is at energy E
  !! in ceNeutronCache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateTotalNucXS(self, E, nucIdx, kT, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: nucIdx
    real(defReal), intent(in)             :: kT
    class(RNG), optional, intent(inout)   :: rand

    associate (nucCache => cache_nuclideCache(nucIdx), &
               nuc      => self % nuclides(nucIdx)     )

      ! Check if the nuclide needs ures probability tables or S(a,b) at this energy
      if (nuc % needsUrr(E) .or. nuc % needsSabEl(E) .or. nuc % needsSabInel(E)) then
        call self % updateMicroXSs(E, nucIdx, kT, rand)

      else
        nucCache % E_tot  = E
        call nuc % search(nucCache % idx, nucCache % f, E)
        nucCache % xss % total = nuc % totalXS(nucCache % idx, nucCache % f)

      end if

    end associate

  end subroutine updateTotalNucXS

  !!
  !! Make sure that the microscopic XSs for the nuclide with nucIdx are set
  !! to energy E in ceNeutronCache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateMicroXSs(self, E, nucIdx, kT, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: nucIdx
    real(defReal), intent(in)             :: kT
    class(RNG), optional, intent(inout)   :: rand

    associate (nucCache => cache_nuclideCache(nucIdx), &
               nuc      => self % nuclides(nucIdx)     )

      nucCache % E_tail = E

      ! In case the total XS hasn't been retrieved before (during tracking)
      if (nucCache % E_tot /= E) then
        nucCache % E_tot  = E
        call nuc % search(nucCache % idx, nucCache % f, E)
      end if

      ! Overwrites all the micro cross sections in cache
      ! Check if probability tables should be read
      if (nuc % needsUrr(E)) then
        associate(zaidCache => cache_zaidCache(self % nucToZaid(nucIdx)))

          if (zaidCache % E /= E) then
            ! Save random number for temperature correlation
            zaidCache % xi = rand % get()
            zaidCache % E = E
          end if

          call nuc % getUrrXSs(nucCache % xss, nucCache % idx, nucCache % f, E, zaidCache % xi)

        end associate

      ! Check if S(a,b) should be read
      elseif (nuc % needsSabEl(E) .or. nuc % needsSabInel(E)) then
        call nuc % getThXSs(nucCache % xss, nucCache % idx, nucCache % f, E, kT, rand)

      else
        call nuc % microXSs(nucCache % xss, nucCache % idx, nucCache % f)

      end if

    end associate

  end subroutine updateMicroXSs

  !!
  !! Subroutine to retrieve the nuclide total majorant cross section
  !! over a calculated energy range (needed for TMS) and update the nuclide
  !! cache with majorant value, energy, deltakT and Doppler correction factor
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateTotalTempNucXS(self, E, kT, nucIdx)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    real(defReal), intent(in)             :: kT
    integer(shortInt), intent(in)         :: nucIdx
    real(defReal)                         :: eUpper, eLower, eMin, eMax, nuckT, &
                                             alpha, deltakT, A
    character(100), parameter :: Here = 'updateTotalTempNucXS (aceNeutronDatabase_class.f90)'

    associate (nuc => self % nuclides(nucIdx) , &
               nucCache => cache_nuclideCache(nucIdx))

      nuckT   = nuc % getkT()
      A       = nuc % getMass()
      deltakT = kT - nuckT

      ! Check if an update is required
      if (nucCache % E_maj /= E .or. nucCache % deltakT /= deltakT) then

        ! Find energy limits to define majorant calculation range
        alpha = 3.0_defReal * sqrt( deltakT / (E * A) )
        eUpper = E * (ONE + alpha) * (ONE + alpha)
        eLower = E * (ONE - alpha) * (ONE - alpha)

        ! Find system minimum and maximum energies
        call self % energyBounds(eMin, eMax)

        ! Avoid energy limits being outside system range
        if (eLower < eMin .or. ONE < alpha) eLower = eMin
        if (eUpper > eMax) eUpper = eMax

        ! Doppler g correction factor for low energies
        nucCache % doppCorr = dopplerCorrectionFactor(E, A, deltakT)

        ! Get nuclide total majorant cross section in the energy range
        nucCache % tempMajXS = nuc % getMajXS(eLower, eUpper, N_TOTAL)

        ! Save additional info
        nucCache % deltakT = deltakT
        nucCache % E_maj   = E

      end if

    end associate

  end subroutine updateTotalTempNucXS

  !!
  !! Initialise Database from dictionary and pointer to self
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine init(self, dict, ptr, silent )
    class(aceNeutronDatabase), target, intent(inout) :: self
    class(dictionary), intent(in)                    :: dict
    class(nuclearDatabase), pointer, intent(in)      :: ptr
    logical(defBool), optional, intent(in)           :: silent
    logical(defBool)                                 :: loud
    type(materialItem), pointer                      :: mat
    class(ceNeutronDatabase), pointer                :: ptr_ceDatabase
    type(charMap)                                    :: nucSet
    type(aceCard)                                    :: ACE
    type(aceSabCard)                                 :: ACE_Sab1, ACE_Sab2
    character(pathLen)                               :: aceLibPath
    character(nameLen)                               :: name, name_file1, name_file2, nucDBRC_temp
    integer(shortInt)                                :: i, j, envFlag, nucIdx, idx, idx1, idx2
    integer(shortInt)                                :: maxNuc
    logical(defBool)                                 :: isFissileMat
    integer(shortInt),dimension(:),allocatable       :: nucIdxs, zaidDBRC
    character(nameLen),dimension(:),allocatable      :: nucDBRC
    real(defReal)                                    :: A, nuckT, eUpSab, eUpSabNuc, &
                                                        eLowURR, eLowUrrNuc, alpha, &
                                                        deltakT, eUpper, eLower, kT, &
                                                        temp
    real(defReal), dimension(2)                      :: sabT
    integer(shortInt), parameter :: IN_SET = 1, NOT_PRESENT = 0
    character(100), parameter :: Here = 'init (aceNeutronDatabase_class.f90)'

    ! Set build console output flag
    if(present(silent)) then
      loud = .not.silent
    else
      loud = .true.
    end if

    ! Verify pointer
    if (.not.associated(ptr, self)) then
      call fatalError(Here,"Pointer needs to be associated with the self")
    end if

    ! Cast pointer to ceNeutronDatabase
    ptr_ceDatabase => ceNeutronDatabase_CptrCast(ptr)
    if(.not.associated(ptr_ceDatabase)) call fatalError(Here,"Should not happen. WTF?!")

    ! Create list of all nuclides. Loop over materials
    ! Find maximum number of nuclides: maxNuc
    do i = 1, mm_nMat()
      mat => mm_getMatPtr(i)
      maxNuc = max(maxNuc, size(mat % nuclides))

      ! Add all nuclides in material to the map
      do j = 1, size(mat % nuclides)
        name = self % makeNuclideName(mat % nuclides(j))
        call nucSet % add(name, IN_SET)
      end do
    end do

    ! Check for a minimum average collision distance
    if (dict % isPresent('avgDist')) then
      call dict % get(temp, 'avgDist')

      if (temp <= ZERO) then
        call fatalError(Here, 'Must have a finite, positive minimum average collision distance')
      end if

      self % collisionXS = ONE / temp

    end if

    ! Get path to ACE library
    call dict % get(aceLibPath,'aceLibrary')

    ! Check if probability tables are on in the input file
    call dict % getOrDefault(self % hasUrr, 'ures', .false.)

    if(aceLibPath == '$SCONE_ACE') then
      ! Get Path from enviromental variable
      call get_environment_variable("SCONE_ACE", aceLibPath, status = envFlag)

      ! Process potential errors
      if(envFlag == -1) then
        call fatalError(Here,'$SCONE_ACE EnVar must have length smaller then: '//numToChar(pathLen))

      else if(envFlag == 1) then
        call fatalError(Here,"EnVar $SCONE_ACE does not exist! Need to point to ACE Library")

      else if(envFlag == 2) then
        call fatalError(Here,"Compiler does not support EnVariables. &
                              &Replace $SCONE_ACE with path in input file!")
      else if(envFlag /= 0) then
        call fatalError(Here,"Impossible value of envFlag:"//numToChar(envFlag))

      end if
    end if

    ! Load library
    call aceLib_load(aceLibPath)

    ! Check if DBRC is listed in the input file
    if (dict % isPresent('DBRC')) then

      ! Set flag to true if DBRC nucs are in input file
      self % hasDBRC = .true.

      ! Call through list of DBRC nuclides
      call dict % get(zaidDBRC, 'DBRC')
      allocate(nucDBRC(size(zaidDBRC)))

      ! Add all DBRC nuclides to nucSet for initialisation
      do i = 1, size(zaidDBRC)
        nucDBRC(i)  = numToChar(zaidDBRC(i))
        nucDBRC_temp = trim(nucDBRC(i))//'.00'
        call nucSet % add(nucDBRC_temp, IN_SET)
      end do

    end if

    ! Build nuclide definitions
    allocate(self % nuclides(nucSet % length()))
    i = nucSet % begin()
    nucIdx = 1
    do while (i /= nucSet % end())

      idx1 = index(nucSet % atKey(i),'+')
      idx2 = index(nucSet % atKey(i),'#')
      if (idx1 /= 0) then
        name = trim(nucSet % atKey(i))
        if (idx2 == 0) then
          name_file1 = trim(name(idx1+1:nameLen))
        else
          name_file1 = trim(name(idx1+1:idx2-1))
          name_file2 = trim(name(idx2+1:nameLen))
        end if
        name = name(1:idx1-1)
      else
        name = nucSet % atKey(i)
      end if

      if(loud) then
        print '(A)', "Building: "// trim(name)// " with index: " //numToChar(nucIdx)
        if (idx1 /= 0 .and. idx2 == 0) &
                print '(A)', "including S(alpha,beta) table with file: " //trim(name_file1)
        if (idx1 /= 0 .and. idx2 /= 0) &
                print '(A)', "including S(alpha,beta) tables with files: " //trim(name_file1)//' '//trim(name_file2)
      end if

      call new_neutronACE(ACE, name)
      call self % nuclides(nucIdx) % init(ACE, nucIdx, ptr_ceDatabase)

      ! Initialise S(alpha,beta) tables
      if (idx1 /= 0 ) then
        call new_moderACE(ACE_Sab1, name_file1)
        if (idx2 /= 0) then
          call new_moderACE(ACE_Sab2, name_file2)
          call self % nuclides(nucIdx) % initSab(ACE_Sab1, ACE_Sab2)
        else
          call self % nuclides(nucIdx) % initSab(ACE_Sab1)
        end if
      end if

      ! Initialise probability tables
      if (self % hasUrr) call self % nuclides(nucIdx) % initUrr(ACE)

      ! Store nucIdx in the dictionary
      call nucSet % atSet(nucIdx, i)
      nucIdx = nucIdx + 1
      i = nucSet % next(i)
    end do

    ! Calculate energy bounds
    self % eBounds(1) = self % nuclides(1) % eGrid(1)
    j = size(self % nuclides(1) % eGrid)
    self % eBounds(2) = self % nuclides(1) % eGrid(j)

    do i = 2, size(self % nuclides)
      self % eBounds(1) = max(self % eBounds(1), self % nuclides(i) % eGrid(1))
      j = size(self % nuclides(i) % eGrid)
      self % eBounds(2) = min(self % eBounds(2), self % nuclides(i) % eGrid(j))
    end do

    ! Build Material definitions
    allocate(self % materials(mm_nMat()))
    allocate(nucIdxs(maxNuc))
    do i = 1, mm_nMat()
      mat => mm_getMatPtr(i)

      ! Load nuclide indices on storage space
      ! Find if material is fissile and if stochastic
      ! mixing temperature bounds are respected
      isFissileMat = .false.
      ! Loop over nuclides
      do j = 1, size(mat % nuclides)
        name = self % makeNuclideName(mat % nuclides(j))
        
        ! Find nuclide definition to see if fissile
        ! Also used for checking stochastic mixing bounds
        nucIdxs(j) = nucSet % get(name)
        isFissileMat = isFissileMat .or. self % nuclides(nucIdxs(j)) % isFissile()
          
        ! Check to ensure stochastic mixing temperature 
        ! is bounded by Sab temperatures
        if (mat % nuclides(j) % sabMix) then
          sabT = self % nuclides(nucIdxs(j)) % getSabTBounds()
          kT = mat % T * kBoltzmann / joulesPerMeV
          if ((kT < sabT(1)) .or. (kT > sabT(2))) call fatalError(Here,&
                'Material temperature must be bounded by the provided S(alpha,beta) data. '//&
                'The material temperature is '//numToChar(kT * joulesPerMeV / kBoltzmann)//&
                'K while the data bounds are '//numToChar(sabT(1) * joulesPerMeV / kBoltzmann)//&
                'K and '//numToChar(sabT(2) * joulesPerMeV / kBoltzmann)//'K.')
        end if

      end do

      ! Load data into material
      call self % materials(i) % set( name     = mat % name,     &
                                      matIdx   = i,              &
                                      database = ptr_ceDatabase, &
                                      temp     = mat % T,        &
                                      hasTMS   = mat % hasTMS,   &
                                      fissile  = isFissileMat )
      call self % materials(i) % setComposition( mat % dens, nucIdxs(1:size(mat % nuclides)))

      eUpSab  = self % eBounds(1)
      eLowURR = self % eBounds(2)

      if (mat % hasTMS) then

        ! Loop again to find energy limits of S(a,b) and URES for TMS applicability
        do j = 1, size(mat % nuclides)

          ! Find nuclide information
          idx     = nucIdxs(j)
          nuckT   = self % nuclides(idx) % getkT()
          A       = self % nuclides(idx) % getMass()
          deltakT = self % materials(i) % kT - nuckT

          ! Call fatal error if material temperature is lower then base nuclide temperature
          if (deltakT < ZERO) then
            call fatalError(Here, "Material temperature must be greater than the nuclear data temperature.")
          end if

          ! Find nuclide upper S(a,b) energy
          eUpSabNuc = max(self % nuclides(idx) % SabEl(2), self % nuclides(idx) % SabInel(2))

          ! Find energy limits to define majorant calculation range
          if (eUpSabNuc > ZERO) then
            alpha = 4.0_defReal * sqrt( deltakT / (eUpSabNuc * A) )
            eUpper = eUpSabNuc * (ONE + alpha) * (ONE + alpha)
          else
            eUpper = ZERO
          end if

          eLowUrrNuc = self % nuclides(idx) % urrE(1)

          if (eLowUrrNuc /= ZERO) then
            alpha = 4.0_defReal * sqrt( deltakT / (eLowUrrNuc * A) )
            eLower = eLowUrrNuc * (ONE - alpha) * (ONE - alpha)
          else
            eLower = self % eBounds(2)
          end if

          eUpSab  = max(eUpSab, eUpper)
          eLowURR = min(eLowURR, eLower)

        end do

      end if

      ! Load data into material
      call self % materials(i) % set(eUpperSab = eUpSab, eLowerURR = eLowURR)

    end do

    ! Read unionised majorant flag
    call dict % getOrDefault(self % hasMajorant, 'majorant', .true.)

    ! If on, initialise probability tables for ures
    if (self % hasUrr) then
       call self % initUrr()
    end if

    ! If on, initialise DBRC
    if (self % hasDBRC) then
      call self % initDBRC(nucDBRC, nucSet, self % mapDBRCnuc)
    end if

    !! Clean up
    call aceLib_kill()

  end subroutine init

  !!
  !! Makes a nuclide's name 
  !! Uniquely identifies nuclides with S(alpha,beta) data
  !! variants, including stochastic mixing
  !!
  function makeNuclideName(self, nuclide) result(name)
    class(aceNeutronDatabase), intent(in) :: self
    type(nuclideInfo), intent(in)         :: nuclide
    character(nameLen)                    :: name
    character(:), allocatable             :: file
        
    name = trim(nuclide % toChar())

    ! Name is extended if there is S(alpha,beta) to 
    ! uniquely identify from data without thermal
    ! scattering
    if (nuclide % hasSab) then
 
      file = trim(nuclide % file_Sab1)
      name = trim(name) // '+' // file
      deallocate(file)
     
      ! Attach second Sab file for stochastic mixing
      if (nuclide % sabMix) then
        file = trim(nuclide % file_Sab2)
        name = trim(name) // '#' // file
        deallocate(file)
      end if
 
    end if

  end function makeNuclideName

  !!
  !!  Create list of nuclides with same ZAID, but possibly different temperatures
  !!
  !!  NOTE: compares the first 5 letters of the ZAID.TT. It would be wrong with isotopes
  !!        with Z > 99
  !!
  subroutine initUrr(self)
    class(aceNeutronDatabase), intent(inout) :: self
    integer(shortInt)                        :: i, j
    character(nameLen)                       :: zaid
    type(charMap)                            :: map
    integer(shortInt), parameter :: BEGIN=1, NOT_PRESENT = -12

    ! Allocate array to map ZAIDs
    allocate(self % nucToZaid(size(self % nuclides)))

    ! Initialise ZAID map
    write(zaid,'(A5)') adjustl(self % nuclides(1) % ZAID)

    call map % add(zaid, BEGIN)
    self % nucToZaid(1) = 1

    j = 1
    ! Loop over all nuclides
    do i = 2,size(self % nuclides)
      ! Get the ZAID without temperature -> only compares the first 5 letters of the ZAID.TT string
      write(zaid,'(A5)') adjustl(self % nuclides(i) % ZAID)
      ! Create ZAID list and mapping
      if (map % getOrDefault(zaid, NOT_PRESENT) == NOT_PRESENT) then
        j = j+1
        call map % add(zaid,j)
        self % nucToZaid(i) = j
      else
        self % nucToZaid(i) = map % get(zaid)
      end if
    end do

  end subroutine initUrr

  !!
  !!  Checks through all nuclides, creates map with nuclides present and corresponding 0K nuclide
  !!
  subroutine initDBRC(self, nucDBRC, nucSet, map)
    class(aceNeutronDatabase), intent(inout)     :: self
    character(nameLen), dimension(:), intent(in) :: nucDBRC
    type(charMap), intent(in)                    :: nucSet
    type(intMap), intent(out)                    :: map
    integer(shortInt)                            :: i, j, idx0K, last
    character(nameLen)                           :: nuc0K, nucTemp

    idx0K = 1

    ! Loop through DBRC nuclides
    do i = 1, size(nucDBRC)

      ! Get ZAID with 0K temperature code
      nuc0K = trim(nucDBRC(i))//'.00'

      ! Find the nucIdxs of the 0K DBRC nuclides
      idx0K = nucSet % get(nuc0K)

      ! Loop through nucSet to find the nucIdxs of the DBRC nuclides with
      ! temperature different from 0K
      j = nucSet % begin()
      do while (j /= nucSet % end())

        ! Get ZAIDs in the nuclide set without temperature code
        nucTemp = nucSet % atKey(j)
        ! Remove temperature code (.TT)
        last = len_trim(nucTemp)
        nucTemp = nucTemp(1:last-3)

        ! If the ZAID of the DBRC nuclide matches one in nucSet, save them in the map
        if (nucTemp == nucDBRC(i)) then
          call map % add(nucSet % atVal(j), idx0K)
          ! Set nuclide DBRC flag on
          call self % nuclides(nucSet % atVal(j)) % set(dbrc=.true.)
        end if

        ! Increment index
        j = nucSet % next(j)

      end do

    end do

  end subroutine initDBRC

  !!
  !! Activate this nuclearDatabase
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine activate(self, activeMat, silent)
    class(aceNeutronDatabase), intent(inout)    :: self
    integer(shortInt), dimension(:), intent(in) :: activeMat
    logical(defBool), optional, intent(in)      :: silent
    logical(defBool)                            :: loud

    ! Load active materials
    if(allocated(self % activeMat)) deallocate(self % activeMat)
    self % activeMat = activeMat

    ! Configure Cache
    if (self % hasUrr) then
      call cache_init(size(self % materials), size(self % nuclides), nZaid = maxval(self % nucToZaid))
    else
      call cache_init(size(self % materials), size(self % nuclides))
    end if

    ! If unionised majorant cross section is requested, build it
    if (self % hasMajorant) then

      ! Set build console output flag
      if (present(silent)) then
        loud = .not. silent
      else
        loud = .true.
      end if

      ! Precompute majorant cross section
      call self % initMajorant(loud)

    end if

  end subroutine activate

  !!
  !! Precomputes majorant cross section
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine initMajorant(self, loud)
    class(aceNeutronDatabase), intent(inout) :: self
    logical(defBool), intent(in)             :: loud
    real(defReal), dimension(:), allocatable :: tmpGrid
    integer(shortInt)                        :: i, j, k, matIdx, nNuc, nucIdx, isDone, &
                                                sizeGrid, eIdx, nucIdxLast, eIdxLast, &
                                                urrIdx
    type(intMap)                             :: nucSet
    real(defReal)                            :: eRef, eNuc, E, maj, trackXS, dens, urrMaj, &
                                                nucXS, f, eMax, eMin
    class(RNG), allocatable                  :: rand
    integer(shortInt), parameter :: IN_SET = 1, NOT_PRESENT = 0
    real(defReal), parameter     :: NUDGE = 1.0e-06_defReal

    ! Find the size of the unionised energy grid (with duplicates)
    ! Initialise size
    sizeGrid = 0

    ! Loop over active materials
    do i = 1, size(self % activeMat)

      ! Get current material index and number of nuclides in that material
      matIdx = self % activeMat(i)
      nNuc = size(self % materials(matIdx) % nuclides)

      ! Loop over nuclides present in that material
      do j = 1, nNuc

        ! Get index and check if it's already been added to the set
        nucIdx = self % materials(matIdx) % nuclides(j)
        isDone = nucSet % getOrDefault(nucIdx, NOT_PRESENT)

        ! If it's a new nuclide, add it to the set and find the size of its energy grid
        if (isDone /= IN_SET) then

          ! Add nuclide to the set
          call nucSet % add(nucIdx, IN_SET)

          ! Update energy grid size
          sizeGrid = sizeGrid + size(self % nuclides(nucIdx) % eGrid)

          ! If URR probability tables or S(a,b) tables are used, add their energy
          ! boundary values to the grid to minimise interpolation errors
          if (self % nuclides(nucIdx) % hasProbTab) sizeGrid = sizeGrid + 2
          if (self % nuclides(nucIdx) % hasThData)  sizeGrid = sizeGrid + 3

        end if

      end do
    end do

    ! Allocate temporary grid vector and initialise to largest value allowed
    allocate(tmpGrid(sizeGrid))
    tmpGrid = self % eBounds(2)

    ! Loop over the energy grid
    i = 1
    do while (i < sizeGrid)

      ! Loop over all nuclides in the set - here the value of the intMap is used as an energy index
      j = nucSet % begin()
      do while (j /= nucSet % end())

        ! Retrieve energy in the grid and nuclide information
        eRef    = tmpGrid(i)
        nucIdx  = nucSet % atKey(j)
        eIdx    = nucSet % atVal(j)

        ! Check if we already added all the energy values for this nuclide
        if (eIdx > size(self % nuclides(nucIdx) % eGrid)) then
          j = nucSet % next(j)
          cycle
        end if

        ! Get energy from nuclide grid
        eNuc = self % nuclides(nucIdx) % eGrid(eIdx)

        ! Check if the energy from the nuclide grid is out of bounds
        if (eNuc < self % eBounds(1) .or. eNuc > self % eBounds(2)) then
          j = nucSet % next(j)
          cycle
        end if

        ! Add energy value in the sorted grid, and save index of current nuclide
        if (eNuc <= eRef) then
          tmpGrid(i) = eNuc
          nucIdxLast = nucIdx
          eIdxLast   = eIdx
        end if

        j = nucSet % next(j)

      end do

      ! Increment the energy index saved in the intMap for the nuclides whose energy was added
      call nucSet % add(nucIdxLast, eIdxLast + 1)

      ! Loop over all nuclides again to add S(a,b) and ures energy boundaries to grid
      j = nucSet % begin()
      do while (j /= nucSet % end())

        ! Retrieve energy in the grid and nuclide information
        if (i /= 1) then
          eMin = tmpGrid(i - 1)
        else
          eMin = ZERO
        end if
        eMax = tmpGrid(i)
        nucIdx  = nucSet % atKey(j)

        ! Check for URR probability tables
        if (self % nuclides(nucIdx) % hasProbTab) then

          ! Lower energy boundary
          E = self % nuclides(nucIdx) % urrE(1)
          if (E >= eMin .and. E < eMax) then
            tmpGrid(i) = E
            tmpGrid(i + 1) = eMax
            ! Update counter
            i = i + 1
          end if

          ! Upper energy boundary
          E = self % nuclides(nucIdx) % urrE(2)
          if (E >= eMin .and. E < eMax) then
            tmpGrid(i) = E
            tmpGrid(i + 1) = eMax
            ! Update counter
            i = i + 1
          end if

        end if

        ! Check for Sab tables
        if (self % nuclides(nucIdx) % hasThData) then

          ! Elastic upper energy boundary (NOTE: lower boundary is fixed)
          E = self % nuclides(nucIdx) % SabEl(2)
          if (E >= eMin .and. E < eMax ) then
            tmpGrid(i) = E
            tmpGrid(i + 1) = eMax
            ! Update counter
            i = i + 1
          end if

          ! Inelastic lower energy boundary
          E = self % nuclides(nucIdx) % SabInel(1)
          if (E >= eMin .and. E < eMax ) then
            tmpGrid(i) = E
            tmpGrid(i + 1) = eMax
            ! Update counter
            i = i + 1
          end if

          ! Inelastic upper energy boundary
          E = self % nuclides(nucIdx) % SabInel(2)
          if (E >= eMin .and. E < eMax ) then
            tmpGrid(i) = E
            tmpGrid(i + 1) = eMax
            ! Update counter
            i = i + 1
          end if

        end if

        j = nucSet % next(j)

      end do

      i = i + 1

    end do

    ! Save final grid and remove duplicates
    self % eGridUnion = removeDuplicatesSorted(tmpGrid)

    if (loud) then
      print '(A)', 'CE unionised energy grid has size: '//numToChar(size(self % eGridUnion))
    end if

    ! Allocate unionised majorant
    allocate(self % majorant(size(self % eGridUnion)))

    ! Initialise RNG needed to call update XS routines. The initial seed doesn't
    ! matter because the RNG is only used to sample from probability tables, which
    ! are corrected for the majorant anyway
    allocate(rand)
    call rand % init(1_longInt)

    ! Loop over all the energies
    do i = 1, size(self % eGridUnion)

      ! Retrieve current energy
      E = self % eGridUnion(i)

      ! Correct for energies higher or lower than the allowed boundaries
      if (E < self % eBounds(1)) E = self % eBounds(1)
      if (E > self % eBounds(2)) E = self % eBounds(2)

      ! Initialise majorant value for this energy
      maj = ZERO

      ! Loop over active materials
      do j = 1, size(self % activeMat)

        ! Get material index
        matIdx = self % activeMat(j)

        ! Get material tracking cross section
        call self % updateTrackMatXS(E, matIdx, rand)
        trackXS = cache_materialCache(matIdx) % trackXS

        ! Loop over nuclides to check and correct for ures
        do k = 1, size(self % materials(matIdx) % nuclides)
          dens     = self % materials(matIdx) % dens(k)
          nucIdx   = self % materials(matIdx) % nuclides(k)

          associate (nuc => self % nuclides(nucIdx))

            ! Check if present nuclide uses URR tables
            if (nuc % needsUrr(E)) then

              ! Find maximum URR table total XS
              urrIdx = binarySearch(nuc % probTab % eGrid, E)
              urrMaj = nuc % probTab % majorant(urrIdx)

              ! Check if URR tables contain xs or multiplicative factor
              if (nuc % IFF == 1) then
                call nuc % search(eIdx, f, E)
                nucXS  = nuc % totalXS(eIdx, f) * urrMaj
              else
                nucXS = urrMaj
              end if

            ! Update total material cross section
            trackXS = trackXS + dens * (nucXS - cache_nuclideCache(nucIdx) % xss % total)

            end if

          end associate

        end do

        ! Select majorant cross section
        maj = max(maj, trackXS)

      end do

      ! Save majorant for this energy. Nudge it up to avoid small discrepancies
      self % majorant(i) = maj * (ONE + NUDGE)

    end do

    if (loud) print '(A)', 'CE unionised majorant cross section calculation completed'

  end subroutine initMajorant

  !!
  !! Cast nuclearDatabase pointer to aceNeutronDatabase type pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of aceNeutronDatabase type
  !!   Target points to source if source is aceNeutronDatabase type
  !!
  pure function aceNeutronDatabase_TptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    type(aceNeutronDatabase), pointer           :: ptr

    select type(source)
      type is(aceNeutronDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function aceNeutronDatabase_TptrCast

  !!
  !! Cast nuclearDatabase pointer to aceNeutronDatabase class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of aceNeutronDatabase class
  !!   Target points to source if source is aceNeutronDatabase class
  !!
  pure function aceNeutronDatabase_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(aceNeutronDatabase), pointer          :: ptr

    select type(source)
      class is(aceNeutronDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function aceNeutronDatabase_CptrCast



end module aceNeutronDatabase_class
