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
  !!        majorant < 1 or 0 >; aceLibrary <nuclear data path> ;} }
  !!
  !! Public Members:
  !!   nuclides    -> array of aceNeutronNuclides with data
  !!   materials   -> array of ceNeutronMaterials with data
  !!   Ebounds     -> array with bottom (1) and top (2) energy bound
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
    real(defReal), dimension(2)                  :: Ebounds   = ZERO
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
    procedure :: updateTotalMatXS
    procedure :: updateMajorantXS
    procedure :: updateMacroXSs
    procedure :: updateTotalNucXS
    procedure :: updateMicroXSs
    procedure :: getScattMicroMajXS

    ! class Procedures
    procedure :: initUrr
    procedure :: initDBRC
    procedure :: initMajorant

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

    self % EBounds = ZERO

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
    integer(shortInt)                     :: idxMT

    ! Catch case of invalid reaction
    !   MT < 0 -> material reaction
    !   MT = 0 -> does not exist
    !   MT = 1 -> N_total has no reaction object
    if ( MT <= 1 ) then
      reac => null()
      return
    end if

    ! Detect invalid indices
    if( idx < 1 .or. idx > size(self % nuclides)) then
      reac => null()
      return
    end if

    ! Get nuclide reaction
    if( MT == N_N_elastic) then
      if (cache_nuclideCache(idx) % needsSabEl) then
        reac => self % nuclides(idx) % thData % elasticOut
      else
        reac => self % nuclides(idx) % elasticScatter
      end if
    else if ( MT == N_fission) then
      reac => self % nuclides(idx) % fission
    else if (MT == N_N_ThermINEL) then
      reac => self % nuclides(idx) % thData % inelasticOut
    else
      ! Find index of MT reaction
      idxMT = self % nuclides(idx) % idxMT % getOrDefault(MT, 0)
      ! See if the MT is present or not
      if(idxMT == 0) then
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
    real(defReal)                         :: E_upper, E_lower, E_min, E_max
    real(defReal)                         :: alpha

    ! Find energy limits to define majorant calculation range
    alpha = 4 / (sqrt( E * A / kT ))
    E_upper = E * (1 + alpha) * (1 + alpha)
    E_lower = E * (1 - alpha) * (1 - alpha)

    ! Find system minimum and maximum energies
    call self % energyBounds(E_min, E_max)

    ! Avoid energy limits being outside system range
    if (E_lower < E_min) E_lower = E_min
    if (E_upper > E_max) E_upper = E_max

    ! Find largest elastic scattering xs in energy range given by E_lower and E_upper
    maj = self % nuclides(nucIdx) % elScatteringMaj(E_lower, E_upper)

  end function getScattMicroMajXS

  !!
  !! Return energy bounds for data in the database
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine energyBounds(self, E_min, E_max)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(out)            :: E_min
    real(defReal), intent(out)            :: E_max

    E_min = self % Ebounds(1)
    E_max = self % Ebounds(2)

  end subroutine energyBounds

  !!
  !! Make sure that totalXS of material with matIdx is at energy E
  !! in ceNeutronChache
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

    associate (mat => cache_materialCache(matIdx))
      ! Set new energy
      mat % E_tot  = E

      ! Clean current total XS
      mat % xss % total = ZERO

      ! Construct total macro XS
      do i = 1, size(self % materials(matIdx) % nuclides)
        dens   = self % materials(matIdx) % dens(i)
        nucIdx = self % materials(matIdx) % nuclides(i)

        ! Update if needed
        if(cache_nuclideCache(nucIdx) % E_tot /= E) then
          call self % updateTotalNucXS(E, nucIdx, rand)
        end if

        ! Add microscopic XSs
        mat % xss % total = mat % xss % total + dens * cache_nuclideCache(nucIdx) % xss % total
      end do

    end associate

  end subroutine updateTotalMatXS

  !!
  !! Make sure that the majorant of ALL Active materials is at energy E
  !! in ceNeutronChache
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

    associate (maj => cache_majorantCache(1) )
      maj % E  = E

      ! Get majorant via the precomputed unionised cross section
      if (self % hasMajorant) then
        idx = binarySearch(self % eGridUnion, E)

        if(idx <= 0) then
          call fatalError(Here,'Failed to find energy: '//numToChar(E)//&
                               ' in unionised majorant grid')
        end if

        associate(E_top => self % eGridUnion(idx + 1), E_low  => self % eGridUnion(idx))
          f = (E - E_low) / (E_top - E_low)
        end associate

        maj % xs = self % majorant(idx+1) * f + (ONE - f) * self % majorant(idx)

      else ! Compute majorant on the fly

        maj % xs = ZERO

        ! Loop over materials
        do i = 1, size(self % activeMat)
          matIdx = self % activeMat(i)

          ! Update if needed
          if( cache_materialCache(matIdx) % E_tot /= E) then
            call self % updateTotalMatXS(E, matIdx, rand)
          end if

          maj % xs = max(maj % xs, cache_materialCache(matIdx) % xss % total)
        end do

      end if

    end associate

  end subroutine updateMajorantXS

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

    associate (mat => cache_materialCache(matIdx))
      ! Set new energy
      mat % E_tot  = E
      mat % E_tail = E

      ! Clean current xss
      call mat % xss % clean()

      ! Construct microscopic XSs
      do i = 1, size(self % materials(matIdx) % nuclides)
        dens   = self % materials(matIdx) % dens(i)
        nucIdx = self % materials(matIdx) % nuclides(i)

        ! Update if needed
        if (cache_nuclideCache(nucIdx) % E_tail /= E .or. cache_nuclideCache(nucIdx) % E_tot /= E) then
          call self % updateMicroXSs(E, nucIdx, rand)
        end if

        ! Add microscopic XSs
        call mat % xss % add(cache_nuclideCache(nucIdx) % xss, dens)
      end do

    end associate

  end subroutine updateMacroXSs

  !!
  !! Make sure that totalXS of nuclide with nucIdx is at energy E
  !! in ceNeutronChache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateTotalNucXS(self, E, nucIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: nucIdx
    class(RNG), optional, intent(inout)   :: rand
    logical(defBool)                      :: needsSab

    associate (nucCache => cache_nuclideCache(nucIdx), &
               nuc      => self % nuclides(nucIdx)     )

      ! Check if the nuclide needs ures probability tables at this energy
      nucCache % needsUrr = (nuc % hasProbTab .and. E >= nuc % urrE(1) .and. E <= nuc % urrE(2))
      ! Check if the nuclide needs S(a,b) at this energy
      nucCache % needsSabEl = (nuc % hasThData .and. E >= nuc % SabEl(1) .and. E <= nuc % SabEl(2))
      nucCache % needsSabInel = (nuc % hasThData .and. E >= nuc % SabInel(1) .and. E <= nuc % SabInel(2))
      needsSab = (nucCache % needsSabEl .or. nucCache % needsSabInel)

      if (nucCache % needsUrr .or. needsSab) then
        call self % updateMicroXSs(E, nucIdx, rand)

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
  subroutine updateMicroXSs(self, E, nucIdx, rand)
    class(aceNeutronDatabase), intent(in) :: self
    real(defReal), intent(in)             :: E
    integer(shortInt), intent(in)         :: nucIdx
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
      if (nucCache % needsUrr) then
        associate(zaidCache => cache_zaidCache(self % nucToZaid(nucIdx)))

          if (zaidCache % E /= E) then
            ! Save random number for temperature correlation
            zaidCache % xi = rand % get()
            zaidCache % E = E
          end if

          call nuc % getUrrXSs(nucCache % xss, nucCache % idx, nucCache % f, E, zaidCache % xi)

        end associate

      elseif (nucCache % needsSabEl .or. nucCache % needsSabInel) then
        call nuc % getThXSs(nucCache % xss, nucCache % idx, nucCache % f, E)

      else
        call nuc % microXSs(nucCache % xss, nucCache % idx, nucCache % f)

      end if

    end associate

  end subroutine updateMicroXSs

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
    type(aceSabCard)                                 :: ACE_Sab
    character(pathLen)                               :: aceLibPath
    character(nameLen)                               :: name, name_file, nucDBRC_temp
    character(:), allocatable                        :: zaid, file
    integer(shortInt)                                :: i, j, envFlag, nucIdx, idx
    integer(shortInt)                                :: maxNuc
    logical(defBool)                                 :: isFissileMat
    integer(shortInt),dimension(:),allocatable       :: nucIdxs, zaidDBRC
    character(nameLen),dimension(:),allocatable      :: nucDBRC
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
        name = trim(mat % nuclides(j) % toChar())
        if (mat % nuclides(j) % hasSab) then
          zaid = trim(mat % nuclides(j) % toChar())
          file = trim(mat % nuclides(j) % file_Sab)
          name = zaid // '+' // file
          deallocate(zaid, file)
        end if
        call nucSet % add(name, IN_SET)
      end do
    end do

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

      idx = index(nucSet % atKey(i),'+')
      if (idx /= 0) then
        name = trim(nucSet % atKey(i))
        name_file = trim(name(idx+1:nameLen))
        name = name(1:idx-1)
      else
        name = nucSet % atKey(i)
      end if

      if(loud) then
        print '(A)', "Building: "// trim(name)// " with index: " //numToChar(nucIdx)
        if (idx /= 0) print '(A)', "including S(alpha,beta) tables with file: " //trim(name_file)
      end if

      call new_neutronACE(ACE, name)
      call self % nuclides(nucIdx) % init(ACE, nucIdx, ptr_ceDatabase)

      ! Initialise S(alpha,beta) tables
      if (idx /= 0 ) then
        call new_moderACE(ACE_Sab, name_file)
        call self % nuclides(nucIdx) % initSab(ACE_Sab)
      end if

      ! Initialise probability tables
      if (self % hasUrr) call self % nuclides(nucIdx) % initUrr(ACE)

      ! Store nucIdx in the dictionary
      call nucSet % atSet(nucIdx, i)
      nucIdx = nucIdx + 1
      i = nucSet % next(i)
    end do

    ! Build Material definitions
    allocate(self % materials(mm_nMat()))
    allocate(nucIdxs(maxNuc))
    do i = 1, mm_nMat()
      mat => mm_getMatPtr(i)

      ! Load nuclide indices on storage space
      isFissileMat = .false.
      do j = 1, size(mat % nuclides)
        name = trim(mat % nuclides(j) % toChar())
        if (mat % nuclides(j) % hasSab) then
          zaid = trim(mat % nuclides(j) % toChar())
          file = trim(mat % nuclides(j) % file_Sab)
          name = zaid // '+' // file
          deallocate(zaid, file)
        end if
        nucIdxs(j) = nucSet % get(name)
        isFissileMat = isFissileMat .or. self % nuclides(nucIdxs(j)) % isFissile()
      end do

      ! Load data into material
      call self % materials(i) % set( matIdx = i, &
                                      database = ptr_ceDatabase, &
                                      fissile = isFissileMat )
      call self % materials(i) % setComposition( mat % dens, nucIdxs(1:size(mat % nuclides)))
    end do

    ! Calculate energy bounds
    self % Ebounds(1) = self % nuclides(1) % eGrid(1)
    j = size(self % nuclides(1) % eGrid)
    self % Ebounds(2) = self % nuclides(1) % eGrid(j)

    do i = 2, size(self % nuclides)
      self % Ebounds(1) = max(self % Ebounds(1), self % nuclides(i) % eGrid(1))
      j = size(self % nuclides(i) % eGrid)
      self % Ebounds(2) = min(self % Ebounds(2), self % nuclides(i) % eGrid(j))
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
      call cache_init(size(self % materials), size(self % nuclides), 1, maxval(self % nucToZaid))
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
    real(defReal)                            :: eRef, eNuc, E, maj, total, dens, urrMaj, &
                                                nucXS, f, eMax, eMin
    logical(defBool)                         :: needsUrr
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
    tmpGrid = self % EBounds(2)

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
        if (eNuc < self % EBounds(1) .or. eNuc > self % EBounds(2)) then
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

    ! Loop over all the energies
    do i = 1, size(self % eGridUnion)

      ! Retrieve current energy
      E = self % eGridUnion(i)

      ! Correct for energies higher or lower than the allowed boundaries
      if (E < self % EBounds(1)) E = self % EBounds(1)
      if (E > self % EBounds(2)) E = self % EBounds(2)

      ! Initialise majorant value for this energy
      maj = ZERO

      ! Loop over active materials
      do j = 1, size(self % activeMat)

        ! Get material index
        matIdx = self % activeMat(j)
        total  = ZERO

        ! Loop over nuclides
        do k = 1, size(self % materials(matIdx) % nuclides)
          dens     = self % materials(matIdx) % dens(k)
          nucIdx   = self % materials(matIdx) % nuclides(k)

          associate (nuc => self % nuclides(nucIdx))

            needsUrr = (nuc % hasProbTab .and. E >= nuc % urrE(1) .and. E <= nuc % urrE(2))

            ! Check if present nuclide uses URR tables
            if (needsUrr) then

              ! Find maximum URR table total XS
              urrIdx = binarySearch(nuc % probTab % eGrid, E)
              urrMaj = nuc % probTab % majorant(urrIdx)

              ! Check if URR tables contain xs or multiplicative factor
              if (nuc % IFF == 1) then
                call nuc % search(eIdx, f, E)
                nucXS = nuc % totalXS(eIdx, f) * urrMaj
              else
                nucXS = urrMaj
              end if

            else
              call self % updateTotalNucXS(E, nucIdx)
              nucXS = cache_nuclideCache(nucIdx) % xss % total

            end if

          end associate

          ! Update total material cross section
          total = total + dens * nucXS

        end do

        maj = max(maj, total)

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
