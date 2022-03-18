module aceNeutronDatabase_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures, only : fatalError, numToChar
  use dictionary_class,  only : dictionary
  use RNG_class,         only : RNG
  use charMap_class,     only : charMap
  use intMap_class,      only : intMap

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
  !!   ce {type aceNeutronDatabase; ures <1 or 0>; aceLibrary <nuclear data path> ;} }
  !!
  !! Public Members:
  !!   nuclides  -> array of aceNeutronNuclides with data
  !!   materials -> array of ceNeutronMaterials with data
  !!   Ebounds   -> array with bottom (1) and top (2) energy bound
  !!   activeMat -> array of materials present in the geometry
  !!   nucToZaid -> map to link nuclide index to zaid index
  !!   hasUrr    -> ures probability tables flag, it's false by default
  !!
  !! Interface:
  !!   nuclearData Interface
  !!   ceNeutronDatabase Interface
  !!
  type, public, extends(ceNeutronDatabase) :: aceNeutronDatabase
    type(aceNeutronNuclide),dimension(:),pointer :: nuclides  => null()
    type(ceNeutronMaterial),dimension(:),pointer :: materials => null()
    real(defReal), dimension(2)                  :: Ebounds   = ZERO
    integer(shortInt),dimension(:),allocatable   :: activeMat

    ! Probability tables data
    integer(shortInt),dimension(:),allocatable   :: nucToZaid
    logical(defBool)                             :: hasUrr = .false.

  contains
    ! nuclearData Procedures
    procedure :: kill
    procedure :: matNamesMap
    procedure :: getMaterial
    procedure :: getNuclide
    procedure :: getReaction
    procedure :: init
    procedure :: init_urr
    procedure :: activate

    ! ceNeutronDatabase Procedures
    procedure :: energyBounds
    procedure :: updateTotalMatXS
    procedure :: updateMajorantXS
    procedure :: updateMacroXSs
    procedure :: updateTotalNucXS
    procedure :: updateMicroXSs

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
    else if (MT == 800) then
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
    class(RNG), intent(inout)             :: rand
    integer(shortInt)                     :: i, nucIdx
    real(defReal)                         :: dens
    integer(shortInt), parameter :: IN_SET = 1, NOT_PRESENT = 0

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
    class(RNG), intent(inout)             :: rand
    integer(shortInt)                     :: i, matIdx

    associate (maj => cache_majorantCache(1) )
      maj % E  = E
      maj % xs = ZERO

      do i = 1, size(self % activeMat)
        matIdx = self % activeMat(i)

        ! Update if needed
        if( cache_materialCache(matIdx) % E_tot /= E) then
          call self % updateTotalMatXS(E, matIdx, rand)
        end if

        maj % xs = max(maj % xs, cache_materialCache(matIdx) % xss % total)
      end do
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
    class(RNG), intent(inout)             :: rand
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
        if(cache_nuclideCache(nucIdx) % E_tail /= E) then
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
    class(RNG), intent(inout)             :: rand
    logical(defBool)                      :: needsUrr, needsSab

    associate (nucCache => cache_nuclideCache(nucIdx), &
               nuc      => self % nuclides(nucIdx)     )

      needsUrr = (nuc % hasProbTab .and. E >= nuc % urrE(1) .and. E <= nuc % urrE(2))
      nucCache % needsSabInel = .false.
      nucCache % needsSabEl = .false.
      ! Check if the nuclide needs S(a,b) at this energy
      if (nuc % hasThData) then
        nucCache % needsSabEl = (E >= nuc % SabEl(1) .and. E < nuc % SabEl(2))
        nucCache % needsSabInel = (E >= nuc % SabInel(1) .and. E < nuc % SabInel(2))
      end if

      needsSab = (nucCache % needsSabEl .or. nucCache % needsSabInel)

      if (needsUrr .or. needsSab) then
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
    class(RNG), intent(inout)             :: rand
    logical(defBool)                      :: needsUrr

    associate (nucCache => cache_nuclideCache(nucIdx), &
               nuc      => self % nuclides(nucIdx)     )

      nucCache % E_tail = E
      needsUrr = (nuc % hasProbTab .and. E >= nuc % urrE(1) .and. E <= nuc % urrE(2))

      ! In case the total XS hasn't been retrieved before (during tracking)
      if (nucCache % E_tot /= E) then
        nucCache % E_tot  = E
        call nuc % search(nucCache % idx, nucCache % f, E)
      end if

      ! Overwrites all the micro cross sections in cache
      ! Check if probability tables should be read
      if (needsUrr) then
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
    character(nameLen)                               :: name, name_file
    character(:), allocatable                        :: zaid, file
    integer(shortInt)                                :: i, j, envFlag, nucIdx, idx
    integer(shortInt)                                :: maxNuc
    logical(defBool)                                 :: isFissileMat
    integer(shortInt),dimension(:),allocatable       :: nucIdxs
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
        call self % nuclides(nucIdx) % init_Sab(ACE_Sab)
      end if

      ! Initialise probability tables
      if (self % hasUrr) call self % nuclides(nucIdx) % init_urr(ACE)

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

    if (self % hasUrr) then
       call self % init_urr()
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
  subroutine init_urr(self)
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

  end subroutine init_urr

  !!
  !! Activate this nuclearDatabase
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine activate(self, activeMat)
    class(aceNeutronDatabase), intent(inout)    :: self
    integer(shortInt), dimension(:), intent(in) :: activeMat

    ! Load active materials
    if(allocated(self % activeMat)) deallocate(self % activeMat)
    self % activeMat = activeMat

    ! Configure Cache
    if (self % hasUrr) then
      call cache_init(size( self % materials), size(self % nuclides), 1, maxval(self % nucToZaid))
    else
      call cache_init(size( self % materials), size(self % nuclides))
    end if

  end subroutine activate

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
