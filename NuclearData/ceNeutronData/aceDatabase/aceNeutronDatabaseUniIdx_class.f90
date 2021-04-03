module aceNeutronDatabaseUniIdx_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures, only : fatalError, numToChar, concatenate, removeDuplicatesSorted, &
                                hasDuplicatesSorted, findDuplicatesSorted, quickSort, binarySearch
  use dictionary_class,  only : dictionary
  use RNG_class,         only : RNG
  use charMap_class,     only : charMap

  ! Nuclear Data Interfaces
  use nuclearDatabase_inter,        only : nuclearDatabase
  use materialHandle_inter,         only : materialHandle
  use nuclideHandle_inter,          only : nuclideHandle
  use reactionHandle_inter,         only : reactionHandle
  use ceNeutronDatabase_inter,      only : ceNeutronDatabase, ceNeutronDatabase_CptrCast
  use neutronXSPackages_class,      only : neutronMicroXSs
  use ceNeutronMaterialUni_class,   only : ceNeutronMaterialUni

  ! Material Menu
  use materialMenu_mod,             only : materialItem, nuclideInfo, mm_nMat => nMat, &
                                           mm_getMatPtr => getMatPtr, mm_nameMap => nameMap

  ! ACE CE Nuclear Data Objects
  use aceLibrary_mod,            only : new_neutronAce, aceLib_load => load, aceLib_kill => kill
  use aceCard_class,             only : aceCard
  use aceNeutronNuclide_class,   only : aceNeutronNuclide


  ! CE NEUTRON CACHE
  use ceNeutronCache_mod,           only : cache_nuclideCache => nuclideCache, &
                                           cache_materialCache => materialCache, &
                                           cache_majorantCache => majorantCache, &
                                           cache_init => init

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: aceNeutronDatabaseUniIdx_TptrCast
  public :: aceNeutronDatabaseUniIdx_CptrCast

  !!
  !! A CE Neutron Database based on ACE file format
  !!
  !! As aceNeutronDatabaseUni, this Database builds unionised cross-section grids
  !! per each material. Then, vectors that link the energy index in each unionised
  !! grid to the corresponding nuclide energy grid index are constructed.
  !!
  !! This has similar performance of aceNeutronDatabaseUni, but with better
  !! memory usage.
  !!
  !! Sample input:
  !!   nuclearData {
  !!   handles {
  !!   ce {type aceNeutronDatabaseUniIdx; aceLibrary <nuclear data path> ;} }
  !!
  !! Public Members:
  !!   nuclides  -> array of aceNeutronNuclides with data
  !!   materials -> array of ceNeutronMaterials with data
  !!   Ebounds   -> array with bottom (1) and top (2) energy bound
  !!
  !! Interface:
  !!   nuclearData Interface
  !!   ceNeutronDatabase Interface
  !!
  type, public, extends(ceNeutronDatabase) :: aceNeutronDatabaseUniIdx
    type(aceNeutronNuclide),dimension(:),pointer    :: nuclides  => null()
    type(ceNeutronMaterialUni),dimension(:),pointer :: materials => null()
    real(defReal), dimension(2)                     :: Ebounds   = ZERO
    integer(shortInt),dimension(:),allocatable      :: activeMat
  contains
    ! nuclearData Procedures
    procedure :: kill
    procedure :: matNamesMap
    procedure :: getMaterial
    procedure :: getNuclide
    procedure :: getReaction
    procedure :: init
    procedure :: activate
    procedure :: unionise
    procedure :: unioniseEnergy
    procedure :: interpolateXSs

    ! ceNeutronDatabase Procedures
    procedure :: energyBounds
    procedure :: updateTotalMatXS
    procedure :: updateMajorantXS
    procedure :: updateMacroXSs
    procedure :: updateTotalNucXS
    procedure :: updateMicroXSs

  end type aceNeutronDatabaseUniIdx



contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(aceNeutronDatabaseUniIdx), intent(inout) :: self

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
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    type(charMap), pointer                      :: map

    map => mm_nameMap

  end function matNamesMap

  !!
  !! Return pointer to material in a database
  !!
  !! See nuclearData_inter for  more details
  !!
  function getMaterial(self, matIdx) result(mat)
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    integer(shortInt), intent(in)               :: matIdx
    class(materialHandle), pointer              :: mat

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
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    integer(shortInt), intent(in)               :: nucIdx
    class(nuclideHandle), pointer               :: nuc

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
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    integer(shortInt), intent(in)               :: MT
    integer(shortInt), intent(in)               :: idx
    class(reactionHandle),pointer               :: reac
    integer(shortInt)                           :: idxMT

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
      reac => self % nuclides(idx) % elasticScatter
    else if ( MT == N_fission) then
      reac => self % nuclides(idx) % fission
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
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    real(defReal), intent(out)                  :: E_min
    real(defReal), intent(out)                  :: E_max

    E_min = self % Ebounds(1)
    E_max = self % Ebounds(2)

  end subroutine energyBounds


  !!
  !! Make sure that totalXS of material with matIdx is at energy E
  !! in ceNeutronCache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateTotalMatXS(self, E, matIdx, rand)
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    real(defReal), intent(in)                   :: E
    integer(shortInt), intent(in)               :: matIdx
    class(RNG), intent(inout)                   :: rand
    integer(shortInt)                           :: i, nucIdx
    character(100), parameter :: Here = 'updateTotalMatXS (aceNeutronDatabaseUniIdx_class.f90)'

    associate (matCache => cache_materialCache(matIdx), &
               thisMat => self % materials(matIdx))
      ! Set new energy
      matCache % E_tot  = E

      ! Find index and interpolation factor
      call thisMat % search(matCache % idx, matCache % f, E)

      ! Find current total XS
      matCache % xss % total = thisMat % totalXS(matCache % idx + 1) * matCache % f + &
                          (ONE - matCache % f) * thisMat % totalXS(matCache % idx)

      ! Save the corresponding energy index in cache per each nuclide in the material
      do i = 1,size(thisMat % nuclides)
        nucIdx = thisMat % nuclides(i)
        associate (nucCache => cache_nuclideCache(nucIdx), &
                   nuc      => self % nuclides(nucIdx)     )

          nucCache % E_mat = E
          nucCache % idx = thisMat % gridIdx (i, matCache % idx)

          if ((nucCache % idx == 0) .or. (nucCache % idx == size(nuc % eGrid))) then

            call fatalError (Here, "Index "//numToChar(nucCache % idx)//" out of range &
            in energy grid of nuclide: "//numToChar(nucIdx)//" for energy "//numToChar(E))

          end if

          associate(E_top => nuc % eGrid(nucCache % idx + 1), E_low  => nuc % eGrid(nucCache % idx))
            nucCache % f = (E - E_low) / (E_top - E_low)
          end associate

        end associate
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
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    real(defReal), intent(in)                   :: E
    class(RNG), intent(inout)                   :: rand
    integer(shortInt)                           :: i, matIdx

    associate (maj => cache_majorantCache(1) )
      maj % E  = E
      maj % xs = ZERO

      do i=1,size(self % activeMat)
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
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    real(defReal), intent(in)                   :: E
    integer(shortInt), intent(in)               :: matIdx
    class(RNG), intent(inout)                   :: rand
    integer(shortInt)                           :: i, nucIdx
    real(defReal)                               :: dens

    associate (matCache => cache_materialCache(matIdx))
      ! Set new energy
      matCache % E_tot  = E
      matCache % E_tail = E

      ! Clean current xss
      call matCache % xss % clean()

      ! Construct microscopic XSs
      do i = 1,size(self % materials(matIdx) % nuclides)
        dens   = self % materials(matIdx) % dens(i)
        nucIdx = self % materials(matIdx) % nuclides(i)

        ! Update if needed
        if(cache_nuclideCache(nucIdx) % E_tail /= E) then
          call self % updateMicroXSs(E, nucIdx, rand)
        end if

        ! Add microscopic XSs
        call matCache % xss % add(cache_nuclideCache(nucIdx) % xss, dens)
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
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    real(defReal), intent(in)                   :: E
    integer(shortInt), intent(in)               :: nucIdx
    class(RNG), intent(inout)                   :: rand

    associate (nucCache => cache_nuclideCache(nucIdx), &
               nuc      => self % nuclides(nucIdx)     )

      if (nucCache % E_mat /= E) then
        call nuc % search(nucCache % idx, nucCache % f, E)
      end if

      nucCache % E_tot = E
      nucCache % xss % total = nuc % totalXS(nucCache % idx, nucCache % f)

    end associate

  end subroutine updateTotalNucXS

  !!
  !! Make sure that the microscopic XSs for the nuclide with nucIdx are set
  !! to energy E in ceNeutronCache
  !!
  !! See ceNeutronDatabase for more details
  !!
  subroutine updateMicroXSs(self, E, nucIdx, rand)
    class(aceNeutronDatabaseUniIdx), intent(in) :: self
    real(defReal), intent(in)                   :: E
    integer(shortInt), intent(in)               :: nucIdx
    class(RNG), intent(inout)                   :: rand

    associate (nucCache => cache_nuclideCache(nucIdx), &
               nuc      => self % nuclides(nucIdx)     )

      if (nucCache % E_mat /= E .and. nucCache % E_tot /= E) then
        call nuc % search(nucCache % idx, nucCache % f, E)
      end if

      nucCache % E_tail = E
      ! Overwrites all the micro cross sections in cache
      call nuc % microXSs(nucCache % xss, nucCache % idx, nucCache % f)

    end associate

  end subroutine updateMicroXSs

  !!
  !! Initialise Database from dictionary and pointer to self
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine init(self, dict, ptr, silent )
    class(aceNeutronDatabaseUniIdx), target, intent(inout) :: self
    class(dictionary), intent(in)                          :: dict
    class(nuclearDatabase), pointer, intent(in)            :: ptr
    logical(defBool), optional, intent(in)                 :: silent
    logical(defBool)                                       :: loud
    type(materialItem), pointer                            :: mat
    class(ceNeutronDatabase), pointer                      :: ptr_ceDatabase
    type(charMap)                                          :: nucSet
    type(aceCard)                                          :: ACE
    character(pathLen)                                     :: aceLibPath
    integer(shortInt)                                      :: i, j, envFlag, nucIdx
    integer(shortInt)                                      :: maxNuc
    logical(defBool)                                       :: isFissileMat
    integer(shortInt),dimension(:),allocatable             :: nucIdxs
    integer(shortInt), parameter :: IN_SET = 1, NOT_PRESENT = 0
    character(100), parameter :: Here = 'init (aceNeutronDatabaseUniIdx_class.f90)'

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
    maxNuc = 0
    do i=1,mm_nMat()
      mat => mm_getMatPtr(i)
      maxNuc = max(maxNuc, size(mat % nuclides))

      ! Add all nuclides in material to the map
      do j=1,size(mat % nuclides)
        call nucSet % add(mat % nuclides(j) % toChar(), IN_SET)

      end do
    end do

    ! Get path to ACE library
    call dict % get(aceLibPath,'aceLibrary')

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
      if(loud) then
        print '(A)', "Building: "// trim(nucSet % atKey(i))// " with index: " //numToChar(nucIdx)
      end if

      call new_neutronACE(ACE, nucSet % atKey(i))
      call self % nuclides(nucIdx) % init(ACE, nucIdx, ptr_ceDatabase)

      ! Store nucIdx in the dictionary
      call nucSet % atSet(nucIdx, i)
      nucIdx = nucIdx + 1
      i = nucSet % next(i)
    end do

    ! Build Material definitions
    allocate(self % materials(mm_nMat()))
    allocate(nucIdxs(maxNuc))
    do i=1,mm_nMat()
      mat => mm_getMatPtr(i)

      ! Load nuclide indices on storage space
      isFissileMat = .false.
      do j=1,size(mat % nuclides)
        nucIdxs(j) = nucSet % get( mat % nuclides(j) % toChar())
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

    do i=2,size(self % nuclides)
      self % Ebounds(1) = max(self % Ebounds(1), self % nuclides(i) % eGrid(1))
      j = size(self % nuclides(i) % eGrid)
      self % Ebounds(2) = min(self % Ebounds(2), self % nuclides(i) % eGrid(j))
    end do

    call self % unionise()

    !! Clean up
    call aceLib_kill()

  end subroutine init

  !!
  !! Create the material unionised grid given the individual nuclide grids.
  !!
  !! It concatenates nuclide energy grids, sort the grid and remove duplicate points.
  !!
  !! Args:
  !!   mat [inout] -> pointer to material, its output includes the unionised grid
  !!
  subroutine unioniseEnergy(self,mat)
    class(aceNeutronDatabaseUniIdx),intent(inout) :: self
    type(ceNeutronMaterialUni),intent(inout)      :: mat
    real(defReal),dimension(:),allocatable        :: tmp_grid
    integer(shortInt),dimension(:),allocatable    :: duplicatesIdx
    logical(defBool)                              :: doesIt
    integer(shortInt)                             :: i, j, nucIdx, E_idx

    ! Initialise grid
    tmp_grid = self % nuclides(mat % nuclides(1)) % eGrid

    ! Concatenate energy grids
    do i=2,size(mat % nuclides)
      nucIdx = mat % nuclides(i)

      ! Make sure that nuclide energy grid doesn't have repeated energies
      doesIt = hasDuplicatesSorted(self % nuclides(nucIdx) % eGrid)
      ! If there are duplicate energies, move them by 1.0+100*floatTol
      if (doesIt) then
        duplicatesIdx = findDuplicatesSorted(self % nuclides(nucIdx) % eGrid)

        do j=1,size(duplicatesIdx)
          E_idx = duplicatesIdx(j)
          self % nuclides(nucIdx) % eGrid(E_idx) = &
          (1.0+100*floatTol) * self % nuclides(nucIdx) % eGrid(E_idx-1)
        end do

      end if
      ! Add energy grid to temporary (unsorted) array
      tmp_grid = concatenate(tmp_grid,self % nuclides(nucIdx) % eGrid)
    end do

    ! Sort energies
    call quickSort(tmp_grid)

    ! Remove duplicate energy values
    mat % unionGrid = removeDuplicatesSorted(tmp_grid)

  end subroutine unioniseEnergy

  !!
  !! Interpolate nuclides XSs to match material unionised grids
  !!
  !! Args:
  !!   mat [in]    -> pointer to current material
  !!   nucIdx [in] -> index of the nuclise whose XSs are modified
  !!
  !! Fatal Errors:
  !!   - Grid boundaries have been breached when expanding the nuclide enegy grids to match
  !!     the material unionised one
  !!   - Grid boundaries have been breached when interpolating the nuclide cross-section values
  !!     that correspons to the added energy points
  !!
  subroutine interpolateXSs(self,mat,nucIdx,nucNumber,unionData)
    class(aceNeutronDatabaseUniIdx),intent(inout) :: self
    type(ceNeutronMaterialUni),intent(inout)      :: mat
    integer(shortInt),intent(in)                  :: nucIdx
    integer(shortInt),intent(in)                  :: nucNumber
    real(defReal),dimension(:,:),intent(inout)    :: unionData
    real(defReal),dimension(:),allocatable        :: XS_low, XS_top
    real(defReal)                                 :: f, E_low, den
    integer(shortInt)                             :: h, k, E_idx
    character(100), parameter :: Here = 'interpolateXSs (aceNeutronDatabaseUniIdx_class.f90)'

    associate (nuc => self % nuclides(nucIdx))

      ! STEP 1
      ! Fill XS arrays with either the value corresponding to the energy, or zero
      h = 1
      preLoop: do k=1,size(mat % unionGrid)

        if (abs(nuc % eGrid(h)/mat % unionGrid(k) - ONE) < 1.0e-8) then

          unionData(:,k) = nuc % mainData(:,h)
          if (h < size(nuc % eGrid)) then
            if (nuc % eGrid(h+1) == nuc % eGrid(h)) then
              h = h + 1
            end if
            h = h + 1
          else
            unionData(:,k+1:size(mat % unionGrid)) = ZERO
            mat % gridIdx(nucNumber,k) = h-1
            mat % gridIdx(nucNumber,k+1:size(mat % unionGrid)) = h
            exit preLoop
          end if

        else
          unionData(:,k) = ZERO
        end if
        ! Create double indexing matrix per each material, referring to nuclides' grids
        mat % gridIdx(nucNumber,k) = h - 1

      end do preLoop

      ! Ensure that the indexes are what they are expected to be
      if (h /= size(nuc % eGrid)) then
        call fatalError(Here,'Failed to create unionised grid for isotope '//numToChar(nucIdx))
      end if

      ! STEP 2
      ! The filling starts from the first energy of the nuclide grid
      if (mat % unionGrid(1) /= nuc % eGrid(1)) then
        h = binarySearch(mat % unionGrid,nuc % eGrid(1))
      else
        h = 1
      end if
      k = h + 1
      ! Interpolate to fill all the values where XSs are ZERO
      fillLoop: do while (k <= size(mat % unionGrid))

        do while (unionData(1,k) == ZERO)
          k = k + 1
          if (k > size(mat % unionGrid)) exit fillLoop
        end do

        ! Get interpolation constants
        E_low = mat % unionGrid(h)
        den = (mat % unionGrid(k) - E_low)
        XS_low = unionData(:, h)
        XS_top = unionData(:, k)

        ! Perform interpolation
        do while (h < k-1)
          f = (mat % unionGrid(h+1) - E_low) / den
          unionData(:, h+1) = XS_top * f + (ONE-f) * XS_low
          h = h + 1
        end do
        h = k
        k = k + 1

      end do fillLoop

      ! Ensure that the indexes are what they are expected to be
      if (mat % unionGrid(size(mat % unionGrid)) /= nuc % eGrid(size(nuc % eGrid))) then
        E_idx = binarySearch(mat % unionGrid,nuc % eGrid(size(nuc % eGrid)))
        if (E_idx /= h) then
          call fatalError(Here,'Failed to interpolate unionised grid for isotope '//numToChar(nucIdx))
        end if
      end if

    end associate

    end subroutine interpolateXSs

  !!
  !! Unionise energy grid per material and generate XS
  !!
  subroutine unionise(self)
    class(aceNeutronDatabaseUniIdx),intent(inout) :: self
    type(ceNeutronMaterialUni),pointer            :: mat
    real(defReal),dimension(:,:),allocatable      :: unionData
    integer(shortInt)                             :: i, j, nucIdx


    ! UNIONISE MATERIAL ENERGY GRIDS
    matLoop: do i=1,size(self % materials)

      print*, '...................................................................'
      print*, 'Building energy grid for material '//numToChar(i)

      mat => self % materials(i)

      ! Uninise material energy grids
      call self % unioniseEnergy(mat)

      ! Allocate memory for index matrix and TOTAL XS
      allocate(mat % gridIdx(size(mat % nuclides), size(mat % unionGrid)))
      allocate(mat % totalXS(size(mat % unionGrid)))
      mat % totalXS = ZERO

      ! Fill XS arrays with either the value corresponding to the energy, or zero
      nucLoop: do j=1,size(mat % nuclides)
        nucIdx = mat % nuclides(j)
        allocate(unionData(size(self % nuclides(nucIdx) % mainData(:,1)),size(mat % unionGrid)))

        ! Interpolates the XSs of each nuclide to match the unionised grid
        call self % interpolateXSs(mat,nucIdx,j,unionData)

        ! Precalculate TOTAL XS of each material
        mat % totalXS =  mat % totalXS + mat % dens(j) * unionData(1,:)
        deallocate(unionData)

      end do nucLoop

    end do matLoop

  end subroutine unionise

  !!
  !! Activate this nuclearDatabase
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine activate(self, activeMat)
    class(aceNeutronDatabaseUniIdx), intent(inout) :: self
    integer(shortInt), dimension(:), intent(in)    :: activeMat

    ! Load active materials
    if(allocated(self % activeMat)) deallocate(self % activeMat)
    self % activeMat = activeMat

    ! Configure Cache
    call cache_init(size( self % materials), size(self % nuclides))

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
  pure function aceNeutronDatabaseUniIdx_TptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in)  :: source
    type(aceNeutronDatabaseUniIdx), pointer      :: ptr

    select type(source)
    type is(aceNeutronDatabaseUniIdx)
        ptr => source

      class default
        ptr => null()
    end select

  end function aceNeutronDatabaseUniIdx_TptrCast

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
  pure function aceNeutronDatabaseUniIdx_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(aceNeutronDatabaseUniIdx), pointer    :: ptr

    select type(source)
    class is(aceNeutronDatabaseUniIdx)
        ptr => source

      class default
        ptr => null()
    end select

  end function aceNeutronDatabaseUniIdx_CptrCast



end module aceNeutronDatabaseUniIdx_class
