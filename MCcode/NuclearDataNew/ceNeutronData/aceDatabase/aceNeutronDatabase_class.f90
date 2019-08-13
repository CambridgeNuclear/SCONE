module aceNeutronDatabase_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures, only : fatalError, numToChar
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
  use ceNeutronMaterial_class,      only : ceNeutronMaterial

  ! Material Menu
  use materialMenu_mod,             only : materialItem, nuclideInfo, mm_nMat => nMat, &
                                           mm_getMatPtr => getMatPtr

  ! ACE CE Nuclear Data Objects
  use aceLibrary_mod,               only : new_neutronAce, aceLib_load => load, aceLib_kill => kill
  use aceCard_class,                only : aceCard
  use aceNeutronNuclide_class,      only : aceNeutronNuclide


  ! CE NEUTRON CACHE
  use ceNeutronCache_mod,           only : nuclideCache

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
  !! For now the simplest possible implementation. 
  !!
  !! Public Members: 
  !!
  !!
  !! Interface:
  !!   nuclearData Interface
  !!   ceNeutronDatabase Interface
  !!   
  type, public, extends(ceNeutronDatabase) :: aceNeutronDatabase
    type(aceNeutronNuclide),dimension(:),pointer :: nuclides  => null()
    type(ceNeutronMaterial),dimension(:),pointer :: materials => null()

  contains
    ! nuclearData Procedures 
    procedure :: kill
    procedure :: matNamesMap 
    procedure :: getMaterial 
    procedure :: getNuclide 
    procedure :: getReaction

    ! ceNeutronDatabase Procedures 
    procedure :: energyBounds 
    procedure :: updateTotalMatXS
    procedure :: updateMajorantXS
    procedure :: updateMacroXSs
    procedure :: updateTotalNucXS
    procedure :: updateMicroXSs

    ! This type procedures 
    procedure :: init      
  end type aceNeutronDatabase



contains 

  !!
  !! Return to uninitialised state 
  !!
  elemental subroutine kill(self) 
    class(aceNeutronDatabase), intent(inout) :: self 

    !! IMPLEMENT


  end subroutine kill

  !!
  !! Return pointer to material names map
  !!
  !! See nuclearData_inter for  more details 
  !!
  function matNamesMap(self) result(map)
    class(aceNeutronDatabase), intent(in) :: self
    type(charMap), pointer                :: map
  
    map => null()

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

    mat => null()

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
  
    nuc => null()

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
    class(reactionHandle),pointer         :: reac

    reac => null()

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
  
    E_min = ONE
    E_max = ONE

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
  end subroutine updateMicroXSs

  !!
  !! Initialise Database from dictionary and pointer to self 
  !!
  !!
  !! Args: 
  !!   dict [in] -> Dictionary with the settings
  !!   ptr  [in] -> Pointer to self (instance of the aceNeutronDatabase beeing build) 
  !!     of type nuclearData 
  !!
  !! Errors 
  !!   FatalError is ptr is not assosiated with self 
  !!
  subroutine init(self, dict, ptr) 
    class(aceNeutronDatabase), target, intent(inout) :: self 
    class(dictionary), intent(in)                    :: dict 
    class(nuclearDatabase), pointer, intent(in)      :: ptr
    type(materialItem), pointer                      :: mat
    class(ceNeutronDatabase), pointer                :: ptr_ceDatabase
    type(charMap)                                    :: nucSet
    type(aceCard)                                    :: ACE
    character(pathLen)                               :: aceLibPath
    integer(shortInt)                                :: i, j, envFlag, nucIdx
    integer(shortInt)                                :: maxNuc
    logical(defBool)                                 :: isFissileMat
    integer(shortInt),dimension(:),allocatable       :: nucIdxs
    integer(shortInt), parameter :: IN_SET = 1, NOT_PRESENT = 0
    character(100), parameter :: Here = 'init (aceNeutronDatabase_class.f90)'

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
      print '(A)', "Building: "// trim(nucSet % atKey(i))// " with index: " //numToChar(nucIdx)

      call new_neutronACE(ACE, nucSet % atKey(i))
      call self % nuclides(nucIdx) % init(ACE, nucIdx, ptr_ceDatabase)

      ! Store nucIdx in the dictionary
      call nucSet % add(nucSet % atKey(i), nucIdx)

      nucIdx = nucIdx + 1
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

  end subroutine init 


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
