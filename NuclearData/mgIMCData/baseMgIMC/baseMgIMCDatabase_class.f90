module baseMgIMCDatabase_class

  use numPrecision
  use endfConstants
  use universalVariables, only : VOID_MAT, OUTSIDE_MAT
  use genericProcedures,  only : fatalError, numToChar
  use particle_class,     only : particle
  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : fileToDict
  use RNG_class,          only : RNG

  ! Nuclear Data Interfaces
  use nuclearDatabase_inter,   only : nuclearDatabase
  use mgIMCDatabase_inter,     only : mgIMCDatabase
  use materialHandle_inter,    only : materialHandle
  use nuclideHandle_inter,     only : nuclideHandle
  use reactionHandle_inter,    only : reactionHandle
  use materialMenu_mod,        only : materialItem, mm_getMatPtr => getMatPtr, mm_nMat => nMat, &
                                      mm_nameMap => nameMap, mm_matName => matName

  ! baseMgIMC Objects
  use baseMgIMCMaterial_class, only : baseMgIMCMaterial
  use materialEquations,       only : mgEnergyGrid

  ! Energy grid for multi-frequency
  use energyGrid_class,         only : energyGrid
  use energyGridRegistry_mod,   only : get_energyGrid

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: baseMgIMCDatabase_TptrCast
  public :: baseMgIMCDatabase_CptrCast

  !!
  !! Database for MG IMC and ISMC calculations
  !!
  !! All materials in a problem are of class baseMgIMCMaterial. See its documentation for
  !! details on how the physics is handled
  !!
  !! Public Members:
  !!   mats       -> array containing all defined materials (by matIdx)
  !!   activeMats -> list of matIdxs of materials active in the problem
  !!   nG         -> number of energy groups
  !!
  !! Interface:
  !!   nuclearDatabase interface
  !!
  type, public, extends(mgIMCDatabase) :: baseMgIMCDatabase
    type(baseMgIMCMaterial), dimension(:), pointer     :: mats => null()
    integer(shortInt), dimension(:), allocatable       :: activeMats
    integer(shortInt)                                  :: nG = 0
    type(energyGrid), private                          :: eGrid

  contains
    ! Superclass Interface
    procedure :: getTransMatXS
    procedure :: getTotalMatXS
    procedure :: getMajorantXS
    procedure :: matNamesMap
    procedure :: getMaterial
    procedure :: getNuclide
    procedure :: getReaction
    procedure :: getEmittedRad
    procedure :: getMaterialEnergy
    procedure :: updateProperties
    procedure :: setCalcType
    procedure :: sampleTransformTime
    procedure :: sampleEnergyGroup
    procedure :: kill
    procedure :: init
    procedure :: activate

    ! Local interface
    procedure :: nGroups

  end type baseMgIMCDatabase

contains

  !!
  !! Get Transport XS given a particle
  !!
  !! See nuclearDatabase documentation for details
  !!
  !! Note:
  !!   DOES NOT check if particle is MG. Will refer to G in the particle and give error
  !!   if the value is invalid
  !!
  function getTransMatXS(self, p, matIdx) result(xs)
    class(baseMgIMCDatabase), intent(inout)     :: self
    class(particle), intent(in)                 :: p
    integer(shortInt), intent(in)               :: matIdx
    real(defReal)                               :: xs

    xs = self % getTotalMatXS(p, matIdx)

  end function getTransMatXS

  !!
  !! Get Total XS given a particle
  !!
  !! See nuclearDatabase documentation for details
  !!
  !! Note:
  !!   DOES NOT check if particle is MG. Will refer to G in the particle and give error
  !!   if the value is invalid
  !!
  function getTotalMatXS(self, p, matIdx) result(xs)
    class(baseMgIMCDatabase), intent(inout)     :: self
    class(particle), intent(in)                 :: p
    integer(shortInt), intent(in)               :: matIdx
    real(defReal)                               :: xs
    character(100), parameter :: Here = 'getTotalMatXS (baseMgIMCDatabase_class.f90)'

    if (matIdx == OUTSIDE_MAT) call fatalError(Here, 'Requesting XS in OUTSIDE_MAT')

    if (matIdx == VOID_MAT) then
      xs = ZERO
    else
      xs = self % mats(matIdx) % getTotalXS(p % G, p % pRNG)
    end if

  end function getTotalMatXS

  !!
  !! Get Majorant XS given a particle
  !!
  !! See nuclearDatabase documentation for details
  !!
  !! Note:
  !!   DOES NOT check if particle is MG. Will refer to G in the particle and give error
  !!   if the value is invalid
  !!
  function getMajorantXS(self, p) result(xs)
    class(baseMgIMCDatabase), intent(inout)     :: self
    class(particle), intent(in)                 :: p
    real(defReal)                               :: xs
    integer(shortInt)                           :: i, idx

    xs = ZERO
    do i=1,size(self % activeMats)
      idx = self % activeMats(i)
      xs = max(xs, self % getTotalMatXS(p, idx))
    end do

  end function getMajorantXS

  !!
  !! Return pointer to material names map
  !!
  !! See nuclearDatabase documentation for details
  !!
  function matNamesMap(self) result(map)
    class(baseMgIMCDatabase), intent(in)     :: self
    type(charMap), pointer                   :: map

    map => mm_nameMap

  end function matNamesMap

  !!
  !! Return pointer to a material in the database
  !!
  !! See nuclearDatabase documentation for details
  !!
  function getMaterial(self, matIdx) result(mat)
    class(baseMgIMCDatabase), intent(in)     :: self
    integer(shortInt), intent(in)            :: matIdx
    class(materialHandle), pointer           :: mat

    if(matIdx < 1 .or. matIdx > size(self % mats)) then
      mat => null()
    else
      mat => self % mats(matIdx)
    end if

  end function getMaterial

  !!
  !! Return pointer to a nuclide in the database
  !!
  !! See nuclearDatabase documentation for details
  !!
  !! Note:
  !!   This database has no nucldie. Returns NULL always!
  !!
  function getNuclide(self, nucIdx) result(nuc)
    class(baseMgIMCDatabase), intent(in)     :: self
    integer(shortInt), intent(in)            :: nucIdx
    class(nuclideHandle), pointer            :: nuc

    nuc => null()

  end function getNuclide

  !!
  !! Return pointer to a reaction
  !!
  !! See nuclearDatabase documentation for details
  !!
  function getReaction(self, MT, idx) result(reac)
    class(baseMgIMCDatabase), intent(in)     :: self
    integer(shortInt), intent(in)            :: MT
    integer(shortInt), intent(in)            :: idx
    class(reactionHandle), pointer           :: reac
    character(100), parameter                :: Here = 'getReaction (baseMgIMCDatabase_class.f90)'

    reac => null()
    call fatalError(Here, "Pointless function call")

  end function getReaction

  !!
  !! Return energy to be emitted during current time step
  !!
  !! Args:
  !!   matIdx [in] [optional] -> If provided, return the energy to be emitted from only matIdx
  !!                             Otherwise, return the total energy to be emitted from all mats
  !!
  function getEmittedRad(self, matIdx) result(energy)
    class(baseMgIMCDatabase), intent(in)    :: self
    integer(shortInt), intent(in), optional :: matIdx
    real(defReal)                           :: energy
    integer(shortInt)                       :: i

    ! If matIdx provided, return radiation emitted from only that material
    if (present(matIdx)) then
      energy = self % mats(matIdx) % getEmittedRad()
      return
    end if

    ! Otherwise, return total energy emitted from all materials
    energy = 0

    do i=1, size(self % mats)
      energy = energy + self % mats(i) % getEmittedRad()
    end do

  end function getEmittedRad

  !!
  !! Return material energy
  !!
  !! Args:
  !!   matIdx [in] [optional] -> If provided, return the energy of only matIdx
  !!                             Otherwise, return the total energy of all mats
  !!
  function getMaterialEnergy(self, matIdx) result(energy)
    class(baseMgIMCDatabase), intent(in)    :: self
    integer(shortInt), intent(in), optional :: matIdx
    real(defReal)                           :: energy
    integer(shortInt)                       :: i

    ! If matIdx provided, return radiation emitted from only that material
    if (present(matIdx)) then
      energy = self % mats(matIdx) % getMatEnergy()
      return
    end if

    ! Otherwise, return total energy emitted from all materials
    energy = 0

    do i=1, size(self % mats)
      energy = energy + self % mats(i) % getMatEnergy()
    end do

  end function getMaterialEnergy

  !!
  !! Update material properties based on energy absorbed during the time step
  !!
  subroutine updateProperties(self, tallyEnergy, printUpdates)
    class(baseMgIMCDatabase), intent(inout) :: self
    real(defReal), dimension(:), intent(in) :: tallyEnergy
    integer(shortInt), intent(in)           :: printUpdates
    integer(shortInt)                       :: i
    character(100), parameter :: Here = 'updateProperties (baseMgIMCDatabase_class.f90)'

    ! Check for valid inputs
    if (size(tallyEnergy) /= size(self % mats)) call fatalError(Here, &
                                &'Energy tally array must have size nMats')
    if (printUpdates > size(self % mats)) call fatalError(Here, &
                                &'printUpdates must be <= nMats')

    ! Update mats to be printed (if any), not in parallel to allow correct order of console output
    do i = 1, printUpdates
      print *
      print *, '  Material update for '//trim(mm_matName(i))//':'
      call self % mats(i) % updateMat(tallyEnergy(i), .true.)
      if (i == printUpdates) print *
    end do

    ! Update remaining mats
    !$omp parallel do
    do i = (printUpdates+1), size(tallyEnergy)
      call self % mats(i) % updateMat(tallyEnergy(i))
    end do
    !$omp end parallel do

  end subroutine updateProperties

  !!
  !! Tell each material if we are using IMC or ISMC
  !!
  subroutine setCalcType(self, type)
    class(baseMgIMCDatabase), intent(inout) :: self
    integer(shortInt), intent(in)           :: type
    integer(shortInt)                       :: i

    do i=1, size(self % mats)
      call self % mats(i) % setCalcType(type)
    end do

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
    class(baseMgIMCDatabase), intent(inout) :: self
    integer(shortInt), intent(in)           :: matIdx
    class(RNG), intent(inout)               :: rand
    integer(shortInt)                       :: G

    G = self % mats(matIdx) % sampleEnergyGroup(rand)

  end function sampleEnergyGroup

  !!
  !! Sample the time taken for a material particle to transform into a photon
  !! Used for ISMC only
  !!
  function sampleTransformTime(self, matIdx, rand) result(t)
    class(baseMgIMCDatabase), intent(inout) :: self
    integer(shortInt), intent(in)           :: matIdx
    class(RNG), intent(inout)               :: rand
    real(defReal)                           :: t

    t = self % mats(matIdx) % sampleTransformTime(rand)

  end function sampleTransformTime

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(baseMgIMCDatabase), intent(inout) :: self

    if(associated(self % mats)) then
      call self % mats % kill()
      deallocate(self % mats)
    end if

    if(allocated(self % activeMats)) deallocate (self % activeMats)
    self % nG = 0

  end subroutine kill

  !!
  !! Initialise Database from dictionary and pointer to self
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine init(self, dict, ptr, silent)
    class(baseMgIMCDatabase), target,intent(inout)     :: self
    class(dictionary), intent(in)                      :: dict
    class(nuclearDatabase), pointer,intent(in)         :: ptr
    logical(defBool), intent(in), optional             :: silent
    logical(defBool)                                   :: loud
    integer(shortInt)                                  :: i, nMat
    type(materialItem), pointer                        :: matDef
    character(pathLen)                                 :: path
    type(dictionary)                                   :: tempDict
    logical(defBool)                                   :: err
    character(nameLen)                                 :: gridName
    character(100), parameter :: Here = 'init (baseMgIMCDatabase_class.f90)'
 
    ! Prevent reallocations
    call self % kill()

    ! Set build console output flag
    if(present(silent)) then
      loud = .not.silent
    else
      loud = .true.
    end if

    ! Get energy grid in multi-frequency case
    gridName = 'mg'
    call get_energyGrid(self % eGrid, gridName, err)
    if (err .eqv. .false.) then
      self % nG = self % eGrid % getSize()
      mgEnergyGrid => self % eGrid
    else
      self % nG = 1
    end if

    ! Find number of materials and allocate space
    nMat = mm_nMat()

    allocate(self % mats(nMat))

    ! Build materials
    do i=1,nMat

      ! Get Path to the xsFile
      matDef => mm_getMatPtr(i)
      call matDef % extraInfo % get(path,'xsFile')

      ! Print status
      if(loud) then
        print '(A)', "Building material: " // trim(matDef % name) // " From: " // trim(path)
      end if

      ! Load dictionary
      call fileToDict(tempDict, path)

      ! Add temperature and volume into dictionary
      call tempDict % store('T', matDef % T)
      call tempDict % store('V', matDef % V)

      ! Initialise material
      call self % mats(i) % init(tempDict)

      ! Verify number of groups
      if(self % nG /= self % mats(i) % nGroups()) then
        call fatalError(Here,'Inconsistent # of groups in materials in matIdx'//numToChar(i))
      end if

    end do

  end subroutine init

  !!
  !! Activate this nuclearDatabase
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine activate(self, activeMat)
    class(baseMgIMCDatabase), intent(inout) :: self
    integer(shortInt), dimension(:), intent(in) :: activeMat

    if(allocated(self % activeMats)) deallocate(self % activeMats)
    self % activeMats = activeMat

  end subroutine activate

  !!
  !! Return number of energy groups in this database
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  pure function nGroups(self) result(nG)
    class(baseMgIMCDatabase), intent(in)     :: self
    integer(shortInt)                        :: nG

    nG = self % nG

  end function nGroups

  !!
  !! Cast nuclearDatabase pointer to baseMgIMCDatabase type pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of baseMgIMCDatabase type
  !!   Target points to source if source is baseMgIMCDatabasetype
  !!
  pure function baseMgIMCDatabase_TptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    type(baseMgIMCDatabase), pointer            :: ptr

    select type(source)
      type is(baseMgIMCDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgIMCDatabase_TptrCast

  !!
  !! Cast nuclearDatabase pointer to baseMgIMCDatabase class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of baseMgIMCDatabase class
  !!   Target points to source if source is baseMgIMCDatabase class
  !!
  pure function baseMgIMCDatabase_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(baseMgIMCDatabase), pointer           :: ptr

    select type(source)
      class is(baseMgIMCDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgIMCDatabase_CptrCast


end module baseMgIMCDatabase_class
