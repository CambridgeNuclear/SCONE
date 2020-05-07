module baseMgNeutronDatabase_class

  use numPrecision
  use endfConstants
  use genericProcedures,  only : fatalError, numToChar
  use particle_class,     only : particle
  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : fileToDict

  ! Nuclear Data Interfaces
  use nuclearDatabase_inter,   only : nuclearDatabase
  use mgNeutronDatabase_inter, only : mgNeutronDatabase
  use materialHandle_inter,    only : materialHandle
  use nuclideHandle_inter,     only : nuclideHandle
  use reactionHandle_inter,    only : reactionHandle
  use materialMenu_mod,        only : materialItem, mm_getMatPtr => getMatPtr, mm_nMat => nMat, &
                                      mm_nameMap => nameMap

  ! baseMgNeutron Objects
  use baseMgNeutronMaterial_class, only : baseMgNeutronMaterial

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: baseMgNeutronDatabase_TptrCast
  public :: baseMgNeutronDatabase_CptrCast

  !!
  !! Basic type of MG nuclear Data for neutrons
  !!
  !! All materials in aproblem are baseMgMaterials. See its documentation for
  !! details on how the physics is handled
  !!
  !! Public Members:
  !!   mats       -> array containing all defined materials (by matIdx)
  !!   activeMats -> list of matIdxs of materials active in the problem
  !!
  !! Interface:
  !!   nuclearDatabase interface
  !!
  type, public, extends(mgNeutronDatabase) :: baseMgNeutronDatabase
    type(baseMgNeutronMaterial), dimension(:), pointer :: mats => null()
    integer(shortInt), dimension(:), allocatable       :: activeMats
    integer(shortInt)                                  :: nG = 0

  contains
    ! Superclass Interface
    procedure :: getTransMatXS
    procedure :: getTotalMatXS
    procedure :: getMajorantXS
    procedure :: matNamesMap
    procedure :: getMaterial
    procedure :: getNuclide
    procedure :: getReaction
    procedure :: kill
    procedure :: init
    procedure :: activate

    ! Local interface
    procedure :: nGroups

  end type baseMgNeutronDatabase

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
  !! Sample input dictionary:
  !!   nucData {
  !!     type baseMgNeutronDatabase;
  !!     PN P0;                        // or P1
  !!   }
  !!
  function getTransMatXS(self, p, matIdx) result(xs)
    class(baseMgNeutronDatabase), intent(inout) :: self
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
    class(baseMgNeutronDatabase), intent(inout) :: self
    class(particle), intent(in)                 :: p
    integer(shortInt), intent(in)               :: matIdx
    real(defReal)                               :: xs

    xs = self % mats(matIdx) % getTotalXS(p % G, p % pRNG)

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
    class(baseMgNeutronDatabase), intent(inout) :: self
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
    class(baseMgNeutronDatabase), intent(in) :: self
    type(charMap), pointer                   :: map

    map => mm_nameMap

  end function matNamesMap

  !!
  !! Return pointer to a material in the database
  !!
  !! See nuclearDatabase documentation for details
  !!
  function getMaterial(self, matIdx) result(mat)
    class(baseMgNeutronDatabase), intent(in) :: self
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
    class(baseMgNeutronDatabase), intent(in) :: self
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
    class(baseMgNeutronDatabase), intent(in) :: self
    integer(shortInt), intent(in)            :: MT
    integer(shortInt), intent(in)            :: idx
    class(reactionHandle), pointer           :: reac

    ! Catch Invalid index
    if(idx < 1 .or. idx > size(self % mats)) then
      reac => null()
      return
    end if

    ! Select correct reaction
    select case(MT)
      case(macroFission)
        ! Point to null if material is not fissile
        if (self % mats(idx) % isFissile()) then
          reac => self % mats(idx) % fission
        else
          reac => null()
        end if

      case(macroIEScatter)
        reac => self % mats(idx) % scatter

      case default
        reac => null()

    end select

  end function getReaction

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(baseMgNeutronDatabase), intent(inout) :: self

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
    class(baseMgNeutronDatabase), target,intent(inout) :: self
    class(dictionary), intent(in)                      :: dict
    class(nuclearDatabase), pointer,intent(in)         :: ptr
    logical(defBool), intent(in), optional             :: silent
    logical(defBool)                                   :: loud
    integer(shortInt)                                  :: i, nMat
    type(materialItem), pointer                        :: matDef
    character(pathLen)                                 :: path
    character(nameLen)                                 :: scatterKey
    type(dictionary)                                   :: tempDict
    character(100), parameter :: Here = 'init (baseMgNeutronDatabase_class.f90)'

    ! Prevent reallocations
    call self % kill()

    ! Set build console output flag
    if(present(silent)) then
      loud = .not.silent
    else
      loud = .true.
    end if

    ! Find number of materials and allocate space
    nMat = mm_nMat()

    allocate(self % mats(nMat))

    ! Read scatterKey
    call dict % get(scatterKey, 'PN')

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
      call self % mats(i) % init(tempDict, scatterKey)

    end do

    ! Load and verify number of groups
    self % nG = self % mats(1) % nGroups()
    do i=2,nMat
      if(self % nG /= self % mats(i) % nGroups()) then
        call fatalError(Here,'Inconsistant # of groups in materials in matIdx'//numToChar(i))
      end if
    end do

  end subroutine init

  !!
  !! Activate this nuclearDatabase
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine activate(self, activeMat)
    class(baseMgNeutronDatabase), intent(inout) :: self
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
    class(baseMgNeutronDatabase), intent(in) :: self
    integer(shortInt)                        :: nG

    nG = self % nG

  end function nGroups

  !!
  !! Cast nuclearDatabase pointer to baseMgNeutronDatabase type pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of baseMgNeutronDatabase type
  !!   Target points to source if source is baseMgNeutronDatabasetype
  !!
  pure function baseMgNeutronDatabase_TptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    type(baseMgNeutronDatabase), pointer           :: ptr

    select type(source)
      type is(baseMgNeutronDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgNeutronDatabase_TptrCast

  !!
  !! Cast nuclearDatabase pointer to baseMgNeutronDatabase class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class nuclearDatabase
  !!
  !! Result:
  !!   Null if source is not of baseMgNeutronDatabase class
  !!   Target points to source if source is baseMgNeutronDatabase class
  !!
  pure function baseMgNeutronDatabase_CptrCast(source) result(ptr)
    class(nuclearDatabase), pointer, intent(in) :: source
    class(baseMgNeutronDatabase), pointer          :: ptr

    select type(source)
      class is(baseMgNeutronDatabase)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgNeutronDatabase_CptrCast


end module baseMgNeutronDatabase_class
