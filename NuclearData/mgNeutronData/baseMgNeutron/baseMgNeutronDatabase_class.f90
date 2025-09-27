module baseMgNeutronDatabase_class

  use numPrecision
  use endfConstants
  use universalVariables
  use errors_mod,         only : fatalError
  use genericProcedures,  only : numToChar
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

  ! Cache
  use mgNeutronCache_mod,           only : materialCache, trackingCache, &
                                           cache_init => init

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
  !! All materials in a problem are baseMgMaterials. See its documentation for
  !! details on how the physics is handled
  !!
  !! Sample input dictionary:
  !!   nucData {
  !!     type baseMgNeutronDatabase;
  !!     PN P0;                        // or P1
  !!     #avgDist 2.718; #
  !!   }
  !!
  !! Public Members:
  !!   mats       -> array containing all defined materials (by matIdx)
  !!   majorant   -> majorant xs for delta tracking
  !!   activeMats -> list of matIdxs of materials active in the problem
  !!
  !! Interface:
  !!   nuclearDatabase interface
  !!
  type, public, extends(mgNeutronDatabase) :: baseMgNeutronDatabase
    type(baseMgNeutronMaterial), dimension(:), pointer :: mats => null()
    real(defReal), dimension(:), allocatable           :: majorant
    integer(shortInt), dimension(:), allocatable       :: activeMats

  contains
    ! Superclass Interface
    procedure :: getTrackingXS
    procedure :: getTrackMatXS
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
    procedure :: initMajorant
    procedure :: nGroups

  end type baseMgNeutronDatabase

contains

  !!
  !! Get tracking XS requested
  !!
  !! See nuclearDatabase documentation for details
  !!
  !! Note:
  !!   DOES NOT check if particle is MG. Will refer to G in the particle and give error
  !!   if the value is invalid
  !!
  function getTrackingXS(self, p, matIdx, what) result(xs)
    class(baseMgNeutronDatabase), intent(inout) :: self
    class(particle), intent(in)                 :: p
    integer(shortInt), intent(in)               :: matIdx
    integer(shortInt), intent(in)               :: what
    real(defReal)                               :: xs
    character(100),parameter :: Here = 'getTrackingXS (baseMgNeutronDatabase_class.f90)'

    ! Process request
    select case(what)

      case (MATERIAL_XS)
        if (matIdx == VOID_MAT) then
          xs = self % collisionXS
        else
          xs = max(self % getTrackMatXS(p, matIdx), self % collisionXS)
        end if

      case (MAJORANT_XS)
        xs = max(self % getMajorantXS(p), self % collisionXS)

      case (TRACKING_XS)

        ! READ ONLY - read from previously updated cache
        if (p % G == trackingCache(1) % G) then
          xs = trackingCache(1) % xs
          return
        else
          call fatalError(Here, 'Tracking cache failed to update during tracking')
        end if

      case default
        call fatalError(Here, 'Neither material nor majorant xs was asked')

    end select

    ! Update Cache
    trackingCache(1) % G  = p % G
    trackingCache(1) % xs = xs

  end function getTrackingXS

  !!
  !! Get tracking XS given a particle. In MG, it is always identical to the material
  !! total XS.
  !!
  !! See nuclearDatabase documentation for details
  !!
  function getTrackMatXS(self, p, matIdx) result(xs)
    class(baseMgNeutronDatabase), intent(inout) :: self
    class(particle), intent(in)                 :: p
    integer(shortInt), intent(in)               :: matIdx
    real(defReal)                               :: xs
    character(100),parameter :: Here = 'getTrackMatXS (baseMgNeutronDatabase_class.f90)'

    ! Check that matIdx exists
    if (matIdx < 1 .or. matIdx > mm_nMat()) then 
      print *,'Particle location: ', p % rGlobal()
      call fatalError(Here, 'Particle is in an undefined material with index: '&
              //numToChar(matIdx))
    end if
    
    xs = self % getTotalMatXS(p, matIdx)

  end function getTrackMatXS

  !!
  !! Get total XS given a particle
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
    character(100),parameter :: Here = 'getTotalMatXS (baseMgNeutronDatabase_class.f90)'
    
    ! Check that matIdx exists
    if (matIdx < 1 .or. matIdx > mm_nMat()) then 
      print *,'Particle location: ', p % rGlobal()
      call fatalError(Here, 'Particle is in an undefined material with index: '&
              //numToChar(matIdx))
    end if

    associate (matCache => materialCache(matIdx))

      if (matCache % G_tot /= p % G) then
        ! Get cross section
        xs = self % mats(matIdx) % getTotalXS(p % G, p % pRNG)
        ! Update cache
        matCache % xss % total = xs
        matCache % G_tot = p % G

      else
        ! Retrieve cross section from cache
        xs = matCache % xss % total

      end if

    end associate

  end function getTotalMatXS

  !!
  !! Get majorant XS given a particle
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
    character(100), parameter :: Here = ' getMajorantXS (baseMgNeutronDatabase_class.f90)'

    ! Verify bounds
    if (p % G < 1 .or. self % nG < p % G) then
      call fatalError(Here,'Invalid group number: '//numToChar(p % G)// &
                           ' Data has only: ' // numToChar(self % nG))
      xs = ZERO ! Avoid warning
    end if

    xs = self % majorant(p % G)

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

    if (matIdx < 1 .or. matIdx > size(self % mats)) then
      mat => null()
    else
      ! Retrieve pointer from cache
      mat => materialCache(matIdx) % mat
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
    real(defReal)                                      :: temp
    character(100), parameter :: Here = 'init (baseMgNeutronDatabase_class.f90)'

    ! Prevent reallocations
    call self % kill()

    ! Set build console output flag
    if(present(silent)) then
      loud = .not.silent
    else
      loud = .true.
    end if
    
    ! Check for a minimum average collision distance
    if (dict % isPresent('avgDist')) then
      call dict % get(temp, 'avgDist')

      if (temp <= ZERO) then
        call fatalError(Here, 'Must have a finite, positive minimum average collision distance')
      end if

      self % collisionXS = ONE / temp

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
    do i = 2,nMat
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
  subroutine activate(self, activeMat, silent)
    class(baseMgNeutronDatabase), intent(inout) :: self
    integer(shortInt), dimension(:), intent(in) :: activeMat
    logical(defBool), optional, intent(in)      :: silent
    logical(defBool)                            :: loud
    integer(shortInt)                           :: idx

    if(allocated(self % activeMats)) deallocate(self % activeMats)
    self % activeMats = activeMat

    ! Initialies cross section cache
    call cache_init(size(self % mats))

    ! Store the material pointer in the material cache
    !$omp parallel
    do idx = 1,size(self % mats)
      materialCache(idx) % mat => self % mats(idx)
    end do
    !$omp end parallel

    ! Set build console output flag
    if (present(silent)) then
      loud = .not. silent
    else
      loud = .true.
    end if

    ! Build unionised majorant
    call self % initMajorant(loud)

  end subroutine activate

  !!
  !! Precomputes majorant cross section
  !!
  !! See nuclearDatabase documentation for details
  !!
  subroutine initMajorant(self, loud)
    class(baseMgNeutronDatabase), intent(inout) :: self
    logical(defBool), intent(in)                :: loud
    integer(shortInt)                           :: g, i, idx
    real(defReal)                               :: xs
    integer(shortInt), parameter                :: TOTAL_XS = 1

    ! Allocate majorant
    allocate (self % majorant(self % nG))

    ! Loop over energy groups
    do g = 1,self % nG
      xs = ZERO
      do i = 1,size(self % activeMats)
        idx = self % activeMats(i)
        xs = max(xs, self % mats(idx) % data(TOTAL_XS, g))
      end do
      self % majorant(g) = xs
    end do

    if (loud) print '(A)', 'MG unionised majorant cross section calculation completed'

  end subroutine initMajorant

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
    type(baseMgNeutronDatabase), pointer        :: ptr

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
