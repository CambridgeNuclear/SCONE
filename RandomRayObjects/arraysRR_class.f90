module arraysRR_class

  use numPrecision
  use universalVariables
  use constantsRR
  use genericProcedures,              only : fatalError, numToChar, rotateVector, printFishLineR
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Data
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use dataRR_class,                   only : dataRR

  ! Geometry
  use coord_class,                    only : coordList
  use geometryStd_class,              only : geometryStd

  ! Visualisation
  use tallyMap_inter,                 only : tallyMap
  use visualiser_class,               only : visualiser
  use particle_class,                 only : particleState
  
  ! For locks
  use omp_lib

  implicit none
  private
  
  
  !!
  !! Object to store all arrays in random ray
  !! By default, will have the current flux, previous flux, accumulated flux,
  !! source arrays, and geometric arrays.
  !!
  !! Can be extended to having corresponding flux and source moments, as well
  !! as the fixed source and additional geometric info.
  !!
  !! Private Members
  !!   nG             -> Number of energy groups, kept for convenience.
  !!   nCells         -> Number of unique cells in the geometry, kept for convenience.
  !!   lengthPerIt    -> RR active length per iteration, kept for convenience
  !!   XSData         -> Pointer to nuclear data, for convenience.
  !!   geom           -> Pointer to geometry, for convenience.
  !!   rho            -> Stabilisation factor: 0 is no stabilisation, 1 is aggressive stabilisation
  !!   ani            -> Order of anisotropic flux moments to be stored
  !!   simulationType -> Identifies which simulation to perform: flat/linear, isotropic/anisotropic
  !!
  !!   scalarFlux     -> Array of scalar flux values of length [nG * nCells]
  !!   prevFlux       -> Array of previous scalar flux values of length [nG * nCells]
  !!   fluxScores     -> Array of scalar flux values and squared values to be reported 
  !!                     in results [nG * nCells, 2]
  !!
  !!   source      -> Array of sources [nG * nCells]
  !!   fixedSource -> Array of fixed sources [nG * nCells]
  !!   sourceIdx   -> Array of material indices containing fixed sources
  !!
  !!   volumeTracks  -> Array of sum of track lengths for computing volumes [nCells]
  !!   lengthSquared -> Array of sum of lengths^2 for use in cells with zero XS [nCells]
  !!   volume        -> Array of dimensionless cell volumes [nCells]
  !!
  !!   cellHit   -> Array of ints whether a cell was visited this iteration [nCells]
  !!   cellFound -> Array of logicals whether a cell was ever visited [nCells]
  !!   cellPos   -> Array of cell positions [3 * nCells]
  !!
  !!   locks -> Array of OpenMP locks for each geometric cell
  !!
  type, public :: arraysRR
    private
    ! Components
    class(geometryStd), pointer           :: geom        => null()
    type(dataRR)                          :: XSData      
    integer(shortInt)                     :: nG          = 0
    integer(shortInt)                     :: nCells      = 0
    real(defReal)                         :: lengthPerIt = ZERO
    real(defFlt)                          :: rho         = 0.0_defFlt
    integer(shortInt)                     :: ani         = 0
    integer(shortInt)                     :: simulationType = 0
    
    ! Flux arrays
    real(defFlt), dimension(:), allocatable    :: scalarFlux
    real(defFlt), dimension(:), allocatable    :: prevFlux
    real(defReal), dimension(:,:), allocatable :: fluxScores
    
    ! Source arrays
    real(defFlt), dimension(:), allocatable      :: source
    real(defFlt), dimension(:), allocatable      :: fixedSource
    integer(shortInt), dimension(:), allocatable :: sourceIdx

    ! Geometry arrays
    real(defReal), dimension(:), allocatable     :: volumeTracks
    real(defReal), dimension(:), allocatable     :: lengthSquared
    real(defReal), dimension(:), allocatable     :: volume
    integer(shortInt), dimension(:), allocatable :: cellHit
    logical(defBool), dimension(:), allocatable  :: cellFound
    real(defReal), dimension(:,:), allocatable   :: cellPos
    
    ! OMP locks
    integer(kind=omp_lock_kind), dimension(:), allocatable :: locks

  contains
    
    ! Public procedures
    procedure :: init
    procedure :: kill

    ! Access procedures
    procedure :: getDataPointer
    procedure :: getGeomPointer
    procedure :: getSourcePointer
    procedure :: getFluxPointer
    procedure :: getSource
    procedure :: getPrevFlux
    procedure :: getFluxScore
    procedure :: getFluxSD
    procedure :: getNG
    procedure :: getCellPos
    procedure :: hasHit
    procedure :: getCellHitRate
    procedure :: getSimulationType
    procedure :: found

    ! Change individual elements of the type
    ! Predominantly for use in the transport sweep
    procedure :: incrementVolume
    procedure :: incrementLengthSquared
    procedure :: hitCell
    procedure :: wipeCellHits
    procedure :: newFound
    procedure :: setLock
    procedure :: unsetLock

    ! Basic RR procedures
    procedure :: resetFluxes
    procedure :: normaliseFluxAndVolume
    procedure :: updateSource
    procedure :: accumulateFluxScores
    procedure :: finaliseFluxScores
    procedure :: calculateKeff

    ! Output procedures
    procedure :: outputMap
    procedure :: outputToVTK

    ! Private procedures
    procedure, private :: initialiseFixedSource
    
    procedure, private :: resetFluxesFlatIso
    procedure, private :: resetFluxesLinIso
    procedure, private :: resetFluxesLIFA
    procedure, private :: resetFluxesFlatAni
    
    procedure, private :: normaliseFluxAndVolumeFlatIso
    procedure, private :: normaliseFluxAndVolumeLinIso
    procedure, private :: normaliseFluxAndVolumeLIFA
    procedure, private :: normaliseFluxAndVolumeFlatAni
    
    procedure, private :: sourceUpdateKernelFlatIso
    procedure, private :: sourceUpdateKernelLinIso
    procedure, private :: sourceUpdateKernelLIFA
    procedure, private :: sourceUpdateKernelFlatAni
    
    procedure, private :: calculateKeffKernel

  end type arraysRR

contains

  !!
  !! Initialise the arrays object
  !!
  !! The object is fed sizes and requirements by the physics package.
  !! This will allocate the necessary arrays
  !!
  subroutine init(self, db, geom, lengthPerIt, rho, lin, ani, doKinetics, loud, dictFS)
    class(arraysRR), intent(inout)                      :: self
    class(baseMgNeutronDatabase), pointer, intent(in)   :: db
    class(geometryStd), pointer, intent(in)             :: geom
    real(defReal), intent(in)                           :: lengthPerIt
    real(defReal), intent(in)                           :: rho
    logical(defBool), intent(in)                        :: lin
    integer(shortInt), intent(in)                       :: ani
    logical(defBool), intent(in)                        :: doKinetics
    logical(defBool), intent(in)                        :: loud
    class(dictionary), pointer, intent(inout), optional :: dictFS
    integer(shortInt)                                   :: i
    character(100), parameter :: Here = 'init (arraysRR_class.f90)'

    call self % XSData % init(db, doKinetics, ani, loud)
    self % nG          = self % XSdata % getNG()
    self % geom        => geom
    self % nCells      = self % geom % numberOfCells() 
    
    self % lengthPerIt = lengthPerIt
    self % rho = real(rho, defFlt)
    
    ! Set simulation type
    if (.not. lin .and. ani == 0) then
      self % simulationType = flatIso
    elseif (lin .and. ani == 0) then
      self % simulationType = linearIso
    elseif (.not. lin .and. ani > 0) then
      self % simulationType = flatAni
    else
      self % simulationType = linearAni
    end if
    if (ani >= 0) self % ani = ani

    ! Allocate and initialise arrays
    allocate(self % scalarFlux(self % nG * self % nCells))
    allocate(self % prevFlux(self % nG * self % nCells))
    allocate(self % fluxScores(2, self % nG * self % nCells))
    allocate(self % source(self % nG * self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % lengthSquared(self % nCells))
    allocate(self % volume(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(3, self % nCells))
    
    self % scalarFlux    = 0.0_defFlt
    self % prevFlux      = 1.0_defFlt
    self % fluxScores    = ZERO
    self % source        = 0.0_defFlt
    self % volumeTracks  = ZERO
    self % lengthSquared = ZERO
    self % volume        = ZERO
    self % cellHit       = 0
    self % cellFound     = .false.
    self % cellPos       = -INFINITY

    ! Initialise the fixed source if present
    if (present(dictFS)) then
      allocate(self % fixedSource(self % nG * self % nCells))
      call self % initialiseFixedSource(dictFS)
    end if

    ! TODO: allocate linear and anisotropic components, if present
    if (lin) then

    end if

    if (ani > 0) then

    end if
    
    ! Initialise OMP locks
    allocate(self % locks(self % nCells))
    do i = 1, self % nCells
      call OMP_init_lock(self % locks(i))
    end do

  end subroutine init
  
  !!
  !! Initialises fixed sources to be used in the simulation.
  !! Takes a dictionary containing names of materials in the geometry and
  !! source strengths in each energy group and places these in the appropriate
  !! elements of the fixed source vector.
  !!
  !! Also sets source material identities for future use with uncollided calculations. 
  !!
  subroutine initialiseFixedSource(self, dict)
    class(arraysRR), intent(inout)               :: self
    class(dictionary), intent(inout)             :: dict
    character(nameLen),dimension(:), allocatable :: names
    real(defReal), dimension(:), allocatable     :: sourceStrength
    integer(shortInt)                            :: i, nSource, cIdx
    integer(shortInt), save                      :: g, matIdx, idx, id
    logical(defBool)                             :: found
    character(nameLen)                           :: sourceName
    character(nameLen), save                     :: localName
    character(100), parameter :: Here = 'initialiseFixedSource (arraysRR_class.f90)'
    !$omp threadprivate(matIdx, localName, idx, g, id)

    call dict % keys(names)

    nSource = size(names)

    ! Use for uncollided flux sampling
    allocate(self % sourceIdx(nSource))

    ! Cycle through entries of the dictionary
    do i = 1, nSource

      sourceName = names(i)
      call dict % get(sourceStrength, sourceName)

      ! Ensure correct number of energy groups
      if (size(sourceStrength) /= self % nG) call fatalError(Here,'Source '//sourceName//&
              ' has '//numToChar(size(sourceStrength))//' groups rather than '//numToChar(self % nG))

      ! Make sure that the source corresponds to a material present in the geometry
      found = .false.
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells

        id        = cIdx
        matIdx    = self % geom % geom % graph % getMatFromUID(id)
        localName = self % XSData % getName(matIdx)

        if (localName == sourceName) then

          if (.not. found) then
            !$omp critical
            self % sourceIdx(i) = matIdx         
            !$omp end critical
          end if

          found = .true.
          do g = 1, self % nG
            idx = (cIdx - 1) * self % nG + g
            self % fixedSource(idx) = real(sourceStrength(g),defFlt)
          end do

        end if

      end do
      !$omp end parallel do

      if (.not. found) call fatalError(Here,'The source '//trim(sourceName)//' does not correspond to '//&
              'any material found in the geometry.')

    end do

  end subroutine initialiseFixedSource


  !!
  !! Return a pointer to the flux vector for a given cell
  !!
  subroutine getFluxPointer(self, cIdx, fluxVec)
    class(arraysRR), intent(in), target              :: self
    integer(shortInt), intent(in)                    :: cIdx
    real(defFlt), dimension(:), pointer, intent(out) :: fluxVec
    integer(shortInt)                                :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    fluxVec => self % scalarFlux(baseIdx1:baseIdx2)

  end subroutine getFluxPointer
  
  !!
  !! Return a pointer to the nuclear data object
  !!
  function getDataPointer(self) result(dataPtr)
    class(arraysRR), intent(in), target :: self
    class(dataRR), pointer              :: dataPtr

    dataPtr => self % XSData

  end function getDataPointer
  
  !!
  !! Return a pointer to the geometry object
  !!
  function getGeomPointer(self) result(geomPtr)
    class(arraysRR), intent(in), target :: self
    class(geometryStd), pointer         :: geomPtr

    geomPtr => self % geom

  end function getGeomPointer
  
  !!
  !! Return a pointer to the source vector for a given cell
  !!
  subroutine getSourcePointer(self, cIdx, sourceVec)
    class(arraysRR), intent(in), target              :: self
    integer(shortInt), intent(in)                    :: cIdx
    real(defFlt), dimension(:), pointer, intent(out) :: sourceVec
    integer(shortInt)                                :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    sourceVec => self % source(baseIdx1:baseIdx2)

  end subroutine getSourcePointer
  
  !!
  !! Return source value given cell and group
  !!
  function getSource(self, cIdx, g) result(src)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: src

    src = self % source(self % nG * (cIdx - 1) + g)

  end function getSource
  
  !!
  !! Return previous flux value given cell and group
  !!
  function getPrevFlux(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: flux

    flux = self % prevFlux(self % nG * (cIdx - 1) + g)

  end function getPrevFlux
  
  !!
  !! Return final flux value given cell and group
  !!
  function getFluxScore(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal)                 :: flux

    flux = self % fluxScores(1, self % nG * (cIdx - 1) + g)

  end function getFluxScore
  
  !!
  !! Return final flux standard deviation given cell and group
  !! Will return square of flux scores if called before finaliseFluxScores
  !!
  function getFluxSD(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal)                 :: flux

    flux = self % fluxScores(2, self % nG * (cIdx - 1) + g)

  end function getFluxSD

  !!
  !! Return cell position given cell ID
  !!
  function getCellPos(self, cIdx, g) result(x)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal), dimension(3)   :: x

    x = self % cellPos(1:3, cIdx)

  end function getCellPos
  
  !!
  !! Return the simulation type
  !!
  function getSimulationType(self) result(simType)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: simType

    simType = self % simulationType

  end function getSimulationType

  !!
  !! Increment the local volume estimate in cell cIdx.
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementVolume(self, cIdx, length)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx
    real(defReal), intent(in)      :: length     
    
    self % volumeTracks(cIdx) = self % volumeTracks(cIdx) + length
  
  end subroutine incrementVolume
  
  !!
  !! Increment the sum of length squared in cell cIdx.
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementLengthSquared(self, cIdx, length)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx
    real(defReal), intent(in)      :: length     
    
    self % lengthSquared(cIdx) = self % lengthSquared(cIdx) + length * length
  
  end subroutine incrementLengthSquared
  
  
  !!
  !! Check if a cell has been hit
  !!
  elemental function hasHit(self, cIdx) result (hit)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt)             :: hit
    
    hit = self % cellHit(cIdx)
  
  end function hasHit
  
  !!
  !! Hit a cell 
  !!
  subroutine hitCell(self, cIdx)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx
    
    self % cellHit(cIdx) = 1
  
  end subroutine hitCell

  !!
  !! Return the cell hit rate for the given iteration
  !!
  function getCellHitRate(self) result(hitRate)
    class(arraysRR), intent(in)  :: self
    integer(shortInt)            :: totalHit
    real(defReal)                :: hitRate

    totalHit = sum(self % cellHit)
    hitRate  = real(totalHit, defReal) / self % nCells

  end function getCellHitRate

  !! 
  !! Wipe cell hits
  !!
  subroutine wipeCellHits(self)
    class(arraysRR), intent(inout) :: self

    self % cellHit = 0

  end subroutine wipeCellHits
  
  !!
  !! Has a cell ever been found?
  !!
  function found(self, cIdx) result(wasFound)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    logical(defBool)              :: wasFound
    
    wasFound = self % cellFound(cIdx)
  
  end function found

  !!
  !! Note that a new cell has been found
  !!
  subroutine newFound(self, cIdx, r)
    class(arraysRR), intent(inout)          :: self
    integer(shortInt), intent(in)           :: cIdx
    real(defReal), dimension(3), intent(in) :: r

    !$omp critical 
    self % cellFound(cIdx) = .true.
    self % cellPos(:,cIdx) = r
    !$omp end critical

  end subroutine newFound

  !!
  !! Return number of energy groups used
  !!
  function getNG(self) result(nG)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: nG

    nG = self % nG

  end function getNG

  !!
  !! Set the OMP lock in a given cell
  !!
  subroutine setLock(self, cIdx) 
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx

    call OMP_set_lock(self % locks(cIdx))

  end subroutine setLock
  
  !!
  !! Unset the OMP lock in a given cell
  !!
  subroutine unsetLock(self, cIdx) 
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx

    call OMP_unset_lock(self % locks(cIdx))

  end subroutine unsetLock

  !!
  !! Calls appropriate normalise flux and volume subroutines
  !!
  subroutine normaliseFluxAndVolume(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it
    character(100), parameter :: Here = 'normaliseFluxAndVolume (arraysRR_class.f90)'

    select case(self % simulationType)
      case(flatIso)
        call self % normaliseFluxAndVolumeFlatIso(it)
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

  end subroutine normaliseFluxAndVolume

  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source
  !!
  subroutine normaliseFluxAndVolumeFlatIso(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it
    real(defReal)                             :: norm, normVol
    real(defReal), save                       :: vol
    real(defFlt), save                        :: sigGG, D, norm_V
    real(defFlt), dimension(:), pointer, save :: total
    integer(shortInt), save                   :: g, matIdx, idx
    integer(shortInt)                         :: cIdx
    !$omp threadprivate(total, vol, norm_V, idx, g, matIdx, sigGG, D)

    norm = ONE / self % lengthPerIt
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do 
    do cIdx = 1, self % nCells
      matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
      
      ! Update volume due to additional rays
      self % volume(cIdx) = self % volumeTracks(cIdx) * normVol
      vol = self % volume(cIdx)

      ! Save effort by skipping normalisation if volume is too small
      if (vol < volume_tolerance) then
        do g = 1, self % nG
          idx   = self % nG * (cIdx - 1) + g
          self % scalarFlux(idx) = 0.0_defFlt
        end do
        cycle
      end if
      norm_V = real(norm / vol, defFlt)

      call self % XSData % getTotalPointer(matIdx, total)

      do g = 1, self % nG

        idx   = self % nG * (cIdx - 1) + g
        self % scalarFlux(idx) = self % scalarFlux(idx) * norm_V 

        ! Apply the standard MoC post-sweep treatment and
        ! stabilisation for negative XSs
        if (matIdx <= self % XSData % getNMat() .and. total(g) > 0) then
          
          self % scalarFlux(idx) = self % scalarFlux(idx) / total(g)
          
          sigGG = self % XSData % getScatterXS(matIdx, g, g)

          ! Presumes non-zero total XS
          if ((sigGG < 0) .and. (total(g) > 0)) then
            D = -self % rho * sigGG / total(g)
          else
            D = 0.0_defFlt
          end if
          self % scalarFlux(idx) =  (self % scalarFlux(idx) + self % source(idx)/total(g) &
                + D * self % prevFlux(idx) ) / (1 + D)
        
        ! Alternatively, handle unidentified/void regions
        else
          self % scalarFlux(idx) = self % scalarFlux(idx) + &
               real(self % source(idx) * self % lengthSquared(cIdx) * norm / (2 * vol), defFlt)
        end if
      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolumeFlatIso
  
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for linear isotropic sources
  !!
  subroutine normaliseFluxAndVolumeLinIso(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it

  end subroutine normaliseFluxAndVolumeLinIso
  
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for flat anisotropic sources
  !!
  subroutine normaliseFluxAndVolumeFlatAni(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it

  end subroutine normaliseFluxAndVolumeFlatAni
  
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for Linear sources with flat
  !! anisotropic sources
  !!
  subroutine normaliseFluxAndVolumeLIFA(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it

  end subroutine normaliseFluxAndVolumeLIFA
  
  !!
  !! Update all sources given a prevFlux
  !! This nesting allows using combined OMP + SIMD
  !!
  subroutine updateSource(self, ONE_KEFF)
    class(arraysRR), intent(inout) :: self
    real(defReal), intent(in)      :: ONE_KEFF
    real(defFlt)                   :: ONE_K
    integer(shortInt)              :: cIdx
    character(100), parameter      :: Here = 'updateSource (arraysRR_class.f90)'

    ONE_K = real(ONE_KEFF, defFlt)

    select case(self % simulationType)
      case(flatIso)
        !$omp parallel do 
        do cIdx = 1, self % nCells
          call self % sourceUpdateKernelFlatIso(cIdx, ONE_K)
        end do
        !$omp end parallel do
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

  end subroutine updateSource

  !!
  !! Kernel to update sources given a cell index
  !!
  subroutine sourceUpdateKernelFlatIso(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF
    real(defFlt)                             :: scatter, fission
    real(defFlt), dimension(:), pointer      :: nuFission, chi, scatterXS, fluxVec 
    integer(shortInt)                        :: matIdx, g, gIn, baseIdx, idx, sIdx1, sIdx2

    ! Identify material
    matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
    
    ! Guard against void cells
    if (matIdx > self % XSData % getNMat()) then
      baseIdx = self % nG * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % source(idx) = 0.0_defFlt
      end do
      return
    end if

    ! Obtain XSs
    call self % XSData % getProdPointers(matIdx, nuFission, scatterXS, chi)

    baseIdx = self % nG * (cIdx - 1)
    fluxVec => self % prevFlux((baseIdx + 1):(baseIdx + self % nG))

    ! Calculate fission source
    fission = 0.0_defFlt
    !$omp simd reduction(+:fission)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF

    do g = 1, self % nG

      sIdx1 = self % nG * (g - 1) + 1
      sIdx2 = self % nG * g
      associate(scatterVec => scatterXS(sIdx1:sIdx2))

        ! Calculate scattering source
        scatter = 0.0_defFlt
        !$omp simd reduction(+:scatter)
        do gIn = 1, self % nG
          scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
        end do

      end associate

      ! Output index
      idx = baseIdx + g

      self % source(idx) = chi(g) * fission + scatter
      if (allocated(self % fixedSource)) then
        self % source(idx) = self % source(idx) + self % fixedSource(idx)
      end if

    end do

  end subroutine sourceUpdateKernelFlatIso
  
  !!
  !! Kernel to update sources given a cell index for linear sources
  !! with isotropic scattering
  !!
  subroutine sourceUpdateKernelLinIso(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF

  end subroutine sourceUpdateKernelLinIso
  
  !!
  !! Kernel to update sources given a cell index for flat sources
  !! with anisotropic scattering
  !!
  subroutine sourceUpdateKernelFlatAni(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF

  end subroutine sourceUpdateKernelFlatAni
  
  !!
  !! Kernel to update sources given a cell index for linear sources
  !! with flat anisotropic scattering
  !!
  subroutine sourceUpdateKernelLIFA(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF

  end subroutine sourceUpdateKernelLIFA
  
  
  !! 
  !! Calculate keff
  !! Wraps the main kernel call to allow for OMP + SIMD (thanks Fortran)
  !!
  function calculateKeff(self, k0) result(k1)
    class(arraysRR), intent(in)           :: self
    real(defReal), intent(in)             :: k0
    real(defReal)                         :: k1
    integer(shortInt)                     :: cIdx
    real(defReal)                         :: fissTotal, prevFissTotal
    real(defReal), save                   :: fissLocal, prevFissLocal
    !$omp threadprivate (fissLocal, prevFissLocal)

    fissTotal     = ZERO
    prevFissTotal = ZERO
    !$omp parallel do reduction(+:fissTotal, prevFissTotal)
    do cIdx = 1, self % nCells
      call self % calculateKeffKernel(cIdx, fissLocal, prevFissLocal)
      fissTotal     = fissTotal + fissLocal
      prevFissTotal = prevFissTotal + prevFissLocal
    end do 
    !$omp end parallel do

    k1 = k0 * fissTotal / prevFissTotal 

  end function calculateKeff
  
  !!
  !! Calculate keff for a single cell
  !!
  subroutine calculateKeffKernel(self, cIdx, fissionRate, prevFissionRate)
    class(arraysRR), target, intent(in)  :: self
    integer(shortInt), intent (in)       :: cIdx
    real(defReal), intent(out)           :: fissionRate, prevFissionRate
    real(defReal)                        :: vol
    integer(shortInt)                    :: g, matIdx
    real(defFlt), dimension(:), pointer  :: nuSigmaF, flux, prevFlux

    fissionRate     = ZERO
    prevFissionRate = ZERO

    ! Identify material
    matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
      
    ! Check whether to continue in this cell
    if (matIdx > self % XSData % getNMat()) return
    if (.not. self % XSData % isFissile(matIdx)) return
    vol = self % volume(cIdx)
    if (vol < volume_tolerance) return

    call self % XSData % getNuFissPointer(matIdx, nuSigmaF)
    flux => self % scalarFlux((self % nG * (cIdx - 1) + 1):(self % nG * cIdx))
    prevFlux => self % prevFlux((self % nG * (cIdx - 1) + 1):(self % nG * cIdx))

    !$omp simd reduction(+: fissionRate, prevFissionRate)
    do g = 1, self % nG
      fissionRate     = fissionRate     + real(flux(g) * nuSigmaF(g), defReal)
      prevFissionRate = prevFissionRate + real(prevFlux(g) * nuSigmaF(g), defReal)
    end do

    fissionRate     = fissionRate * vol
    prevFissionRate = prevFissionRate * vol

  end subroutine calculateKeffKernel
  
  !!
  !! Reset fluxes
  !!
  subroutine resetFluxes(self)
    class(arraysRR), intent(inout) :: self
    character(100), parameter      :: Here = 'resetFluxes (arraysRR_class.f90)'

    select case(self % simulationType)
      case(flatIso)
        call self % resetFluxesFlatIso()
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

  end subroutine resetFluxes
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !!
  subroutine resetFluxesFlatIso(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do 
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesFlatIso
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! for linear sources with isotropic scattering
  !!
  subroutine resetFluxesLinIso(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesLinIso

  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! for flat sources with anisotropic scattering
  !!
  subroutine resetFluxesFlatAni(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesFlatAni
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! for linear sources with flat anisotropic scattering
  !!
  subroutine resetFluxesLIFA(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesLIFA


  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxScores(self)
    class(arraysRR), intent(inout) :: self
    real(defReal), save            :: flux
    integer(shortInt)              :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = real(self % scalarFlux(idx),defReal)
      self % fluxScores(1, idx) = self % fluxScores(1, idx) + flux
      self % fluxScores(2, idx) = self % fluxScores(2, idx) + flux * flux
    end do
    !$omp end parallel do

  end subroutine accumulateFluxScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxScores(self,it)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    integer(shortInt)              :: idx
    real(defReal)                  :: N1, Nm1

    if (it /= 1) then
      Nm1 = 1.0_defReal/(it - 1)
    else
      Nm1 = 1.0_defReal
    end if
    N1 = 1.0_defReal/it

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      self % fluxScores(1, idx) = self % fluxScores(1, idx) * N1
      self % fluxScores(2, idx) = self % fluxScores(2, idx) * N1
      self % fluxScores(2, idx) = Nm1 * (self % fluxScores(2, idx) - &
            self % fluxScores(1, idx) * self % fluxScores(1, idx)) 
      if (self % fluxScores(2, idx) <= ZERO) then
        self % fluxScores(2, idx) = ZERO
      else
        self % fluxScores(2, idx) = sqrt(self % fluxScores(2, idx))
      end if
    end do
    !$omp end parallel do

  end subroutine finaliseFluxScores

  !!
  !! Outputs integral flux or fission rate when
  !! given a tally map
  !!
  subroutine outputMap(self, out, map, doFission)
    class(arraysRR), intent(in)                :: self
    class(outputFile), intent(inout)           :: out
    class(tallyMap), pointer, intent(in)       :: map
    logical(defBool), intent(in)               :: doFission
    character(nameLen)                         :: name
    integer(shortInt)                          :: cIdx
    integer(shortInt),dimension(:),allocatable :: resArrayShape
    type(particleState), save                  :: s
    real(defReal), save                        :: vol
    real(defFlt), save                         :: sig
    integer(shortInt), save                    :: i, matIdx, g
    real(defReal), dimension(:), allocatable   :: res, resSD
    !$omp threadprivate(s, vol, sig, i, matIdx, g)

    resArrayShape = [map % binArrayShape()]
    allocate(res(map % bins(0)))
    allocate(resSD(map % bins(0)))
    res   = ZERO
    resSD = ZERO

    ! Find whether cells are in map and sum their contributions
    !$omp parallel do reduction(+: res, resSD)
    do cIdx = 1, self % nCells
        
      vol    =  self % volume(cIdx)
      if (vol < volume_tolerance) cycle

      ! Fudge a particle state to search tally map
      s % r = self % cellPos(:,cIdx)
      i = map % map(s)

      if (i > 0) then
        matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
        do g = 1, self % nG
          if (doFission) then
            sig = self % XSData % getFissionXS(matIdx, g)
          else
            sig = 1.0_defFlt
          end if
          res(i) = res(i) + vol * self % getFluxScore(cIdx,g) * sig
          ! Neglects uncertainty in volume - assumed small.
          resSD(i) = resSD(i) + &
                  self % getFluxSD(cIdx,g)**2 * vol * vol * sig * sig
        end do
      end if

    end do
    !$omp end parallel do

    do i = 1,size(resSD)
      resSD(i) = sqrt(resSD(i))
      if (res(i) > 0) resSD(i) = resSD(i) / res(i)
    end do

    if (doFission) then
      name = 'fiss1G'
    else
      name = 'flux1G'
    end if
    call out % startBlock(name)
    call out % startArray(name, resArrayShape)
    ! Add all map elements to results
    do i = 1, map % bins(0)
      call out % addResult(res(i), resSD(i))
    end do
    call out % endArray()
    ! Output tally map
    call map % print(out)
    call out % endBlock()
      
    deallocate(res)
    deallocate(resSD)

  end subroutine outputMap

  !!
  !! Send all arrays of interest to VTK output
  !!
  subroutine outputToVTK(self, viz)
    class(arraysRR), intent(in)               :: self
    class(visualiser), intent(inout)          :: viz
    real(defReal), dimension(:), allocatable  :: resVec
    character(nameLen)                        :: name
    integer(shortInt)                         :: cIdx, g

    allocate(resVec(self % nCells))

    ! Output all fluxes (assuming finalisation of scores happened)
    do g = 1, self % nG
      name = 'flux_g'//numToChar(g)
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        resVec(cIdx) = self % getFluxScore(cIdx,g)
      end do
      !$omp end parallel do
      call viz % addVTKData(resVec,name)
    end do

    ! Output all flux uncertainties
    do g = 1, self % nG
      name = 'std_g'//numToChar(g)
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        resVec(cIdx) = self % getFluxSD(cIdx,g) /self % getFluxScore(cIdx,g)
      end do
      !$omp end parallel do
      call viz % addVTKData(resVec,name)
    end do

    ! Output final iteration sources
    do g = 1, self % nG
      name = 'source_'//numToChar(g)
      !$omp parallel do schedule(static)
      do cIdx = 1, self % nCells
        resVec(cIdx) = real(self % getSource(cIdx,g),defReal)
      end do
      !$omp end parallel do
      call viz % addVTKData(resVec,name)
    end do

    ! Output final volume estimates
    ! TODO: scale to be absolute, not relative
    name = 'volume'
    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      resVec(cIdx) = self % volume(cIdx)
    end do
    !$omp end parallel do
    call viz % addVTKData(resVec,name)

    call viz % finaliseVTK

  end subroutine outputToVTK
  
  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt) :: i

    ! Clean contents
    if(allocated(self % scalarFlux)) deallocate(self % scalarFlux)
    if(allocated(self % prevFlux)) deallocate(self % prevFlux)
    if(allocated(self % fluxScores)) deallocate(self % fluxScores)
    if(allocated(self % source)) deallocate(self % source)
    if(allocated(self % fixedSource)) deallocate(self % fixedSource)
    if(allocated(self % sourceIdx)) deallocate(self % sourceIdx)
    if(allocated(self % volumeTracks)) deallocate(self % volumeTracks)
    if(allocated(self % lengthSquared)) deallocate(self % lengthSquared)
    if(allocated(self % volume)) deallocate(self % volume)
    if(allocated(self % cellHit)) deallocate(self % cellHit)
    if(allocated(self % cellFound)) deallocate(self % cellFound)
    if(allocated(self % cellPos)) deallocate(self % cellPos)
    
    if(allocated(self % locks)) then
      do i = 1, self % nCells
        call OMP_destroy_lock(self % locks(i))
      end do
      deallocate(self % locks)
    end if
    
    self % geom   => null()
    call self % XSData % kill()
    self % nG     = 0
    self % nCells = 0
    self % lengthPerIt = ZERO
    self % rho         = 0.0_defFlt
    self % simulationType = 0
    self % ani = 0

  end subroutine kill

end module arraysRR_class
