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
  !!   scalarX   -> Array of x-spatial moments of scalar flux [nG * nCells]
  !!   scalarY   -> Array of y-spatial moments of scalar flux [nG * nCells]
  !!   scalarZ   -> Array of z-spatial moments of scalar flux [nG * nCells]
  !!   prevX     -> Array of previous x-spatial moments of scalar flux [nG * nCells]
  !!   prevY     -> Array of previous y-spatial moments of scalar flux [nG * nCells]
  !!   prevZ     -> Array of previous z-spatial moments of scalar flux [nG * nCells]
  !!   sourceX   -> Array of source x-spatial gradients of scalar flux [nG * nCells]
  !!   sourceY   -> Array of source y-spatial gradients of scalar flux [nG * nCells]
  !!   sourceZ   -> Array of source z-spatial gradients of scalar flux [nG * nCells]
  !!   momMat    -> Array of symmetric spatial moment matrices [nCells * matSize]
  !!   momTracks -> Array of weighted tracks used to computer spatial moment matrices [nCells * matSize]
  !!   centroid       -> Array of cell centroid values [nCells * nDim]
  !!   centroidTracks -> Array of weighted tracks used to computer centroids [nCells * nDim]
  !!
  !!   locks -> Array of OpenMP locks for each geometric cell
  !!
  type, public :: arraysRR
    private
    ! Components
    class(geometryStd), pointer :: geom        => null()
    type(dataRR)                :: XSData      
    integer(shortInt)           :: nG          = 0
    integer(shortInt)           :: nCells      = 0
    real(defReal)               :: lengthPerIt = ZERO
    real(defFlt)                :: rho         = 0.0_defFlt
    integer(shortInt)           :: ani         = 0
    integer(shortInt)           :: simulationType = 0
    real(defReal)               :: totalVolume = ONE
    integer(shortInt)           :: volPolicy   = simAverage
    integer(shortInt)           :: missPolicy  = srcPolicy
    
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
    real(defReal), dimension(:), allocatable     :: allVolumeTracks
    real(defReal), dimension(:), allocatable     :: lengthSquared
    real(defReal), dimension(:), allocatable     :: volume
    integer(shortInt), dimension(:), allocatable :: cellHit
    logical(defBool), dimension(:), allocatable  :: cellFound
    real(defReal), dimension(:,:), allocatable   :: cellPos

    ! Linear source arrays
    real(defFlt), dimension(:), allocatable    :: scalarX
    real(defFlt), dimension(:), allocatable    :: scalarY
    real(defFlt), dimension(:), allocatable    :: scalarZ
    real(defFlt), dimension(:), allocatable    :: prevX
    real(defFlt), dimension(:), allocatable    :: prevY
    real(defFlt), dimension(:), allocatable    :: prevZ
    real(defFlt), dimension(:), allocatable    :: sourceX
    real(defFlt), dimension(:), allocatable    :: sourceY
    real(defFlt), dimension(:), allocatable    :: sourceZ
    real(defFlt), dimension(:), allocatable    :: fixedX
    real(defFlt), dimension(:), allocatable    :: fixedY
    real(defFlt), dimension(:), allocatable    :: fixedZ
    real(defReal), dimension(:), allocatable   :: momMat
    real(defReal), dimension(:), allocatable   :: momTracks
    real(defReal), dimension(:), allocatable   :: centroid
    real(defReal), dimension(:), allocatable   :: centroidTracks
    
    real(defReal), dimension(:,:), allocatable :: xScores
    real(defReal), dimension(:,:), allocatable :: yScores
    real(defReal), dimension(:,:), allocatable :: zScores
    
    ! OMP locks
    integer(kind=omp_lock_kind), dimension(:), allocatable :: locks

  contains
    
    ! Public procedures
    procedure :: init
    procedure :: initAdjoint
    procedure :: kill

    ! Access procedures
    procedure :: getDataPointer
    procedure :: getGeomPointer
    procedure :: getFluxPointer
    procedure :: getSourcePointer
    procedure :: getSource
    procedure :: getFixedSource
    procedure :: getPrevFlux
    procedure :: getFluxScore
    procedure :: getFluxSD
    procedure :: getNG
    procedure :: getVolume
    procedure :: getCellPos
    procedure :: wasHit
    procedure :: getCellHitRate
    procedure :: getSimulationType
    procedure :: wasFound
    procedure :: hasFixedSource
    procedure :: getFluxAtAPoint
    
    procedure :: getFluxXYZPointers
    procedure :: getSourceXYZPointers
    procedure :: getCentroid
    procedure :: getMomentMatrix
    procedure :: getFluxMoments
    procedure :: getFluxMomentSDs

    ! Change individual elements of the type
    ! Predominantly for use in the transport sweep
    procedure :: incrementVolume
    procedure :: incrementLengthSquared
    procedure :: incrementCentroid
    procedure :: incrementMoments
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
    procedure :: outputMaterialIntegrals
    procedure :: outputPointFluxes
    
    ! Private procedures
    procedure, private :: initialiseFixedSource
    
    procedure, private :: resetFluxesFlatIso
    procedure, private :: resetFluxesLinearIso
    procedure, private :: resetFluxesLIFA
    procedure, private :: resetFluxesFlatAni
    
    procedure, private :: normaliseFluxAndVolumeFlatIso
    procedure, private :: normaliseFluxAndVolumeLinearIso
    procedure, private :: normaliseFluxAndVolumeLIFA
    procedure, private :: normaliseFluxAndVolumeFlatAni
    
    procedure, private :: sourceUpdateKernelFlatIso
    procedure, private :: sourceUpdateKernelLinearIso
    procedure, private :: sourceUpdateKernelLIFA
    procedure, private :: sourceUpdateKernelFlatAni
    
    procedure, private :: accumulateFluxScoresFlat
    procedure, private :: accumulateFluxScoresLinear
    
    procedure, private :: finaliseFluxScoresFlat
    procedure, private :: finaliseFluxScoresLinear
    
    procedure, private :: calculateKeffKernel
    procedure, private :: invertMatrix
    procedure, private :: materialIntegral

  end type arraysRR

contains

  !!
  !! Initialise the arrays object
  !!
  !! The object is fed sizes and requirements by the physics package.
  !! This will allocate the necessary arrays
  !!
  subroutine init(self, db, geom, lengthPerIt, rho, lin, doKinetics, loud, &
                  dictFS, volPolicy, missPolicy)
    class(arraysRR), intent(inout)                      :: self
    class(baseMgNeutronDatabase), pointer, intent(in)   :: db
    class(geometryStd), pointer, intent(in)             :: geom
    real(defReal), intent(in)                           :: lengthPerIt
    real(defReal), intent(in)                           :: rho
    logical(defBool), intent(in)                        :: lin
    logical(defBool), intent(in)                        :: doKinetics
    logical(defBool), intent(in)                        :: loud
    class(dictionary), pointer, intent(inout), optional :: dictFS
    integer(shortInt), intent(in), optional             :: volPolicy, missPolicy
    integer(shortInt)                                   :: ani, i
    real(defReal), dimension(6)                         :: bb
    character(100), parameter :: Here = 'init (arraysRR_class.f90)'

    call self % XSData % init(db, doKinetics, loud)
    self % nG          = self % XSdata % getNG()
    self % geom        => geom
    self % nCells      = self % geom % numberOfCells()
    
    self % lengthPerIt = lengthPerIt
    self % rho = real(rho, defFlt)

    if (present(volPolicy)) then
      self % volPolicy = volPolicy
    else
      self % volPolicy = simAverage
    end if
    if (present(missPolicy)) then
      self % missPolicy = missPolicy
    else
      self % missPolicy = srcPolicy
    end if

    ! Assume bounding box of the geometry is filled (and a box)
    ! Can this be relaxed in future?
    bb = self % geom % bounds()
    self % totalVolume = (bb(4) - bb(1)) * (bb(5) - bb(2)) * (bb(6) - bb(3))

    ! Set simulation type
    ! TODO: read ani from nuclear data
    ani = 0
    if (.not. lin .and. ani == 0) then
      self % simulationType = flatIso
    elseif (lin .and. ani == 0) then
      self % simulationType = linearIso
    elseif (.not. lin .and. ani > 0) then
      self % simulationType = flatAni
    else
      self % simulationType = linearAni
    end if
    self % ani = ani

    ! Allocate and initialise arrays
    allocate(self % scalarFlux(self % nG * self % nCells))
    allocate(self % prevFlux(self % nG * self % nCells))
    allocate(self % fluxScores(2, self % nG * self % nCells))
    allocate(self % source(self % nG * self % nCells))
    allocate(self % volumeTracks(self % nCells))
    allocate(self % allVolumeTracks(self % nCells))
    allocate(self % lengthSquared(self % nCells))
    allocate(self % volume(self % nCells))
    allocate(self % cellHit(self % nCells))
    allocate(self % cellFound(self % nCells))
    allocate(self % cellPos(nDim, self % nCells))
    
    self % scalarFlux    = 0.0_defFlt
    self % prevFlux      = 1.0_defFlt
    self % fluxScores    = ZERO
    self % source        = 0.0_defFlt
    self % volumeTracks  = ZERO
    self % allVolumeTracks = ZERO
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

    ! Allocate linear components, if present
    if (lin) then
    
      allocate(self % scalarX(self % nCells * self % nG))
      allocate(self % scalarY(self % nCells * self % nG))
      allocate(self % scalarZ(self % nCells * self % nG))
      allocate(self % prevX(self % nCells * self % nG))
      allocate(self % prevY(self % nCells * self % nG))
      allocate(self % prevZ(self % nCells * self % nG))
      allocate(self % sourceX(self % nCells * self % nG))
      allocate(self % sourceY(self % nCells * self % nG))
      allocate(self % sourceZ(self % nCells * self % nG))
      allocate(self % momMat(self % nCells * matSize))
      allocate(self % momTracks(self % nCells * matSize))
      allocate(self % centroid(self % nCells * nDim))
      allocate(self % centroidTracks(self % nCells * nDim))
      allocate(self % xScores(2, self % nG * self % nCells))
      allocate(self % yScores(2, self % nG * self % nCells))
      allocate(self % zScores(2, self % nG * self % nCells))
      
      self % scalarX        = 0.0_defFlt
      self % scalarY        = 0.0_defFlt
      self % scalarZ        = 0.0_defFlt
      self % prevX          = 0.0_defFlt
      self % prevY          = 0.0_defFlt
      self % prevZ          = 0.0_defFlt
      self % sourceX        = 0.0_defFlt
      self % sourceY        = 0.0_defFlt
      self % sourceZ        = 0.0_defFlt
      self % momMat         = ZERO
      self % momTracks      = ZERO
      self % centroid       = ZERO
      self % centroidTracks = ZERO
      self % xScores        = ZERO
      self % yScores        = ZERO
      self % zScores        = ZERO

    end if

    ! TODO: allocate anisotropic components, if present
    if (ani > 0) then

    end if
    
    ! Initialise OMP locks
    allocate(self % locks(self % nCells))
    do i = 1, self % nCells
#ifdef _OPENMP
      call OMP_init_lock(self % locks(i))
#endif
    end do

  end subroutine init
  
  !!
  !! Initialise the adjoint source and update nuclear data.
  !! For now, assumes the adjoint is for global variance reduction.
  !!
  subroutine initAdjoint(self)
    class(arraysRR), intent(inout)         :: self
    integer(shortInt)                      :: cIdx
    logical(defBool)                       :: doLinear
    integer(shortInt), save                :: i, g
    real(defFlt), dimension(matSize), save :: invM
    real(defFlt), save                     :: xMom, yMom, zMom
    !$omp threadprivate(i, g, invM, xMom, yMom, zMom)

    doLinear = .false.
    invM = 0.0_defFlt

    call self % xsData % setAdjointXS()

    if (.not. allocated(self % fixedSource)) then
      allocate(self % fixedSource(size(self % scalarFlux)))
    end if
    self % fixedSource = 0.0_defFlt

    if (allocated(self % scalarX)) then
      doLinear = .true.
      if (.not. allocated(self % fixedX)) then
        allocate(self % fixedX(size(self % scalarFlux)))
        allocate(self % fixedY(size(self % scalarFlux)))
        allocate(self % fixedZ(size(self % scalarFlux)))
      end if
      self % fixedX = 0.0_defFlt
      self % fixedY = 0.0_defFlt
      self % fixedZ = 0.0_defFlt
    end if 

    ! Create fixed source from the flux scores
    ! Presently assumes the response of interest is global flux
    !$omp parallel do
    do cIdx = 1, self % nCells

      if (.not. self % wasFound(cIdx)) cycle
      if (doLinear) invM = self % invertMatrix(cIdx)

      do g = 1, self % nG

        i = (cIdx - 1) * self % nG + g

        ! Check for inordinately small flux values.
        ! Note, these can have arbitrarily low magnitude.
        ! Maybe should be something more robust.
        if (self % fluxScores(1, i) == 0.0_defFlt) cycle
        self % fixedSource(i) = real(ONE / self % fluxScores(1, i), defFlt)

        ! Linear source treatment relies on performing a Taylor expansion
        ! of q' = 1/phi = 1/(phi_0 + grad phi * (r - r0)) 
        ! = 1/phi_0 - 1/phi^2_0 * gradPhi * (r - r0)
        !if (doLinear) then
        if (1 == 0) then
          xMom = real(self % xScores(1, i), defFlt)
          yMom = real(self % yScores(1, i), defFlt)
          zMom = real(self % zScores(1, i), defFlt)
          self % fixedX(i) = invM(xx) * xMom + invM(xy) * yMom + invM(xz) * zMom 
          self % fixedY(i) = invM(xy) * xMom + invM(yy) * yMom + invM(yz) * zMom 
          self % fixedZ(i) = invM(xz) * xMom + invM(yz) * yMom + invM(zz) * zMom 

          self % fixedX(i) = -self % fixedX(i) * self % fixedSource(i) ** 2
          self % fixedY(i) = -self % fixedY(i) * self % fixedSource(i) ** 2
          self % fixedZ(i) = -self % fixedZ(i) * self % fixedSource(i) ** 2
        end if

      end do

    end do 
    !$omp end parallel do
    
    ! Reinitialise arrays to be used during transport
    self % scalarFlux      = 0.0_defFlt
    self % prevFlux        = 0.0_defFlt
    self % fluxScores      = ZERO
    self % source          = 0.0_defFlt
    
    ! Ideally we would have a way of reusing the volume estimators
    ! No compact ideas at the moment, so these will simply be reinitialised
    self % volumeTracks    = ZERO
    self % allVolumeTracks = ZERO
    self % lengthSquared   = ZERO
    self % volume          = ZERO

    if (doLinear) then
      self % scalarX        = 0.0_defFlt
      self % scalarY        = 0.0_defFlt
      self % scalarZ        = 0.0_defFlt
      self % prevX          = 0.0_defFlt
      self % prevY          = 0.0_defFlt
      self % prevZ          = 0.0_defFlt
      self % sourceX        = 0.0_defFlt
      self % sourceY        = 0.0_defFlt
      self % sourceZ        = 0.0_defFlt
      self % momMat         = ZERO
      self % momTracks      = ZERO
      self % centroid       = ZERO
      self % centroidTracks = ZERO
      self % xScores        = ZERO
      self % yScores        = ZERO
      self % zScores        = ZERO
    end if

  end subroutine initAdjoint

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
  !! Return a pointer to the flux spatial moment vectors for a given cell
  !!
  subroutine getFluxXYZPointers(self, cIdx, xVec, yVec, zVec)
    class(arraysRR), intent(in), target              :: self
    integer(shortInt), intent(in)                    :: cIdx
    real(defFlt), dimension(:), pointer, intent(out) :: xVec
    real(defFlt), dimension(:), pointer, intent(out) :: yVec
    real(defFlt), dimension(:), pointer, intent(out) :: zVec
    integer(shortInt)                                :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    xVec => self % scalarX(baseIdx1:baseIdx2)
    yVec => self % scalarY(baseIdx1:baseIdx2)
    zVec => self % scalarZ(baseIdx1:baseIdx2)

  end subroutine getFluxXYZPointers
  
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
  !! Return pointers to the source gradient vectors for a given cell
  !!
  subroutine getSourceXYZPointers(self, cIdx, xVec, yVec, zVec)
    class(arraysRR), intent(in), target              :: self
    integer(shortInt), intent(in)                    :: cIdx
    real(defFlt), dimension(:), pointer, intent(out) :: xVec
    real(defFlt), dimension(:), pointer, intent(out) :: yVec
    real(defFlt), dimension(:), pointer, intent(out) :: zVec
    integer(shortInt)                                :: baseIdx1, baseIdx2

    baseIdx1 = self % nG * (cIdx - 1) + 1
    baseIdx2 = self % nG * cIdx
    xVec => self % sourceX(baseIdx1:baseIdx2)
    yVec => self % sourceY(baseIdx1:baseIdx2)
    zVec => self % sourceZ(baseIdx1:baseIdx2)

  end subroutine getSourceXYZPointers
  
  !!
  !! Return source value given cell and group
  !!
  elemental function getSource(self, cIdx, g) result(src)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: src

    src = self % source(self % nG * (cIdx - 1) + g)

  end function getSource
  
  !!
  !! Return fixed source value given cell and group
  !!
  elemental function getFixedSource(self, cIdx, g) result(src)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: src

    if (allocated(self % fixedSource)) then
      src = self % fixedSource(self % nG * (cIdx - 1) + g)
    else
      src = 0.0_defFlt
    end if

  end function getFixedSource
  
  
  !!
  !! Return previous flux value given cell and group
  !!
  elemental function getPrevFlux(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defFlt)                  :: flux

    flux = self % prevFlux(self % nG * (cIdx - 1) + g)

  end function getPrevFlux
  
  !!
  !! Return final flux value given cell and group
  !!
  elemental function getFluxScore(self, cIdx, g) result(flux)
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
  elemental function getFluxSD(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal)                 :: flux

    flux = self % fluxScores(2, self % nG * (cIdx - 1) + g)

  end function getFluxSD
  
  !!
  !! Return final flux moment values given cell and group
  !!
  pure function getFluxMoments(self, cIdx, g) result(flux)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal), dimension(3)   :: flux
    integer(shortInt)             :: idx

    idx = self % nG * (cIdx - 1) + g
    flux(1) = self % xScores(1, idx)
    flux(2) = self % yScores(1, idx)
    flux(3) = self % zScores(1, idx)

  end function getFluxMoments
  
  !!
  !! Return final flux moment standard deviations given cell and group
  !! Will return square of moment scores if called before finaliseFluxScores
  !!
  pure function getFluxMomentSDs(self, cIdx, g) result(fluxSD)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    integer(shortInt), intent(in) :: g
    real(defReal), dimension(3)   :: fluxSD
    integer(shortInt)             :: idx

    idx = self % nG * (cIdx - 1) + g
    fluxSD(1) = self % xScores(2, idx)
    fluxSD(2) = self % yScores(2, idx)
    fluxSD(3) = self % zScores(2, idx)

  end function getFluxMomentSDs
  
  !!
  !! Return volume given cell ID
  !!
  elemental function getVolume(self, cIdx) result(vol)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    real(defReal)                 :: vol

    vol = self % volume(cIdx)

  end function getVolume

  !!
  !! Return cell position given cell ID
  !!
  pure function getCellPos(self, cIdx) result(pos)
    class(arraysRR), intent(in)    :: self
    integer(shortInt), intent(in)  :: cIdx
    real(defReal), dimension(nDim) :: pos

    pos = self % cellPos(1:nDim, cIdx)

  end function getCellPos
  
  !!
  !! Return cell centroid given cell ID
  !!
  pure function getCentroid(self, cIdx) result(cent)
    class(arraysRR), intent(in)    :: self
    integer(shortInt), intent(in)  :: cIdx
    real(defReal), dimension(nDim) :: cent
    integer(shortInt)              :: idx0, idx1

    idx0 = nDim * (cIdx - 1) + 1
    idx1 = nDim * cIdx
    cent = self % centroid(idx0:idx1)

  end function getCentroid
  
  !!
  !! Return moment matrix given cell ID
  !!
  pure function getMomentMatrix(self, cIdx) result(mat)
    class(arraysRR), intent(in)       :: self
    integer(shortInt), intent(in)     :: cIdx
    real(defReal), dimension(matSize) :: mat
    integer(shortInt)              :: idx0, idx1

    idx0 = matSize * (cIdx - 1) + 1
    idx1 = matSize * cIdx
    mat = self % momMat(idx0:idx1)

  end function getMomentMatrix
  
  !!
  !! Return the simulation type
  !!
  elemental function getSimulationType(self) result(simType)
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
  !! Increment the local centroid estimate in cell cIdx.
  !! rL is the tracklength-weighted centroid
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementCentroid(self, cIdx, rL)
    class(arraysRR), intent(inout)             :: self
    integer(shortInt), intent(in)              :: cIdx
    real(defReal), dimension(nDim), intent(in) :: rL
    integer(shortInt)                          :: idx0, idx1

    idx0 = nDim * (cIdx - 1) + 1
    idx1 = nDim * cIdx 
    self % centroidTracks(idx0:idx1) = self % centroidTracks(idx0:idx1) + rL

  end subroutine incrementCentroid
  
  !!
  !! Increment the local moment matrix estimate in cell cIdx.
  !! mat is the tracklength-weighted matrix
  !! Assumes this is being called inside a lock for thread privacy.
  !!
  subroutine incrementMoments(self, cIdx, mat)
    class(arraysRR), intent(inout)                :: self
    integer(shortInt), intent(in)                 :: cIdx
    real(defReal), dimension(matSize), intent(in) :: mat
    integer(shortInt)                             :: idx0, idx1

    idx0 = matSize * (cIdx - 1) + 1
    idx1 = matSize * cIdx 
    self % momTracks(idx0:idx1) = self % momTracks(idx0:idx1) + mat
  
  end subroutine incrementMoments
 
  !!
  !! Check if a cell has been hit
  !!
  elemental function wasHit(self, cIdx) result (hit)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    logical(defBool)              :: hit
    
    hit = (self % cellHit(cIdx) == 1)
  
  end function wasHit
  
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
  elemental function getCellHitRate(self, it) result(hitRate)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: it
    integer(shortInt)             :: totalHit, realCells
    real(defReal)                 :: hitRate

    totalHit = sum(self % cellHit)
    if (it > 20) then 
      realCells = count(self % cellFound)
    else
      realCells = self % nCells
    end if
    hitRate = real(totalHit,defReal) / realCells 
    !print *, realCells

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
  elemental function wasFound(self, cIdx) result(found)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    logical(defBool)              :: found
    
    found = self % cellFound(cIdx)
  
  end function wasFound

  !!
  !! Note that a new cell has been found
  !!
  subroutine newFound(self, cIdx, r)
    class(arraysRR), intent(inout)          :: self
    integer(shortInt), intent(in)           :: cIdx
    real(defReal), dimension(3), intent(in) :: r

    ! Remove critical if this is to go in the lock
    !$omp critical 
    self % cellFound(cIdx) = .true.
    self % cellPos(:,cIdx) = r
    !$omp end critical

  end subroutine newFound

  !!
  !! Return number of energy groups used
  !!
  elemental function getNG(self) result(nG)
    class(arraysRR), intent(in) :: self
    integer(shortInt)           :: nG

    nG = self % nG

  end function getNG
  
  !!
  !! Check if a cell has an inhomogeneous source
  !!
  elemental function hasFixedSource(self, cIdx) result (hasSrc)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: cIdx
    logical(defBool)              :: hasSrc
    integer(shortInt)             :: idx1, idx2
    
    if (allocated(self % fixedSource)) then
      idx1 = self % nG * (cIdx - 1) + 1
      idx2 = self % nG * cIdx
      ! Take an absolute value in case of (possibly desirable?) negative sources
      hasSrc = any(abs(self % fixedSource(idx1:idx2)) > 0.0_defFlt)
    else
      hasSrc = .false.
    end if
  
  end function hasFixedSource

  !!
  !! Set the OMP lock in a given cell
  !!
  subroutine setLock(self, cIdx) 
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx

#ifdef _OPENMP
    call OMP_set_lock(self % locks(cIdx))
#endif

  end subroutine setLock
  
  !!
  !! Unset the OMP lock in a given cell
  !!
  subroutine unsetLock(self, cIdx) 
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: cIdx

#ifdef _OPENMP
    call OMP_unset_lock(self % locks(cIdx))
#endif

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
      case(linearIso)
        call self % normaliseFluxAndVolumeLinearIso(it)
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
    real(defReal), save                       :: vol, volAve, volNaive
    real(defFlt), save                        :: sigGG, D, norm_V
    real(defFlt), dimension(:), pointer, save :: total
    integer(shortInt), save                   :: g, matIdx, idx
    integer(shortInt)                         :: cIdx
    logical(defBool), save                    :: hit, isSrc
    character(100), parameter :: Here = 'normaliseFluxAndVolumeLinearIso (arraysRR_class.f90)'
    !$omp threadprivate(total, vol, norm_V, idx, g, matIdx, sigGG, D, hit, isSrc, volAve, volNaive)

    norm = ONE / self % lengthPerIt
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do 
    do cIdx = 1, self % nCells
      matIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
      
      hit = self % wasHit(cIdx)
      isSrc = self % hasFixedSource(cIdx)
      
      self % allVolumeTracks(cIdx) = self % allVolumeTracks(cIdx) + &
              self % volumeTracks(cIdx)
      self % volume(cIdx) = self % allVolumeTracks(cIdx) * normVol
      volAve = self % volume(cIdx)
      volNaive = self % volumeTracks(cIdx) * norm
      self % volumeTracks(cIdx) = ZERO

      ! Decide volume to use
      select case(self % volPolicy)
        case(simAverage)
          vol = volAve
        case(naive)
          vol = volNaive
        case(hybrid)
          if (isSrc) then
            vol = volNaive
          else
            vol = volAve
          end if
        case default
          call fatalError(Here,'Unsupported volume handling requested')
      end select
      norm_V = real(norm / vol, defFlt)
      
      call self % XSData % getTotalPointer(matIdx, total)

      do g = 1, self % nG

        idx   = self % nG * (cIdx - 1) + g

        if (hit) then

          ! Can hit a cell but with a tiny volume, such that 
          ! things break a bit - would rather remove this arbitrary
          ! check in future
          if (vol < volume_tolerance) then
            self % scalarFlux(idx) = 0.0_defFlt
            cycle 
          end if
          
          self % scalarFlux(idx) = self % scalarFlux(idx) * norm_V 
          self % scalarFlux(idx) = self % scalarFlux(idx) / total(g)
          
          ! Presumes non-zero total XS
          sigGG = self % XSData % getScatterXS(matIdx, g, g)
          if ((sigGG < 0) .and. (total(g) > 0)) then
            D = -self % rho * sigGG / total(g)
          else
            D = 0.0_defFlt
          end if
          self % scalarFlux(idx) =  (self % scalarFlux(idx) + self % source(idx)/total(g) &
                + D * self % prevFlux(idx) ) / (1 + D)
        else
          ! Decide flux treatment to use
          select case(self % missPolicy)
            case(srcPolicy)
              self % scalarFlux(idx) = self % source(idx) / total(g)
            case(prevPolicy)
              self % scalarFlux(idx) = self % prevFlux(idx)
            case(hybrid)
              if (isSrc) then
                self % scalarFlux(idx) = self % prevFlux(idx)
              else
                self % scalarFlux(idx) = self % source(idx) / total(g)
              end if
            case default
              call fatalError(Here,'Unsupported miss handling requested')
          end select

        end if
        ! Alternatively, handle unidentified/void regions
        !else
        !  self % scalarFlux(idx) = self % scalarFlux(idx) * norm_V 
        !  self % scalarFlux(idx) = self % scalarFlux(idx) + &
        !       real(self % source(idx) * self % lengthSquared(cIdx) * norm / (2 * vol), defFlt)
        !end if
      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolumeFlatIso
  
  !!
  !! Normalise flux and volume by total track length and increments
  !! the flux by the neutron source for linear isotropic sources
  !!
  subroutine normaliseFluxAndVolumeLinearIso(self, it)
    class(arraysRR), intent(inout)            :: self
    integer(shortInt), intent(in)             :: it
    real(defReal)                             :: norm, normVol
    real(defReal), save                       :: vol, invVol, volNaive, volAve
    real(defFlt), save                        :: norm_V, D, sigGG
    real(defFlt), dimension(:), pointer, save :: total
    integer(shortInt)                         :: cIdx
    integer(shortInt), save                   :: g, matIdx, idx, dIdx, mIdx
    logical(defBool), save                    :: hit, isSrc
    character(100), parameter :: Here = 'normaliseFluxAndVolumeLinearIso (arraysRR_class.f90)'
    !$omp threadprivate(total, vol, idx, mIdx, dIdx, g, matIdx, invVol, norm_V, D, sigGG, hit, isSrc, volNaive, volAve)

    norm = ONE / self % lengthPerIt
    normVol = ONE / (self % lengthPerIt * it)

    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      matIdx =  self % geom % geom % graph % getMatFromUID(cIdx)
      dIdx = (cIdx - 1) * nDim
      mIdx = (cIdx - 1) * matSize

      hit = self % wasHit(cIdx)
      isSrc = self % hasFixedSource(cIdx)

      self % allVolumeTracks(cIdx) = self % allVolumeTracks(cIdx) + &
              self % volumeTracks(cIdx)
      self % volume(cIdx) = self % allVolumeTracks(cIdx) * normVol
      volAve = self % volume(cIdx)
      volNaive = self % volumeTracks(cIdx) * norm
      self % volumeTracks(cIdx) = ZERO

      ! Decide volume to use
      select case(self % volPolicy)
        case(naive)
          vol = volNaive
        case(simAverage)
          vol = volAve
        case(hybrid)
          if (isSrc) then
            vol = volNaive
          else
            vol = volAve
          end if
        case default
          call fatalError(Here,'Unsupported volume handling requested')
      end select
      invVol = ONE / self % allVolumeTracks(cIdx)
       
      ! Update geometric information provided volume has been visited
      if (self % allVolumeTracks(cIdx) > ZERO) then
        ! Update centroids
        self % centroid(dIdx + x) =  self % centroidTracks(dIdx + x) * invVol
        self % centroid(dIdx + y) =  self % centroidTracks(dIdx + y) * invVol
        self % centroid(dIdx + z) =  self % centroidTracks(dIdx + z) * invVol
      
        ! Update spatial moments
        self % momMat(mIdx + xx) = self % momTracks(mIdx + xx) * invVol
        self % momMat(mIdx + xy) = self % momTracks(mIdx + xy) * invVol
        self % momMat(mIdx + xz) = self % momTracks(mIdx + xz) * invVol
        self % momMat(mIdx + yy) = self % momTracks(mIdx + yy) * invVol
        self % momMat(mIdx + yz) = self % momTracks(mIdx + yz) * invVol
        self % momMat(mIdx + zz) = self % momTracks(mIdx + zz) * invVol

      end if

      call self % XSData % getTotalPointer(matIdx, total)
      norm_V = real(norm / vol, defFlt)

      do g = 1, self % nG

        idx = self % nG * (cIdx - 1) + g
        
        if (hit) then
      
          ! Can hit a cell but with a tiny volume, such that 
          ! things break a bit - would rather remove this arbitrary
          ! check in future
          if (vol < volume_tolerance) then
            self % scalarFlux(idx) = 0.0_defFlt
            self % scalarX(idx) = 0.0_defFlt
            self % scalarY(idx) = 0.0_defFlt
            self % scalarZ(idx) = 0.0_defFlt
            cycle
          end if
          
          self % scalarFlux(idx) = self % scalarFlux(idx) * norm_V
          self % scalarX(idx) = self % scalarX(idx) * norm_V
          self % scalarY(idx) = self % scalarY(idx) * norm_V
          self % scalarZ(idx) = self % scalarZ(idx) * norm_V 

          ! Apply the standard MoC post-sweep treatment and
          ! stabilisation for negative XSs
          !if (matIdx <= self % XSData % getNMat() .and. total(g) > 0) then
          
            self % scalarFlux(idx) = self % scalarFlux(idx) / total(g)
            self % scalarX(idx) = self % scalarX(idx) / total(g)
            self % scalarY(idx) = self % scalarY(idx) / total(g)
            self % scalarZ(idx) = self % scalarZ(idx) / total(g)
          
            ! Presumes non-zero total XS
            sigGG = self % XSData % getScatterXS(matIdx, g, g)
            if (sigGG < 0) then
              D = -self % rho * sigGG / total(g)
            else
              D = 0.0_defFlt
            end if
           
            self % scalarFlux(idx) =  (self % scalarFlux(idx) + self % source(idx)/total(g) &
                  + D * self % prevFlux(idx) ) / (1 + D)
          !end if
        else
          ! Decide flux treatment to use
          associate(mat => self % momMat((mIdx + 1):(mIdx + matSize)))
          select case(self % missPolicy)
            ! Note: this is policy to use the source, not policy for hitting a fixed source
            case(srcPolicy)
              self % scalarFlux(idx) = self % source(idx) / total(g)
              ! OPENMC SETS MOMENTS TO ZERO
              !self % scalarX(idx) = 0.0_defFlt
              !self % scalarY(idx) = 0.0_defFlt
              !self % scalarZ(idx) = 0.0_defFlt
              ! Need to multiply source gradients by moment matrix
              self % scalarX(idx) = real(mat(xx) *self % sourceX(idx) + &
                      mat(xy) * self % sourceY(idx) + mat(xz) * self % sourceZ(idx),defFlt)/ total(g)
              self % scalarY(idx) = real(mat(xy) *self % sourceX(idx) + &
                      mat(yy) * self % sourceY(idx) + mat(yz) * self % sourceZ(idx),defFlt)/ total(g)
              self % scalarZ(idx) = real(mat(xz) *self % sourceX(idx) + &
                      mat(yz) * self % sourceY(idx) + mat(zz) * self % sourceZ(idx),defFlt)/ total(g)
            case(prevPolicy)
              self % scalarFlux(idx) = self % prevFlux(idx)
              self % scalarX(idx) = self % prevX(idx)
              self % scalarY(idx) = self % prevY(idx)
              self % scalarZ(idx) = self % prevZ(idx)
            case(hybrid)
              if (isSrc) then
                self % scalarFlux(idx) = self % prevFlux(idx)
                self % scalarX(idx) = self % prevX(idx)
                self % scalarY(idx) = self % prevY(idx)
                self % scalarZ(idx) = self % prevZ(idx)
              else
                self % scalarFlux(idx) = self % source(idx) / total(g)
                ! OPENMC SETS MOMENTS TO ZERO
                !self % scalarX(idx) = 0.0_defFlt
                !self % scalarY(idx) = 0.0_defFlt
                !self % scalarZ(idx) = 0.0_defFlt
                ! Need to multiply source gradients by moment matrix
                self % scalarX(idx) = real(mat(xx) *self % sourceX(idx) + &
                      mat(xy) * self % sourceY(idx) + mat(xz) * self % sourceZ(idx),defFlt)/ total(g)
                self % scalarY(idx) = real(mat(xy) *self % sourceX(idx) + &
                      mat(yy) * self % sourceY(idx) + mat(yz) * self % sourceZ(idx),defFlt)/ total(g)
                self % scalarZ(idx) = real(mat(xz) *self % sourceX(idx) + &
                      mat(yz) * self % sourceY(idx) + mat(zz) * self % sourceZ(idx),defFlt)/ total(g)
              end if
            case default
              call fatalError(Here,'Unsupported miss handling requested')
          end select
          end associate
        end if
        
        if (it < 10) then
          self % scalarX(idx) = 0.0_defFlt
          self % scalarY(idx) = 0.0_defFlt
          self % scalarZ(idx) = 0.0_defFlt
        end if

      end do

    end do
    !$omp end parallel do

  end subroutine normaliseFluxAndVolumeLinearIso
  
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
      case(linearIso)
        !$omp parallel do 
        do cIdx = 1, self % nCells
          call self % sourceUpdateKernelLinearIso(cIdx, ONE_K)
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
        if (allocated(self % fixedSource)) then
          self % source(idx) = self % source(idx) + self % fixedSource(idx)
        end if
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
  subroutine sourceUpdateKernelLinearIso(self, cIdx, ONE_KEFF)
    class(arraysRR), target, intent(inout)   :: self
    integer(shortInt), intent(in)            :: cIdx
    real(defFlt), intent(in)                 :: ONE_KEFF
    real(defFlt)                             :: scatter, xScatter, yScatter, zScatter, &
                                                fission, xFission, yFission, zFission, &
                                                xSource, ySource, zSource
    real(defFlt), dimension(matSize)         :: invM
    real(defFlt), dimension(:), pointer      :: nuFission, chi, scatterXS 
    integer(shortInt)                        :: matIdx, g, gIn, baseIdx, idx, sIdx1, sIdx2
    real(defFlt), pointer, dimension(:)      :: fluxVec, xFluxVec, yFluxVec, zFluxVec

    ! Invert moment matrix
    invM = self % invertMatrix(cIdx)
    
    ! Identify material
    matIdx = self % geom % geom % graph % getMatFromUID(cIdx)
    
    ! Guard against void cells
    if (matIdx > self % XSData % getNMat()) then
      baseIdx = self % nG * (cIdx - 1)
      do g = 1, self % nG
        idx = baseIdx + g
        self % source(idx) = 0.0_defFlt
        if (allocated(self % fixedSource)) then
          self % source(idx) = self % source(idx) + self % fixedSource(idx)
        end if
        self % sourceX(idx) = 0.0_defFlt
        self % sourceY(idx) = 0.0_defFlt
        self % sourceZ(idx) = 0.0_defFlt
      end do
      return
    end if
     
    ! Obtain XSs
    call self % XSData % getProdPointers(matIdx, nuFission, scatterXS, chi)

    baseIdx = self % nG * (cIdx - 1)
    fluxVec  => self % prevFlux((baseIdx + 1):(baseIdx + self % nG))
    xFluxVec => self % prevX((baseIdx + 1):(baseIdx + self % nG))
    yFluxVec => self % prevY((baseIdx + 1):(baseIdx + self % nG))
    zFluxVec => self % prevZ((baseIdx + 1):(baseIdx + self % nG))
    
    ! Calculate fission source
    fission = 0.0_defFlt
    xFission = 0.0_defFlt
    yFission = 0.0_defFlt
    zFission = 0.0_defFlt

    !$omp simd reduction(+:fission, xFission, yFission, zFission)
    do gIn = 1, self % nG
      fission = fission + fluxVec(gIn) * nuFission(gIn)
      xFission = xFission + xFluxVec(gIn) * nuFission(gIn)
      yFission = yFission + yFluxVec(gIn) * nuFission(gIn)
      zFission = zFission + zFluxVec(gIn) * nuFission(gIn)
    end do
    fission = fission * ONE_KEFF
    xFission = xFission * ONE_KEFF
    yFission = yFission * ONE_KEFF
    zFission = zFission * ONE_KEFF

    do g = 1, self % nG

      sIdx1 = self % nG * (g - 1) + 1
      sIdx2 = self % nG * g
      associate(scatterVec => scatterXS(sIdx1:sIdx2))

        ! Calculate scattering source
        scatter = 0.0_defFlt
        xScatter = 0.0_defFlt
        yScatter = 0.0_defFlt
        zScatter = 0.0_defFlt
        !$omp simd reduction(+:scatter, xScatter, yScatter, zScatter)
        do gIn = 1, self % nG
          scatter = scatter + fluxVec(gIn) * scatterVec(gIn)
          xScatter = xScatter + xFluxVec(gIn) * scatterVec(gIn)
          yScatter = yScatter + yFluxVec(gIn) * scatterVec(gIn)
          zScatter = zScatter + zFluxVec(gIn) * scatterVec(gIn)
        end do

      end associate

      ! Output index
      idx = baseIdx + g

      self % source(idx) = chi(g) * fission + scatter
      if (allocated(self % fixedSource)) then
        self % source(idx) = self % source(idx) + self % fixedSource(idx)
      end if
      xSource = chi(g) * xFission + xScatter
      ySource = chi(g) * yFission + yScatter
      zSource = chi(g) * zFission + zScatter
        
      ! Calculate source gradients by inverting the moment matrix
      self % sourceX(idx) = invM(xx) * xSource + &
              invM(xy) * ySource + invM(xz) * zSource
      self % sourceY(idx) = invM(xy) * xSource + &
              invM(yy) * ySource + invM(yz) * zSource
      self % sourceZ(idx) = invM(xz) * xSource + &
           invM(yz) * ySource + invM(zz) * zSource
      if (allocated(self % fixedX)) then
        self % sourceX(idx) = self % sourceX(idx) + self % fixedX(idx)
        self % sourceY(idx) = self % sourceY(idx) + self % fixedY(idx)
        self % sourceZ(idx) = self % sourceZ(idx) + self % fixedZ(idx)
      end if

    end do

  end subroutine sourceUpdateKernelLinearIso
  
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
  !! Inverts the spatial moment matrix for use in linear source calculations.
  !!
  function invertMatrix(self, cIdx) result(invM)
    class(arraysRR), target, intent(in) :: self
    integer(shortInt), intent(in)       :: cIdx
    real(defFlt), dimension(matSize)    :: invM
    integer(shortInt)                   :: condX, condY, condZ, inversionTest
    real(defReal)                       :: det
    real(defFlt)                        :: one_det  

    associate(momVec => self % momMat(((cIdx - 1) * matSize + 1):(cIdx * matSize)))

    ! Pre-invert the moment matrix
    ! Need to check for poor conditioning by evaluating the
    ! diagonal elements of the matrix
    condX = 0
    condY = 0
    condZ = 0

    ! Trying out simpler matrix test
    if (momVec(xx) > condition_tolerance) condX = 1
    if (momVec(yy) > condition_tolerance) condY = 1
    if (momVec(zz) > condition_tolerance) condZ = 1
    !condX = 1
    !condY = 1
    !DEFINITELY DO THIS FOR C5G7 WITH LOW POP
    !condZ = 0

    ! Map conditions to test variable
    inversionTest = condX * 4 + condY * 2 + condZ
    invM = 0.0_defFlt

    select case(inversionTest)
    case(invertXYZ)
      det = momVec(xx) * (momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)) &
            - momVec(yy) * momVec(xz) * momVec(xz) &
            - momVec(zz) * momVec(xy) * momVec(xy) &
            + 2 * momVec(xy) * momVec(xz) * momVec(yz)
      invM(xx) = real(momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz),defFlt)
      invM(xy) = real(momVec(xz) * momVec(yz) - momVec(xy) * momVec(zz),defFlt)
      invM(xz) = real(momVec(xy) * momVec(yz) - momVec(yy) * momVec(xz),defFlt)
      invM(yy) = real(momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz),defFlt)
      invM(yz) = real(momVec(xy) * momVec(xz) - momVec(xx) * momVec(yz),defFlt)
      invM(zz) = real(momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy),defFlt)

    case(invertYZ)
      det = momVec(yy) * momVec(zz) - momVec(yz) * momVec(yz)
      invM(yy) = real(momVec(zz),defFlt)
      invM(yz) = real(-momVec(yz),defFlt)
      invM(zz) = real(momVec(yy),defFlt)

    case(invertXY)
      det = momVec(xx) * momVec(yy) - momVec(xy) * momVec(xy)
      invM(xx) = real(momVec(yy),defFlt)
      invM(xy) = real(-momVec(xy),defFlt)
      invM(yy) = real(momVec(xx),defFlt)

    case(invertXZ)
      det = momVec(xx) * momVec(zz) - momVec(xz) * momVec(xz)
      invM(xx) = real(momVec(zz),defFlt)
      invM(xz) = real(-momVec(xz),defFlt)
      invM(zz) = real(momVec(xx),defFlt)

    case(invertX)
      det = momVec(xx)
      invM(xx) = 1.0_defFlt

    case(invertY)
      det = momVec(yy)
      invM(yy) = 1.0_defFlt

    case(invertZ)
      det = momVec(zz)
      invM(zz) = 1.0_defFlt

    case default
      det = ONE
    end select

    one_det = real(ONE/det, defFlt)
    invM = invM * one_det
    
    ! Check for zero determinant
    if (abs(det) < det_tolerance) invM = 0.0_defFlt

    end associate

  end function invertMatrix

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
    if ((k1 <= 0) .or. (k1 > 5)) k1 = k0

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
    if (.not. self % wasFound(cIdx)) return
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
      case(linearIso)
        call self % resetFluxesLinearIso()
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
  subroutine resetFluxesLinearIso(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
      self % prevX(idx) = self % scalarX(idx)
      self % scalarX(idx) = 0.0_defFlt
      self % prevY(idx) = self % scalarY(idx)
      self % scalarY(idx) = 0.0_defFlt
      self % prevZ(idx) = self % scalarZ(idx)
      self % scalarZ(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesLinearIso

  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! for flat sources with anisotropic scattering
  !! TODO: add additional vectors of interest
  !!
  subroutine resetFluxesFlatAni(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do
    do idx = 1, size(self % scalarFlux)
      self % prevFlux(idx) = self % scalarFlux(idx)
      self % scalarFlux(idx) = 0.0_defFlt
    end do
    !$omp end parallel do

  end subroutine resetFluxesFlatAni
  
  !!
  !! Sets prevFlux to scalarFlux and zero's scalarFlux
  !! for linear sources with flat anisotropic scattering
  !! TODO: add additional vectors of interest
  !!
  subroutine resetFluxesLIFA(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt)              :: idx

    !$omp parallel do
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
    character(100), parameter      :: Here = 'accumulateFLuxScores (arraysRR_class.f90)'

    select case(self % simulationType)
      case(flatIso, flatAni)
        call self % accumulateFluxScoresFlat()
      case(linearIso, linearAni)
        call self % accumulateFluxScoresLinear()
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

  end subroutine accumulateFluxScores

  !!
  !! Accumulate flux scores for stats
  !!
  subroutine accumulateFluxScoresFlat(self)
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

  end subroutine accumulateFluxScoresFlat
  
  !!
  !! Accumulate flux scores for stats
  !! Includes linear components
  !!
  subroutine accumulateFluxScoresLinear(self)
    class(arraysRR), intent(inout) :: self
    real(defReal), save            :: flux
    integer(shortInt)              :: idx
    !$omp threadprivate(flux)

    !$omp parallel do schedule(static)
    do idx = 1, size(self % scalarFlux)
      flux = real(self % scalarFlux(idx),defReal)
      self % fluxScores(1, idx) = self % fluxScores(1, idx) + flux
      self % fluxScores(2, idx) = self % fluxScores(2, idx) + flux * flux
      flux = real(self % scalarX(idx),defReal)
      self % xScores(1, idx) = self % xScores(1, idx) + flux
      self % xScores(2, idx) = self % xScores(2, idx) + flux * flux
      flux = real(self % scalarY(idx),defReal)
      self % yScores(1, idx) = self % yScores(1, idx) + flux
      self % yScores(2, idx) = self % yScores(2, idx) + flux * flux
      flux = real(self % scalarZ(idx),defReal)
      self % zScores(1, idx) = self % zScores(1, idx) + flux
      self % zScores(2, idx) = self % zScores(2, idx) + flux * flux
    end do
    !$omp end parallel do

  end subroutine accumulateFluxScoresLinear
  
  !!
  !! Finalise results
  !!
  subroutine finaliseFluxScores(self, it, totalIt)
    class(arraysRR), intent(inout) :: self
    integer(shortInt), intent(in)  :: it
    integer(shortInt), intent(in)  :: totalIt
    character(100), parameter      :: Here = 'finaliseFLuxScores (arraysRR_class.f90)'

    select case(self % simulationType)
      case(flatIso, flatAni)
        call self % finaliseFluxScoresFlat(it)
      case(linearIso, linearAni)
        call self % finaliseFluxScoresLinear(it)
      case default
        call fatalError(Here,'Unsupported simulation type requested')
    end select

  end subroutine finaliseFluxScores
  
  !!
  !! Finalise flux scores for stats
  !!
  subroutine finaliseFluxScoresFlat(self,it)
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

  end subroutine finaliseFluxScoresFlat
  
  !!
  !! Finalise flux scores for stats
  !! Includes linear components
  !!
  subroutine finaliseFluxScoresLinear(self,it)
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
      
      self % xScores(1, idx) = self % xScores(1, idx) * N1
      self % xScores(2, idx) = self % xScores(2, idx) * N1
      self % xScores(2, idx) = Nm1 * (self % xScores(2, idx) - &
            self % xScores(1, idx) * self % xScores(1, idx)) 
      if (self % xScores(2, idx) <= ZERO) then
        self % xScores(2, idx) = ZERO
      else
        self % xScores(2, idx) = sqrt(self % xScores(2, idx))
      end if
      
      self % yScores(1, idx) = self % yScores(1, idx) * N1
      self % yScores(2, idx) = self % yScores(2, idx) * N1
      self % yScores(2, idx) = Nm1 * (self % yScores(2, idx) - &
            self % yScores(1, idx) * self % yScores(1, idx)) 
      if (self % yScores(2, idx) <= ZERO) then
        self % yScores(2, idx) = ZERO
      else
        self % yScores(2, idx) = sqrt(self % yScores(2, idx))
      end if
      
      self % zScores(1, idx) = self % zScores(1, idx) * N1
      self % zScores(2, idx) = self % zScores(2, idx) * N1
      self % zScores(2, idx) = Nm1 * (self % zScores(2, idx) - &
            self % zScores(1, idx) * self % zScores(1, idx)) 
      if (self % zScores(2, idx) <= ZERO) then
        self % zScores(2, idx) = ZERO
      else
        self % zScores(2, idx) = sqrt(self % zScores(2, idx))
      end if
    end do
    !$omp end parallel do

  end subroutine finaliseFluxScoresLinear

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
        
      vol = self % volume(cIdx)
      !if (vol < volume_tolerance) cycle
      vol = vol * self % totalVolume
      if (.not. self % wasFound(cIdx)) cycle

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
    
    if (allocated(self % fixedSource)) then
      do g = 1, self % nG
        name = 'fixedSource_'//numToChar(g)
        !$omp parallel do schedule(static)
        do cIdx = 1, self % nCells
          resVec(cIdx) = real(self % getFixedSource(cIdx,g),defReal)
        end do
        !$omp end parallel do
        call viz % addVTKData(resVec,name)
      end do
    end if

    ! Output final volume estimates
    name = 'volume'
    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      resVec(cIdx) = self % volume(cIdx) * self % totalVolume
    end do
    !$omp end parallel do
    call viz % addVTKData(resVec,name)

    ! Output material IDs
    name = 'material'
    !$omp parallel do schedule(static)
    do cIdx = 1, self % nCells
      resVec(cIdx) = self % geom % geom % graph % getMatFromUID(cIdx) 
    end do
    !$omp end parallel do
    call viz % addVTKData(resVec,name)

    call viz % finaliseVTK()

  end subroutine outputToVTK

  !!
  !! Output integrals over given materials
  !!
  subroutine outputMaterialIntegrals(self, out, matNames)
    class(arraysRR), intent(in)                  :: self
    class(outputFile), intent(inout)             :: out
    character(nameLen), dimension(:), intent(in) :: matNames
    integer(shortInt)                            :: i, matIdx
    character(nameLen)                           :: name
    real(defReal)                                :: integral, intSD

    name = 'integral'
    call out % startBlock(name)

    do i = 1, size(matNames)
    
      call out % startArray(matNames(i), [1])

      ! Ensure name corresponds to a material present in the geometry
      matIdx = self % XSData % getIdxFromName(matNames(i))

      if (matIdx /= -1) then
        call self % materialIntegral(matIdx, integral, intSD)
      else
        integral = ZERO
        intSD    = ZERO
        print *,'WARNING: Material '//matNames(i)//' not found during integration.'
      end if

      call out % addResult(integral, intSD)
      call out % endArray()

    end do
    call out % endBlock()

  end subroutine outputMaterialIntegrals

  !!
  !! Perform an integral over a given material for all energies
  !!
  subroutine materialIntegral(self, matIdx, integral, intSD)
    class(arraysRR), intent(in)   :: self
    integer(shortInt), intent(in) :: matIdx
    real(defReal), intent(out)    :: integral, intSD
    integer(shortInt)             :: cIdx
    real(defReal), save           :: vol
    integer(shortInt), save       :: mIdx, g
    !$omp threadprivate(vol, mIdx, g)

    integral = ZERO
    intSD    = ZERO

    !$omp parallel do reduction(+: integral, intSD)
    do cIdx = 1, self % nCells

      mIdx = self % geom % geom % graph % getMatFromUID(cIdx) 
      if (mIdx /= matIdx) cycle

      vol = self % volume(cIdx)
      !if (vol < volume_tolerance) cycle
      vol = vol * self % totalVolume
      if (.not. self % wasFound(cIdx)) cycle

      do g = 1, self % nG
        integral = integral + self % getFluxScore(cIdx, g) * vol

        ! Neglects uncertainty in the volume estimate
        intSD = intSD + self % getFluxSD(cIdx, g)**2 * vol * vol
      end do

    end do
    !$omp end parallel do

    if (intSD >= ZERO) intSD = sqrt(intSD)
    if (abs(integral) > ZERO) intSD = intSD / abs(integral)


  end subroutine materialIntegral
  
  !!
  !! Output fluxes at given points
  !!
  subroutine outputPointFluxes(self, out, points, names)
    class(arraysRR), intent(in)                  :: self
    class(outputFile), intent(inout)             :: out
    real(defReal), dimension(:,:), intent(in)    :: points
    character(nameLen), dimension(:), intent(in) :: names
    integer(shortInt)                            :: i
    character(nameLen)                           :: name
    real(defReal), dimension(self % nG)          :: flux, fluxSD
    character(100), parameter                    :: Here = 'outputPointFluxes (arraysRR_class.f90)'

    name = 'pointFlux'
    call out % startBlock(name)

    ! Ensure points and names have the correction dimensions
    if (size(points,1) /= 3) call fatalError(Here, 'Points are not 3D.')
    if (size(points,2) /= size(names)) call fatalError(Here, &
            'Different numbers of sample points to sample names.')

    do i = 1, size(names)
    
      call out % startArray(names(i), [self % nG])
      call self % getFluxAtAPoint(points(:, i), flux, fluxSD)
      call out % addResult(flux, fluxSD)
      call out % endArray()

    end do
    call out % endBlock()

  end subroutine outputPointFluxes

  !!
  !! Returns the flux vector at a point in space
  !!
  subroutine getFluxAtAPoint(self, r, flux, fluxSD)
    class(arraysRR), intent(in)                      :: self
    real(defReal), dimension(3), intent(in)          :: r
    real(defReal), dimension(self % nG), intent(out) :: flux
    real(defReal), dimension(self % nG), intent(out) :: fluxSD
    integer(shortInt)                                :: g, matIdx, cIdx, i
    real(defReal), dimension(3)                      :: mom, momSD, centroid, fluxGrad
    real(defFlt), dimension(matSize)                 :: invM

    ! Identify cell at the given point
    call self % geom % whatIsAt(matIdx, cIdx, r)

    if (cIdx > 0) then
      do g = 1, self % nG
      
        flux(g) = self % getFluxScore(cIdx, g)
        fluxSD(g) = self % getFluxSD(cIdx, g)

        ! Include linear moments if available
        if ((self % simulationType == linearIso) .or. &
                (self % simulationType == linearAni)) then

          fluxSD(g) = fluxSD(g) * fluxSD(g)

          mom = self % getFluxMoments(cIdx, g)
          momSD = self % getFluxMomentSDs(cIdx, g)
          centroid = self % getCentroid(cIdx)
          invM = self % invertMatrix(cIdx)
      
          ! Get flux gradients
          fluxGrad(x) = real(invM(xx) * mom(x) + invM(xy) * mom(y) + invM(xz) * mom(z), defReal)
          fluxGrad(y) = real(invM(xy) * mom(x) + invM(yy) * mom(y) + invM(yz) * mom(z), defReal)
          fluxGrad(z) = real(invM(xz) * mom(x) + invM(yz) * mom(y) + invM(zz) * mom(z), defReal)

          ! Note this will not correctly estimate uncertainty as moments are covariant
          do i = 1, 3
            flux(g) = flux(g) + fluxGrad(i) * (r(i) - centroid(i))
            ! Not sure exactly how to propagate uncertainties - need to do some maths
            fluxSD(g) = fluxSD(g) + momSD(i)**2 * (r(i) - centroid(i))**2
          end do

          if (fluxSD(g) > ZERO) fluxSD(g) = sqrt(fluxSD(g))

        end if

      end do

    else
      print *,'WARNING: No cell found at position '//numToChar(r)
      flux = -ONE
      fluxSD = -ONE
    end if

  end subroutine getFluxAtAPoint
  
  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(arraysRR), intent(inout) :: self
    integer(shortInt) :: i

    ! Clean standard contents
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
    
    ! Clean LS contents
    if(allocated(self % scalarX)) deallocate(self % scalarX)
    if(allocated(self % scalarX)) deallocate(self % scalarY)
    if(allocated(self % scalarX)) deallocate(self % scalarZ)
    if(allocated(self % prevX)) deallocate(self % prevX)
    if(allocated(self % prevY)) deallocate(self % prevY)
    if(allocated(self % prevZ)) deallocate(self % prevZ)
    if(allocated(self % sourceX)) deallocate(self % sourceX)
    if(allocated(self % sourceY)) deallocate(self % sourceY)
    if(allocated(self % sourceZ)) deallocate(self % sourceZ)
    if(allocated(self % momMat)) deallocate(self % momMat)
    if(allocated(self % momTracks)) deallocate(self % momTracks)
    if(allocated(self % centroid)) deallocate(self % centroid)
    if(allocated(self % centroidTracks)) deallocate(self % centroidTracks)
    if(allocated(self % xScores)) deallocate(self % xScores)
    if(allocated(self % yScores)) deallocate(self % yScores)
    if(allocated(self % zScores)) deallocate(self % zScores)

    if(allocated(self % locks)) then
      do i = 1, self % nCells
#ifdef _OPENMP
        call OMP_destroy_lock(self % locks(i))
#endif
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
    self % totalVolume = ZERO

  end subroutine kill

end module arraysRR_class
