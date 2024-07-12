module dataRR_class

  use numPrecision
  use universalVariables
  use rng_class,                      only : RNG
  use genericProcedures,              only : fatalError

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat => nMat, mm_matName => matName
  use materialHandle_inter,           only : materialHandle
  use baseMgNeutronDatabase_class,    only : baseMgNeutronDatabase
  use baseMgNeutronMaterial_class,    only : baseMgNeutronMaterial, baseMgNeutronMaterial_CptrCast
  
  implicit none
  private

  !!
  !! Nuclear data in a random ray-friendly format.
  !!
  !! Stores data and provides access in a manner which is more performant
  !! than is done for MC MG data at present.
  !!
  !! TODO: Add kinetic data and higher-order scattering matrices
  !!
  type, public :: dataRR
    private
    ! Components
    integer(shortInt)                     :: nG     = 0
    integer(shortInt)                     :: nG2    = 0
    integer(shortInt)                     :: nMat   = 0

    ! Data space - absorb all nuclear data for speed
    real(defFlt), dimension(:), allocatable       :: sigmaT
    real(defFlt), dimension(:), allocatable       :: nuSigmaF
    real(defFlt), dimension(:), allocatable       :: sigmaF
    real(defFlt), dimension(:), allocatable       :: sigmaS
    real(defFlt), dimension(:), allocatable       :: chi
    logical(defBool), dimension(:), allocatable   :: fissile
    character(nameLen), dimension(:), allocatable :: names

    ! Optional kinetic data
    logical(defBool)                        :: doKinetics  = .false.
    integer(shortInt)                       :: nP = 0
    real(defFlt), dimension(:), allocatable :: chiD
    real(defFlt), dimension(:), allocatable :: chiP
    real(defFlt), dimension(:), allocatable :: beta
    real(defFlt), dimension(:), allocatable :: invSpeed

    ! Optional higher-order scattering matrices up to P3
    real(defFlt), dimension(:), allocatable :: sigmaS1
    real(defFlt), dimension(:), allocatable :: sigmaS2
    real(defFlt), dimension(:), allocatable :: sigmaS3

  contains
    
    procedure :: init
    procedure :: kill

    ! TODO: add an XS update procedure, e.g., given multiphysics
    ! TODO: add full handling of kinetic data
    ! TODO: add full handling of higher-order anisotropy

    ! Access procedures
    procedure :: getProdPointers
    !procedure :: getAllPointers
    procedure :: getTotalPointer
    procedure :: getNuFissPointer
    procedure :: getChiPointer
    procedure :: getScatterPointer
    procedure :: getScatterVecPointer
    procedure :: getTotalXS
    procedure :: getFissionXS
    procedure :: getScatterXS
    procedure :: getNG
    procedure :: getNMat
    procedure :: getNPrec
    procedure :: getName
    procedure :: getIdxFromName
    procedure :: isFissile

    ! Private procedures
    procedure, private :: getIdxs
    procedure, private :: getScatterIdxs
    !procedure, private :: getKineticIdxs


  end type dataRR

contains

  !!
  !! Initialise necessary nuclear data.
  !! Can optionally include kinetic parameters.
  !!
  subroutine init(self, db, doKinetics, loud)
    class(dataRR), intent(inout)                     :: self
    class(baseMgNeutronDatabase),pointer, intent(in) :: db
    logical(defBool), intent(in)                     :: doKinetics
    logical(defBool), intent(in)                     :: loud
    integer(shortInt)                                :: g, g1, m, matP1, aniOrder
    type(RNG)                                        :: rand
    logical(defBool)                                 :: fiss
    class(baseMgNeutronMaterial), pointer            :: mat
    class(materialHandle), pointer                   :: matPtr
    character(100), parameter :: Here = 'init (dataRR_class.f90)'

    self % doKinetics = doKinetics

    ! Store number of energy groups for convenience
    self % nG = db % nGroups()
    self % nG2 = self % nG * self % nG

    ! Initialise local nuclear data
    ! Allocate nMat + 1 materials to catch any undefined materials
    ! TODO: clean nuclear database afterwards! It is no longer used
    !       and takes up memory.
    self % nMat = mm_nMat()
    matP1 = self % nMat + 1
    allocate(self % sigmaT(matP1 * self % nG))
    self % sigmaT = 0.0_defFlt
    allocate(self % nuSigmaF(matP1 * self % nG))
    self % nuSigmaF = 0.0_defFlt
    allocate(self % sigmaF(matP1 * self % nG))
    self % sigmaF = 0.0_defFlt
    allocate(self % chi(matP1 * self % nG))
    self % chi = 0.0_defFlt
    allocate(self % sigmaS(matP1 * self % nG * self % nG))
    self % sigmaS = 0.0_defFlt
    allocate(self % fissile(matP1))
    self % fissile = .false.
    allocate(self % names(matP1))
    self % names = 'unnamed'

    ! Create a dummy RNG to satisfy the mgDatabase access interface
    call rand % init(1_longInt)

    if (loud) print *,'Initialising random ray nuclear data'
    do m = 1, self % nMat
      matPtr  => db % getMaterial(m)
      mat     => baseMgNeutronMaterial_CptrCast(matPtr)
      fiss = .false.
      do g = 1, self % nG
        self % sigmaT(self % nG * (m - 1) + g) = real(mat % getTotalXS(g, rand),defFlt)
        self % nuSigmaF(self % nG * (m - 1) + g) = real(mat % getNuFissionXS(g, rand),defFlt)
        self % sigmaF(self % nG * (m - 1) + g) = real(mat % getFissionXS(g, rand),defFlt)
        if (self % nuSigmaF(self % nG * (m - 1) + g) > 0) fiss = .true.
        self % chi(self % nG * (m - 1) + g) = real(mat % getChi(g, rand),defFlt)
        ! Include scattering multiplicity
        do g1 = 1, self % nG
          self % sigmaS(self % nG * self % nG * (m - 1) + self % nG * (g - 1) + g1)  = &
                  real(mat % getScatterXS(g1, g, rand) * mat % scatter % prod(g1, g) , defFlt)
        end do
      end do
      self % fissile(m) = fiss
      self % names(m) = mm_matName(m)
    end do

    ! Initialise data necessary for kinetic/noise calculations
    if (self % doKinetics) then
      print *,'Including kinetic data'
      call fatalError(Here,'Kinetic data not yet supported')

    end if

    ! Initialise higher-order scattering matrices
    ! TODO: read anisotropy order from database
    aniOrder = 0
    if (aniOrder > 0) then
      print *,'Including anisotropic scattering data'
      call fatalError(Here,'Anisotropy not yet supported')

    end if

  end subroutine init

  !!
  !! Calculate the lower and upper indices for accessing the XS array
  !! (excluding scattering and kinetic data)
  !!
  pure subroutine getIdxs(self, matIdx, idx1, idx2)
    class(dataRR), intent(in)      :: self
    integer(shortInt), intent(in)  :: matIdx
    integer(shortInt), intent(out) :: idx1, idx2

    idx1 = (matIdx - 1) * self % nG + 1
    idx2 = matIdx * self % nG 

  end subroutine getIdxs

  !!
  !! Calculate the lower and upper indices for accessing the scattering XS array
  !!
  pure subroutine getScatterIdxs(self, matIdx, idx1, idx2)
    class(dataRR), intent(in)      :: self
    integer(shortInt), intent(in)  :: matIdx
    integer(shortInt), intent(out) :: idx1, idx2

    idx1 = (matIdx - 1) * self % nG2 + 1
    idx2 = matIdx * self % nG2 

  end subroutine getScatterIdxs

  !!
  !! Return if a material is fissile
  !!
  elemental function isFissile(self, matIdx) result(isFiss)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx
    logical(defBool)              :: isFiss
    integer(shortInt)             :: mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    isFiss = self % fissile(mIdx)

  end function isFissile
  
  !!
  !! Return the number of groups
  !!
  elemental function getNG(self) result(nG)
    class(dataRR), intent(in) :: self
    integer(shortInt)         :: nG

    nG = self % nG

  end function getNG

  !!
  !! Return the number of materials
  !!
  elemental function getNMat(self) result(nM)
    class(dataRR), intent(in) :: self
    integer(shortInt)         :: nM

    nM = self % nMat

  end function getNMat

  !!
  !! Return the number of precursors
  !!
  elemental function getNPrec(self) result(nP)
    class(dataRR), intent(in) :: self
    integer(shortInt)         :: nP

    nP = self % nP

  end function getNPrec
  
  !!
  !! Return the name of a given material
  !!
  elemental function getName(self, matIdx) result(matName)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx
    integer(shortInt)             :: mIdx
    character(nameLen)            :: matName

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    matName = self % names(mIdx)

  end function getName
  
  !!
  !! Return the index of a material given its name
  !!
  elemental function getIdxFromName(self, matName) result(matIdx)
    class(dataRR), intent(in)      :: self
    character(nameLen), intent(in) :: matName
    integer(shortInt)              :: matIdx

    do matIdx = 1, self % nMat 
      if (self % names(matIdx) == matName) return
    end do
    matIdx = -1

  end function getIdxFromName
  
  !!
  !! Get scatter pointer
  !!
  subroutine getScatterPointer(self, matIdx, sigS)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: sigS
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getScatterIdxs(mIdx, idx1, idx2)
    sigS => self % sigmaS(idx1:idx2)

  end subroutine getScatterPointer
  
  !!
  !! Get scatter vector pointer
  !!
  subroutine getScatterVecPointer(self, matIdx, gOut, sigS)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    integer(shortInt), intent(in)                    :: gOut
    real(defFlt), dimension(:), pointer, intent(out) :: sigS
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    idx1 = (matIdx - 1) * self % nG2 + (gOut - 1) * self % nG + 1
    idx2 = (matIdx - 1) * self % nG2 + gOut * self % nG 
    sigS => self % sigmaS(idx1:idx2)

  end subroutine getScatterVecPointer


  !!
  !! Get chi pointer
  !!
  subroutine getChiPointer(self, matIdx, chi)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: chi
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getIdxs(mIdx, idx1, idx2)
    chi => self % chi(idx1:idx2)

  end subroutine getChiPointer

  !!
  !! Return pointers to all commonly used XSs for neutron production
  !! This is done for a given material, across all energies
  !!
  subroutine getProdPointers(self, matIdx, nuSigF, sigS, chi)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: nuSigF, sigS, chi
    integer(shortInt)                                :: idx1, idx2, idx1s, idx2s, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getIdxs(mIdx, idx1, idx2)
    call self % getScatterIdxs(mIdx, idx1s, idx2s)
    nuSigF => self % nuSigmaF(idx1:idx2)
    chi    => self % chi(idx1:idx2)
    sigS   => self % sigmaS(idx1s:idx2s)

  end subroutine getProdPointers
  
  !!
  !! Return pointers to only the total XS
  !! This is done for a given material, across all energies
  !!
  subroutine getTotalPointer(self, matIdx, sigT)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: sigT
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getIdxs(mIdx, idx1, idx2)
    sigT => self % sigmaT(idx1:idx2)

  end subroutine getTotalPointer
  
  !!
  !! Return pointers to only the nuFission XS
  !! This is done for a given material, across all energies
  !!
  subroutine getNuFissPointer(self, matIdx, nuFiss)
    class(dataRR), target, intent(in)                :: self
    integer(shortInt), intent(in)                    :: matIdx
    real(defFlt), dimension(:), pointer, intent(out) :: nuFiss
    integer(shortInt)                                :: idx1, idx2, mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    call self % getIdxs(mIdx, idx1, idx2)
    nuFiss => self % nuSigmaF(idx1:idx2)

  end subroutine getNuFissPointer

  !!
  !! Return total XS in a given material and group
  !!
  elemental function getTotalXS(self, matIdx, g) result(sigT)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx, g
    real(defFlt)                  :: sigT
    integer(shortInt)             :: mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    sigT = self % sigmaT((mIdx - 1) * self % nG + g)

  end function getTotalXS
  
  !!
  !! Return fission XS in a given material and group
  !!
  elemental function getFissionXS(self, matIdx, g) result(sigF)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx, g
    real(defFlt)                  :: sigF
    integer(shortInt)             :: mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    sigF = self % sigmaF((mIdx - 1) * self % nG + g)

  end function getFissionXS
  
  !!
  !! Return scatter XS in a given material, ingoing group, and outgoing group
  !!
  elemental function getScatterXS(self, matIdx, gIn, gOut) result(sigS)
    class(dataRR), intent(in)     :: self
    integer(shortInt), intent(in) :: matIdx, gIn, gOut
    real(defFlt)                  :: sigS
    integer(shortInt)             :: mIdx

    if (matIdx > self % nMat) then
      mIdx = self % nMat + 1
    else
      mIdx = matIdx
    end if
    sigS = self % sigmaS((mIdx - 1) * self % nG2 + self % nG * (gIn - 1) + gOut)

  end function getScatterXS


  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(dataRR), intent(inout) :: self

    ! Clean contents
    self % nG         = 0
    self % nG2        = 0
    self % nMat       = 0
    self % nP         = 0
    self % doKinetics = .false.
    if(allocated(self % sigmaT)) deallocate(self % sigmaT)
    if(allocated(self % sigmaS)) deallocate(self % sigmaS)
    if(allocated(self % nuSigmaF)) deallocate(self % nuSigmaF)
    if(allocated(self % sigmaF)) deallocate(self % sigmaF)
    if(allocated(self % chi)) deallocate(self % chi)
    if(allocated(self % fissile)) deallocate(self % fissile)
    if(allocated(self % names)) deallocate(self % names)
    if(allocated(self % chiD)) deallocate(self % chiD)
    if(allocated(self % chiP)) deallocate(self % chiP)
    if(allocated(self % beta)) deallocate(self % beta)
    if(allocated(self % invSpeed)) deallocate(self % invSpeed)
    if(allocated(self % sigmaS1)) deallocate(self % sigmaS1)
    if(allocated(self % sigmaS2)) deallocate(self % sigmaS2)
    if(allocated(self % sigmaS3)) deallocate(self % sigmaS3)

  end subroutine kill

end module dataRR_class
