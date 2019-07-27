module isotropicMG_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures,              only : fatalError, linFind, searchError
  use dictionary_class,               only : dictionary
  use IOdictionary_class,             only : IOdictionary
  use RNG_class,                      only : RNG
  use particle_class,                 only : particle

  use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG
  use outscatterCDF_class,            only : outscatterCDF
  use xsMacroSet_class,               only : xsMacroSet, xsMacroSet_ptr
  use releaseMatrixMG_class,          only : releaseMatrixMG

  implicit none
  private

  !!
  !! Small type to store decription of the material
  !!
  type, private :: materialData
    character(nameLen) :: matName
    logical(defBool)   :: isFissile = .false.
    logical(defBool)   :: isActive  = .true.
  end type materialData

  !!
  !! DETAILED EXPLENETAION WILL GO HERE
  !!
  type, public, extends(perMaterialNuclearDataMG) :: isotropicMG
    private
    integer(shortInt)   :: nG
    integer(shortInt)   :: nMat
    type(materialData), dimension(:), allocatable      :: matData         ! (matIdx)

    ! Note that the following are protected members, thus prefix P_ and public attribute
    type(xsMacroSet),dimension(:,:),pointer,public          :: P_XSs => null()   ! (energyGroup,matIdx)
    type(outscatterCDF), dimension(:,:), allocatable,public :: P_transferMatrix  ! (energyGroup,matIdx)
    type(outscatterCDF), dimension(:), allocatable,public   :: P_chiValues       ! (matIdx)
    type(releaseMatrixMG), dimension(:), allocatable,public :: P_releaseData     ! (matIdx)
    real(defReal), dimension(:), allocatable, public        :: P_majorantXS    ! (energyGroup)

  contains
    procedure :: init
    procedure :: kill

    ! Nuclear Data Interface Procedures
    procedure :: getIdx
    procedure :: getName

    ! transportNuclearData Interface Procedures
    procedure :: getTransXS_G
    procedure :: getMajorantXS_G
    procedure :: getTotalMatXS_G
    procedure :: getMatMacroXS_G
    procedure :: isFissileMat
    procedure :: initFissionSite
    procedure :: setActiveMaterials

    ! perMaterialNuclearDataMG InterfaceProcedures
    procedure :: releaseAt
    procedure :: sampleMuGout

    ! Inquiry procedures
    procedure, non_overridable :: numG
    procedure, non_overridable :: numMat

    !* TYPE PROCEDURES *!
    procedure, private :: readMaterial
    procedure, private :: activeIdx
    procedure, private :: calculateMajorant

  end type isotropicMG

  !! Procedures that are used by subclasses of isotropicMG
  public :: init
  public :: kill
  public :: readMaterial
  public :: activeIdx

contains

  !!
  !! Load material data from a dictionary
  !!
  subroutine init(self,dict, matNames)
    class(isotropicMG), intent(inout)           :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen),dimension(:),intent(in)  :: matNames
    class(dictionary), pointer                  :: tempDict => null()
    integer(shortInt)                           :: nG,nMat,matIdx

    ! Read number of energy groups and materials
    call dict % get(nG, 'numberOfGroups')
    nMat = size( matNames )

    self % nG   = nG
    self % nMat = nMat

    ! Allocate space fo XS data
    if (associated( self % P_XSs           )) deallocate (self % P_XSs            )
    if (allocated(  self % P_transferMatrix)) deallocate (self % P_transferMatrix )
    if (allocated(  self % P_chiValues     )) deallocate (self % P_chiValues      )
    if (allocated(  self % P_releaseData   )) deallocate (self % P_releaseData    )
    if (allocated(  self % matData       )) deallocate (self % matData        )

    allocate( self % P_XSs(nG,nMat)           )
    allocate( self % P_transferMatrix(ng,nMat))
    allocate( self % P_chiValues(nMat)        )
    allocate( self % P_releaseData(nMat)      )
    allocate( self % matData(nMat)          )

    ! Read material names
    self % matData(:) % matName = matNames

    ! Read individual material data
    do matIdx=1,nMat
      tempDict => dict % getDictPtr(self % matData(matIdx) % matName)
      call self % readMaterial(tempDict, matIdx, self % nG)
    end do

    call self % calculateMajorant()

  end subroutine init

  !!
  !! Deallocate occupied space
  !!
  elemental subroutine kill(self)
    class(isotropicMG), intent(inout) :: self

    if (associated( self % P_XSs           )) deallocate (self % P_XSs            )
    if (allocated(  self % P_transferMatrix)) deallocate (self % P_transferMatrix )
    if (allocated(  self % P_chiValues     )) deallocate (self % P_chiValues      )
    if (allocated(  self % P_releaseData   )) deallocate (self % P_releaseData    )
    if (allocated(  self % matData       )) deallocate (self % matData        )

  end subroutine kill

  !!
  !! Return matIdx of material with matName
  !!
  function getIdx(self,matName) result(matIdx)
    class(isotropicMG), intent(in)      :: self
    character(*), intent(in)            :: matName
    integer(shortInt)                   :: matIdx
    character(100), parameter           :: Here ='getMatIdx (isotropicMG_class.f90)'

    matIdx = linFind(self % matData % matName, matName)
    call searchError(matIdx,Here)

  end function getIdx

  !!
  !! Returns matName of material with matIdx
  !!
  function getName(self,matIdx) result(matName)
    class(isotropicMG), intent(in)      :: self
    integer(shortInt), intent(in)       :: matIdx
    character(nameLen)                  :: matName

    matName = self % matData(matIdx) % matName

  end function getName

  !!
  !! Return transport XS (in general diffrent from total XS)
  !!
  function getTransXS_G(self,G,matIdx) result(xs)
    class(isotropicMG), intent(inout)  :: self
    integer(shortInt), intent(in)      :: G
    integer(shortInt), intent(in)      :: matIdx
    real(defReal)                      :: xs

    xs = self % P_XSs(G,matIdx) % totalXS

  end function getTransXS_G

  !!
  !! Return majorant XS (in general should be largest of TRANSPORT XSs)
  !!
  function getMajorantXS_G(self,G) result (xs)
    class(isotropicMG), intent(inout)  :: self
    integer(shortInt), intent(in)      :: G
    real(defReal)                      :: xs

    xs = self % P_majorantXS(G)

  end function getMajorantXS_G

  !!
  !! Return total XS of material
  !!
  function getTotalMatXS_G(self,G,matIdx) result (xs)
    class(isotropicMG), intent(inout)  :: self
    integer(shortInt), intent(in)      :: G
    integer(shortInt), intent(in)      :: matIdx
    real(defReal)                      :: xs

    xs = self % P_XSs(G,matIdx) % totalXS

  end function getTotalMatXS_G

  !!
  !! Attach pointer to approperiate XS data
  !!
  subroutine getMatMacroXS_G(self,macroXS,G,matIdx)
    class(isotropicMG), intent(inout)            :: self
    type(xsMacroSet_ptr),intent(inout)           :: macroXS
    integer(shortInt),intent(in)                 :: G
    integer(shortInt),intent(in)                 :: matIdx

    macroXS = self % P_XSs(G,matIdx)

  end subroutine getMatMacroXS_G

  !!
  !! Returns .true. if material is fissile
  !!
  function isFissileMat(self,matIdx) result(isIt)
    class(isotropicMG), intent(in)  :: self
    integer(shortInt), intent(in)   :: matIdx
    logical(defBool)                :: isIt

    select case(matIdx)
      case(OUTSIDE_MAT)
        isIT = .false.

      case(VOID_MAT)
        isIt = .false.

      case default
        isIt = self % matData(matIdx) % isFissile
    end select
  end function isFissileMat

  !!
  !! Procedure to generate a fission site from a fissile material.
  !! Necassary in initialisation of an eigenvalue calculation:
  !!
  subroutine initFissionSite(self,p,r)
    class(isotropicMG), intent(in)         :: self
    class(particle), intent(inout)         :: p
    real(defReal),dimension(3), intent(in) :: r
    integer(shortInt)                      :: G_out, matIdx
    real(defReal)                          :: mu, phi, r1
    integer(shortInt),parameter            :: G_in = 1
    character(100), parameter      :: Here = 'initFissionSite (isotropicMG_class.f90)'

    matIdx = p % matIdx()
    p % isMG = .true.

    ! Determine if material is fissile
    if ( .not.self % isFissileMat(matIdx) ) then
      call fatalError(Here,' Material: '//self % getName(matIdx) //' is not fissile')
    end if

    ! Generate random numbers
    r1 = p % pRNG % get()

    ! Sample outgoing data
    call self % sampleMuGout(mu, G_out, G_in, p % pRNG, macroFission ,matIdx)
    phi = 2*PI * r1

    ! Update particle state
    call p % point([ONE, ZERO, ZERO])
    call p % rotate(mu,phi)
    p % G = G_out
    p % w = ONE
    call p % teleport(r)
    ! *** HAVING FUNCTION TO FOR STTING MATIDX WILL BE GOOD
    p % coords % matIdx = matIdx
    p % isDead = .false.


  end subroutine initFissionSite

  !!
  !! Set all materials that are present in geometry
  !! Active materials will be included in evaluation of majorant XS
  !!
  !! Provide an array of all active material indexes
  !! Order is not significant (assume there is none)
  !!
  subroutine setActiveMaterials(self,matIdxList)
    class(isotropicMG), intent(inout)          :: self
    integer(shortInt),dimension(:), intent(in) :: matIdxList
    integer(shortInt)                          :: i

    ! Switch all materials to inactive
    self % matData % isActive = .false.

    ! Loop through matIdxList and set flags back to active
    do i=1,size(matIdxList)
      self % matData( matIdxList(i) ) % isActive = .true.
    end do

    ! Recalculate Majorant
    call self % calculateMajorant()

  end subroutine setActiveMaterials

  !!
  !! Returns neutron release at G_in for material matIdx and reaction MT
  !!
  function releaseAt(self,G_in,G_out,MT,matIdx) result(nu)
    class(isotropicMG), intent(in)      :: self
    integer(shortInt), intent(in)       :: G_in
    integer(shortInt), intent(in)       :: G_out
    integer(shortInt), intent(in)       :: MT
    integer(shortInt), intent(in)       :: matIdx
    real(defReal)                       :: nu
    character(100), parameter           :: Here = 'releaseAt (isotropicMG_class.f90)'

    select case(MT)
      case(macroAllScatter)
        nu = self % P_releaseData(matIdx) % scatterRelease(G_in,G_out)

      case(macroFission)
        nu = self % P_releaseData(matIdx) % fissionRelease(G_in)

      case default
        call fatalError(Here,'Unrecoginsed MT number')
        nu = ZERO

    end select

  end function releaseAt

  !!
  !! Samples deflection angle in LAB frame and post-colission energy group
  !! This implementation is CRAP. SHOULD BE BRANCHLESS. IMPROVE IT! **##~~##**
  !!
  subroutine sampleMuGout(self,mu,G_out,G_in,rand,MT,matIdx)
    class(isotropicMG), intent(in)     :: self
    real(defReal), intent(out)         :: mu
    integer(shortInt), intent(out)     :: G_out
    integer(shortInt), intent(in)      :: G_in
    class(RNG), intent(inout)          :: rand
    integer(shortInt), intent(in)      :: MT
    integer(shortInt), intent(in)      :: matIdx
    real(defReal)                      :: r1
    character(100), parameter          :: Here ='sampleMuGout (isotropicMG_class.f90)'

    mu = TWO * rand % get() - ONE

    r1 = rand % get()

    select case(MT)
      case(macroAllScatter)
        G_out = self % P_transferMatrix(G_in, matIdx) % invert(r1)

      case(macroFission)
        G_out = self % P_chiValues(matIdx) % invert(r1)

      case default
        call fatalError(Here,'Unrecoginsed MT number')

    end select

  end subroutine sampleMuGout

  !!
  !! Returns number of groups for a data in the isotropicMG
  !!
  pure function numG(self) result(nG)
    class(isotropicMG), intent(in) :: self
    integer(shortInt)              :: nG

    nG = self % nG

  end function numG

  !!
  !! Returns number of materials in the data in the isotropicMG
  !!
  pure function numMat(self) result(nMat)
    class(isotropicMg), intent(in) :: self
    integer(shortInt)              :: nMat

    nMat = self % nMat

  end function numMat

  !!
  !! Read material under index "idx" from dictionary "dict"
  !! For group transfer matrixes indexing convention is (G_out,G_in)
  !! When Matrixes are rank1 then the sequence follows SERPENT(2.1.30) output convention:
  !! G_in->1, G_in->2 ...
  !!
  subroutine readMaterial(self,dict,idx, nG)
    class(isotropicMG), intent(inout) :: self
    class(dictionary), intent(in)     :: dict
    integer(shortInt), intent(in)     :: idx
    integer(shortInt), intent(in)     :: nG
    real(defReal),dimension(:),allocatable       :: tempXS
    real(defReal),dimension(:),allocatable       :: tempXSmatrix_rank1
    real(defReal),dimension(nG,nG)    :: tempXSmatrix
    type(IOdictionary)                :: xsDict
    integer(shortInt)                 :: i
    logical(defBool)                  :: isFissile
    character(pathLen)                :: xsPath
    character(100), parameter         :: Here='readMaterial (isotropicMG_class.f90)'

    ! Obtain path from dict and read into xsDict
    call dict % get(xsPath, 'xsFile')
    call xsDict % initFrom(xsPath)

    ! Check if material is fissile
    isFissile = xsDict % isPresent('fission')
    self % matData(idx) % isFissile = isFissile

    ! Verify size of stored data
    if (xsDict % getSize('capture') /= nG) then
      call fatalError(Here,'capture xs are inconsistant with number of energy groups')

    elseif (xsDict % getSize('scattering_P0') /= nG*nG) then
      call fatalError(Here,'scatter xs are inconsistant with number of energy groups')

    else if (xsDict % getSize('scatteringMultiplicity') /= nG*nG) then
      call fatalError(Here,'scattering production data is inconsistant with number of energy groups')

    end if

    if (isFissile) then
      if (xsDict % getSize('fission') /= nG) then
        call fatalError(Here,'fission xs are inconsistant with number of energy groups')

      else if (xsDict % getSize('chi') /= nG) then
        call fatalError(Here,'chi data is inconsistant with number of energy groups')

      else if (xsDict % getSize('nu') /= nG) then
        call fatalError(Here,'nu data is inconsistant with number of energy groups')

      end if
    end if

    ! Load and store capture XSs
    call xsDict % get(tempXS, 'capture')
    if (any( tempXS < 0.0)) call fatalError(Here,'capture xss are -ve')
    self % P_XSs(:,idx) % captureXS = tempXS

    ! Load and store fission XSs
    if (isFissile) then
      call xsDict % get(tempXS, 'fission')
    else
      tempXS = 0.0
    end if
    if (any( tempXS < 0.0)) call fatalError(Here,'fission xss are -ve')
    self % P_XSs(:,idx) % fissionXS = tempXS

    ! Load and store chi values
    if (isFissile) then
      call xsDict % get(tempXS, 'chi')
    else
      tempXS = 0.0
      tempXS(1) = 1.0 ! Avoid Floating point exception
    end if
    if (any( tempXS < 0.0)) call fatalError(Here,'chi is -ve')
    call self % P_chiValues(idx) % init(tempXS)

    ! Load scattering matrix. Indexing convenction is (G_out,G_in)
    call xsDict % get(tempXSmatrix_rank1, 'scattering_P0')
    if (any( tempXSmatrix_rank1 < 0.0)) call fatalError(Here,'Scattering_P0 is -ve')
    tempXSmatrix = reshape(tempXSmatrix_rank1,[nG, nG])
    do i=1,nG
      call self % P_transferMatrix(i,idx) % init ( tempXSmatrix(:,i) )
    end do
    self % P_XSs(:,idx) % scatterXS = sum(tempXSmatrix,1)

    ! Load production matrix and nu. Indexing convenction is (G_out,G_in).
    call xsDict % get(tempXSmatrix_rank1, 'scatteringMultiplicity')
    if (any( tempXSmatrix_rank1 < 0.0)) call fatalError(Here,'scateringMultiplicity in -ve')
    tempXSmatrix = reshape(tempXSmatrix_rank1,[nG, nG])

    if (isFissile) then
      call xsDict % get(tempXS, 'nu')
    else
      tempXS = 0.0
    end if

    call self % P_releaseData(idx) % init(tempXS, tempXSmatrix)

    ! Calculate nu*Fission
    if (isFissile) then
      tempXS = self % P_XSs(:,idx) % fissionXS * tempXS
    else
      tempXS = 0.0
    end if

    self % P_XSs(:,idx) % nuFissionXS = tempXS

    ! Calculate total XS
    self % P_XSs(:,idx) % totalXS = self % P_XSs(:,idx) % scatterXS + &
                                  self % P_XSs(:,idx) % captureXS + &
                                  self % P_XSs(:,idx) % fissionXS

  end subroutine readMaterial

  !!
  !! Allocates an array with active indices
  !!
  subroutine activeIdx(self,array)
    class(isotropicMG), intent(in)                             :: self
    integer(shortInt), dimension(:),allocatable, intent(inout) :: array
    logical(defBool), dimension(:), allocatable                :: mask
    integer(shortInt)                                          :: i

    if(allocated(array)) deallocate(array)
    mask  = self % matData % isActive
    array = pack([(i,i=1,self % nMat)], mask)

  end subroutine activeIdx


  !!
  !! Recalculates majorant using current active materials
  !!
  subroutine calculateMajorant(self)
    class(isotropicMG), intent(inout)           :: self
    integer(shortInt), dimension(:),allocatable :: activeIdx

    ! Find indexes of active materials
    call self % activeIdx(activeIdx)

    ! Calculate majorantXS
    self % P_majorantXS = maxval( self % P_XSs(:,activeIdx) % totalXS, 2 )

  end subroutine calculateMajorant

end module isotropicMG_class
