module isotropicMG_class

  use numPrecision
  use genericProcedures,     only : fatalError
  use dictionary_class,      only : dictionary
  use IOdictionary_class,    only : IOdictionary

  use perMaterialMgXs_inter, only : perMaterialMgXs
  use outscatterCDF_class,   only : outscatterCDF
  use xsMacroSet_class,      only : xsMacroSet
  use releaseMatrixMG_class, only : releaseMatrixMG

  implicit none
  private

  !!
  !! Small type to store decription of the material
  !!
  type, private :: materialData
    character(nameLen) :: matName
    logical(defBool)   :: isFissile
  end type materialData

  type, public :: isotropicMG
    private
    integer(shortInt)   :: nG
    integer(shortInt)   :: nMat

    type(xsMacroSet),dimension(:,:),pointer            :: XSs => null()  ! (energyGroup,matIdx)
    type(outscatterCDF), dimension(:,:), allocatable   :: transferMatrix ! (energyGroup,matIdx)
    type(outscatterCDF), dimension(:), allocatable     :: chiValues      ! (matIdx)
    type(releaseMatrixMG), dimension(:), allocatable   :: releaseData    ! (matIdx)
    type(materialData), dimension(:), allocatable      :: matData        ! (matIdx

  contains
    procedure  :: init

    ! Procedures to obtain XSs
    !procedure  :: getMatMacroXS
    !procedure  :: getTransXS
    !procedure  :: getMajorantXS
    !procedure  :: getTotalMatXS

    ! Procedures to obtain emission data
    !procedure  :: releaseAt
    !procedure  :: sampleMuGout

    ! Procedures to access material information
    !procedure  :: getMatIdx
    !procedure  :: getMatName
    !procedure  :: isFissileMat

    ! Private procedure
    procedure, private :: readMaterial

  end type isotropicMG

contains

  !!
  !! Load material data from a dictionary
  !!
  subroutine init(self,dict)
    class(isotropicMG), intent(inout) :: self
    type(dictionary), intent(in)      :: dict
    integer(shortInt)                 :: nG,nMat,matIdx

    ! Read number of energy groups and materials
    nG = dict % getInt('numberOfGroups')
    nMat = size( dict % keysDict() )

    self % nG   = nG
    self % nMat = nMat

    ! Allocate space fo XS data
    if (associated( self % XSs           )) deallocate (self % XSs            )
    if (allocated(  self % transferMatrix)) deallocate (self % transferMatrix )
    if (allocated(  self % chiValues     )) deallocate (self % chiValues      )
    if (allocated(  self % releaseData   )) deallocate (self % releaseData    )
    if (allocated(  self % matData       )) deallocate (self % matData        )

    allocate( self % XSs(nG,nMat)           )
    allocate( self % transferMatrix(ng,nMat))
    allocate( self % chiValues(nMat)        )
    allocate( self % releaseData(nMat)      )
    allocate( self % matData(nMat)          )

    ! Read material names
    self % matData(:) % matName = dict % keysDict()

    ! Read individual material data
    do matIdx=1,nMat
      call self % readMaterial(dict % getDict( self % matData(matIdx) % matName ),&
                               matIdx, &
                               self % nG)
    end do

  end subroutine init

  !!
  !! Read material under index "idx" from dictionary "dict"
  !! For group transfer matrixes indexing convention is (G_out,G_in)
  !! When Matrixes are rank1 then the sequence follows SERPENT(2.1.30) output convention:
  !! G_in->1, G_in->2 ...
  !!
  subroutine readMaterial(self,dict,idx, nG)
    class(isotropicMG), intent(inout) :: self
    type(dictionary), intent(in)      :: dict
    integer(shortInt), intent(in)     :: idx
    integer(shortInt), intent(in)     :: nG
    real(defReal),dimension(nG)       :: tempXS
    real(defReal),dimension(nG*nG)    :: tempXSmatrix_rank1
    real(defReal),dimension(nG,nG)    :: tempXSmatrix
    type(IOdictionary)                :: xsDict
    integer(shortInt)                 :: i
    character(100), parameter         :: Here='readMaterial (isotropicMG_class.f90)'

    ! Obtain path from dict and read into xsDict
    call xsDict % initFrom( dict % getChar('xsFile'))

    ! Verify size of stored data
    if (size(xsDict % getRealArray('capture')) /= nG) then
      call fatalError(Here,'capture xs are inconsistant with number of energy groups')

    else if (size(xsDict % getRealArray('fission')) /= nG) then
      call fatalError(Here,'fission xs are inconsistant with number of energy groups')

    elseif (size(xsDict % getRealArray('scattering_P0')) /= nG*nG) then
      call fatalError(Here,'scatter xs are inconsistant with number of energy groups')

    else if (size(xsDict % getRealArray('scatteringMultiplicity')) /= nG*nG) then
      call fatalError(Here,'scattering production data is inconsistant with number of energy groups')

    else if (size(xsDict % getRealArray('chi')) /= nG) then
      call fatalError(Here,'chi data is inconsistant with number of energy groups')

    else if (size(xsDict % getRealArray('nu')) /= nG) then
      call fatalError(Here,'nu data is inconsistant with number of energy groups')

    end if

    ! Load and store capture XSs
    tempXS = xsDict % getRealArray('capture')
    if (any( tempXS < 0.0)) call fatalError(Here,'capture xss are -ve')
    self % XSs(:,idx) % captureXS = tempXS

    ! Load and store fission XSs
    tempXS = xsDict % getRealArray('fission')
    if (any( tempXS < 0.0)) call fatalError(Here,'fission xss are -ve')
    self % XSs(:,idx) % fissionXS = tempXS

    ! Load and store chi values
    tempXS = xsDict % getRealArray('chi')
    if (any( tempXS < 0.0)) call fatalError(Here,'chi is -ve')
    call self % chiValues(idx) % init(tempXS)

    ! Load scattering matrix. Indexing convenction is (G_out,G_in)
    tempXSmatrix_rank1 = xsDict % getRealArray('scattering_P0')
    if (any( tempXSmatrix_rank1 < 0.0)) call fatalError(Here,'Scattering_P0 is -ve')
    tempXSmatrix = reshape(tempXSmatrix_rank1,[nG, nG])
    do i=1,nG
      call self % transferMatrix(i,idx) % init ( tempXSmatrix(:,i) )
    end do
    self % XSs(:,idx) % scatterXS = sum(tempXSmatrix,1)

    ! Load production matrix and nu. Indexing convenction is (G_out,G_in).
    tempXSmatrix_rank1 = xsDict % getRealArray('scatteringMultiplicity')
    if (any( tempXSmatrix_rank1 < 0.0)) call fatalError(Here,'scateringMultiplicity in -ve')
    tempXSmatrix = reshape(tempXSmatrix_rank1,[nG, nG])

    tempXS = xsDict % getRealArray('nu')

    call self % releaseData(idx) % init(tempXS, tempXSmatrix)

    ! Calculate total XS
    self % XSs(:,idx) % totalXS = self % XSs(:,idx) % scatterXS + &
                                  self % XSs(:,idx) % captureXS + &
                                  self % XSs(:,idx) % fissionXS




  end subroutine readMaterial
    
end module isotropicMG_class
