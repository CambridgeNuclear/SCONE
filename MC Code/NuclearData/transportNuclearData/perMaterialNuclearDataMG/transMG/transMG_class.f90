module transMG_class

  use numPrecision
  use genericProcedures,  only : fatalError, numToChar
  use dictionary_class,   only : dictionary
  use IOdictionary_class, only : IOdictionary
  use isotropicMG_class,  only : isotropicMG, init_super         => init,        &
                                              readMaterial_super => readMaterial,&
                                              activeIdx_super    => activeIdx

  implicit none
  private

  !!
  !! Extension of isotropic MG scattering, but with transport correction on the transport XS.
  !! pone stands for P1
  !!
  type, public,extends(isotropicMG) :: transMG
    private
    real(defReal),dimension(:,:), allocatable :: transportXSs
  contains
    ! Overwrite key isotropicMG procedures
    procedure         :: init
    procedure,private :: readMaterial
    procedure         :: calculateMajorant
    procedure         :: getTransXS_G

  end type transMG

contains

  !!
  !! Initialise transMG
  !!
  subroutine init(self, dict)
    class(transMG), intent(inout)  :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt)             :: nG, nMat
    character(nameLen),dimension(:),allocatable :: matNames

    ! Read material names
    call dict % keysDict(matNames)

    ! Read number of energy groups and materials
    call dict % get(nG, 'numberOfGroups')
    nMat = size( matNames )

    ! Allocate space for transport XSs
    allocate(self % transportXSs(nG,nMat))

    ! Call superclass initialisation (note that it will call extended readMaterial in this module)
    call init_super(self,dict)

    print *, self % transportXSs(:,1)
    print *, self % transportXSs(:,2)

  end subroutine init

  !!
  !! Extend superclass subroutine to read additional data to produce P1 corrected
  !! transport XSs
  !!
  subroutine readMaterial(self, dict, idx, nG)
    class(transMG), intent(inout)           :: self
    class(dictionary), intent(in)          :: dict
    integer(shortInt), intent(in)          :: idx
    integer(shortInt), intent(in)          :: nG
    type(IOdictionary)                     :: xsDict
    integer(shortInt)                      :: i
    real(defReal),dimension(nG)            :: tempXS
    real(defReal),dimension(:),allocatable :: tempXSmatrix_rank1
    character(100), parameter :: Here='readMaterial (transMG_class.f90)'

    ! Call superclass procedure
    call readMaterial_super(self, dict, idx, nG)

    ! Obtain path from dict and read into xsDict
    call xsDict % initFrom( dict % getChar('xsFile'))

    ! Obtain array of sum of P1 outscattering XSs
    call xsDict % get(tempXSmatrix_rank1, 'scattering_P1')

    ! Verify size of the matrix
    if (size(tempXSmatrix_rank1) /= nG*nG) then
      call fatalError(Here,'scattering_P1 matrix contains only:' // &
                            numToChar(size(tempXSmatrix_rank1))// ' elements')
    end if

    ! Calculate transport corrected XSs
    do i=1,nG
      tempXS(i) = self % getTotalMatXS_G(i,idx)
    end do
    tempXS = tempXS - sum(reshape(tempXSmatrix_rank1, [nG, nG]), 1)

    ! Check if -ve transport corrected XS was generated
    if (any(tempXS < ZERO)) call fatalError(Here,'-ve Transport corrected XS in material '// &
                                                 self % getName(idx) )

    ! Load Trans XS
    self % transportXSs(:,idx) = tempXS

  end subroutine readMaterial

  !!
  !! Recalculates majorant using current active materials
  !!
  subroutine calculateMajorant(self)
    class(transMG), intent(inout)                :: self
    integer(shortInt), dimension(:),allocatable :: activeIdx

    ! Find indexes of active materials
    call activeIdx_super(self, activeIdx)

    ! Calculate majorantXS
    self % P_majorantXS = maxval( self % transportXSs(:,activeIdx), 2 )

  end subroutine calculateMajorant

  !!
  !! Return transport XS
  !!
  function getTransXS_G(self,G,matIdx) result(xs)
    class(transMG), intent(inout)       :: self
    integer(shortInt), intent(in)      :: G
    integer(shortInt), intent(in)      :: matIdx
    real(defReal)                      :: xs

    xs = self % transportXSs(G,matIdx)

  end function getTransXS_G

end module transMG_class
