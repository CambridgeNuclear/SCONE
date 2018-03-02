module byNucNoMT_class

  ! Global modules
  use numPrecision
  use genericProcedures,        only : fatalError
  use RNG_class,                only : RNG

  ! Modules specific to byNucNoMT type of XS data storage
  use byNucNoMT_Data_class ,    only : byNucNoMT_Data
  use nuclideMemoryNoMT_class,  only : nuclideMemoryNoMT
  use materialMemoryNoMT_class, only : materialMemoryNoMT
  use materialDataNoMT_class,   only : materialDataNoMT
  use aceNoMT_class,            only : aceNoMT

  ! Cross-section packages to interface with Collision Operator
  use xsMainCDF_class,          only : xsMainCDF
  use xsMainSet_class,          only : xsMainSet
  use matNucCDF_class,          only : matNucCDF
  use xsMacroSet_class,         only : xsMacroSet

  implicit none
  private

  type, public :: byNucNoMT
    !private **** Should be uncommented after debug

    ! Storage of static XS and Material data
    type(byNucNoMT_Data),pointer                    :: dataBlock => null()

    ! Dynamic, local handels to store interpolated data and have memory of the last call
    type(nuclideMemoryNoMT), dimension(:), pointer  :: nucShelf  => null()
    type(materialMemoryNoMT), dimension(:),pointer  :: matShelf  => null()

    ! Data for majorant XS calculations
    integer(shortInt), dimension(:), allocatable :: activeMaterials   !! Materials used in majorant calculation
    real(defReal)                                :: majE       = -1.0 !! Energy of current majorant
    real(defReal)                                :: majorantXS = -1.0 !! Current majorant

  contains
    procedure :: readFrom
    ! Procedures to access nuclide data (microscopic xss)
    procedure :: getMainNucCDF
    procedure :: getMainNucXS
    procedure :: sampleMuEout
    procedure :: releaseAt
    procedure :: isInCMframe

    ! Procedures to access material data (Macroscopic XSs)
    procedure :: getTotalMatXS
    procedure :: getMatNucCDF
    procedure :: getMatMacroXS
    procedure :: getMajorantXS

  end type byNucNoMT

contains

  !!
  !! Read material and nuclide data using input files at the provided paths
  !!
  subroutine readFrom(self,matInput,nuclideLib)
    class(byNucNoMT),intent(inout)       :: self
    character(*), intent(in)             :: matInput
    character(*), intent(in)             :: nuclideLib
    type(aceNoMT),pointer                :: nucPtr
    type(materialDataNoMt),pointer       :: matPtr
    integer(shortInt)                    :: numNuclide, numMaterials
    integer(shortInt)                    :: i
    character(100), parameter            :: Here = 'readFrom (byNUcNoMT_class.f90)'

    ! Check if material data was alrady read. If it was return error becouse prodecures
    ! for cleaning memory from material and xs data are not implemented
    if(associated(self % dataBlock)) call fatalError(Here,'It is forbidden to reinitialise XS data')

    ! Read Material data into a shared "fat" object
    allocate (self % dataBlock)
    call self % dataBlock  % readFrom(matInput, nuclideLib)

    ! Allocate space for nuclide and material shelfs
    numNuclide   = size(self % dataBlock % nucXsData)
    numMaterials = size(self % dataBlock % matData  )

    allocate (self % nucShelf (numNuclide  ))
    allocate (self % matShelf (numMaterials))

    ! Attach nuclides to the shelf
    do i=1,numNuclide
      nucPtr => self % dataBlock % nucXsData(i)
      call self % nucShelf(i) % init(i,nucPtr)

    end do

    ! Attach materials to the shelf
    do i=1,numMaterials
      matPtr => self % dataBlock % matData(i)
      call self % matShelf(i) % init(matPtr,self % nucShelf)
    end do

    ! At this point assume all defined materials are present in the geometry
    ! Include all material indexes
    self % activeMaterials = [(i , i=1,numMaterials)]

  end subroutine readFrom

  !!
  !! Subroutine to attach pointer to CDF for the main reaction channel
  !!
  subroutine getMainNucCDF(self,cdfPtr,E,nucIdx)
    class(byNucNoMT), intent(inout)         :: self
    type(xsMainCDF),pointer,intent(inout)   :: cdfPtr
    real(defReal)                           :: E
    integer(shortInt)                       :: nucIdx

    ! Set approperiate nuclide shelf to energy
    call self % nucShelf(nucIdx) % setEnergy(E)

    ! Point to interpolated cdf
    cdfPtr => self % nucShelf(nucIdx) % mainCDF

  end subroutine getMainNucCDF

  !!
  !! Subroutine to attach pointer to the Main xs set of the nuclide
  !!
  subroutine getMainNucXS(self,xsPtr,E,nucIdx)
    class(byNucNoMT), intent(inout)         :: self
    type(xsMainSet),pointer,intent(inout)   :: xsPtr
    real(defReal)                           :: E
    integer(shortInt)                       :: nucIdx

    ! Set approperiate nuclide shelf to energy
    call self % nucShelf(nucIdx) % setEnergy(E)

    ! Point to interpolated xs set
    xsPtr => self % nucShelf(nucIdx) % xs

  end subroutine getMainNucXS

  !!
  !! Subroutine which samples deflecton angle and emission energy for a given MT and nuclide index.
  !!
  subroutine sampleMuEout(self,mu,E_out,E_in,rand,MT,nucIdx)
    class(byNucNoMT), intent(in)  :: self
    real(defReal), intent(out)    :: mu
    real(defReal), intent(out)    :: E_out
    real(defReal), intent(in)     :: E_in
    class(RNG), intent(inout)     :: rand
    integer(shortInt), intent(in) :: MT
    integer(shortInt), intent(in) :: nucIdx

    call self % dataBlock % nucXsData(nucIdx) % sampleMuEout( mu, E_out, E_in, rand, MT )

  end subroutine sampleMuEout

  !!
  !! Function which returns average neutron emission at a given energy
  !!
  function releaseAt(self,E_in,MT,nucIdx) result(nu)
    class(byNucNoMT), intent(in)  :: self
    real(defReal), intent(in)     :: E_in
    integer(shortInt), intent(in) :: MT
    integer(shortInt), intent(in) :: nucIdx
    real(defReal)                 :: nu

    nu = self % dataBlock % nucXsData(nucIdx) % releaseAt(E_in,MT)

  end function releaseAt

  !!
  !! Function which returns .true. if emission data is provided in CM frame for a given MT and
  !! nuclide index.
  !!
  function isInCMframe(self,MT,nucIdx) result(isIt)
    class(byNucNoMT), intent(in)  :: self
    integer(shortInt), intent(in) :: MT
    integer(shortInt), intent(in) :: nucIdx
    logical(defBool)              :: isIt

    isIt = self % dataBlock % nucXSData(nucIdx) % isInCMframe(MT)

  end function isInCMframe

  !!
  !! Function to obtain total XS for material identified by its index
  !!
  function getTotalMatXS(self,E,matIdx) result (xs)
    class(byNucNoMT), intent(inout)  :: self
    real(defReal),intent(in)         :: E
    integer(shortInt), intent(in)    :: matIdx
    real(defReal)                    :: xs

    xs = self % matShelf(matIdx) % getTotal(E)

  end function getTotalMatXS


  !!
  !! Function to obtain majorant XS for all active materials
  !!
  function getMajorantXS(self,E) result (xs)
    class(byNucNoMT), intent(inout)   :: self
    real(defReal), intent(in)         :: E
    real(defReal)                     :: xs
    integer(shortInt)                 :: i, currentMat

    if (self % majE == E) then
      xs = self % majorantXS

    else
      ! Calculate new majorant XS
      xs = 0.0
      do i=1,size(self % activeMaterials)
        currentMat = self % activeMaterials(i)
        xs = max(xs, self % getTotalMatXS(E,currentMat) )

      end do
    end if

  end function getMajorantXS


  !!
  !! Subroutine to attach pointer to a material's cdf to choose collision nuclide
  !!
  subroutine getMatNucCDF(self,nucCDF,E,matIdx)
    class(byNucNoMT),intent(inout)        :: self
    type(matNucCDF),pointer,intent(inout) :: nucCDF
    real(defReal),intent(in)              :: E
    integer(shortInt),intent(in)          :: matIdx

    call self % matShelf(matIdx) % setTotalToEnergy(E)
    nucCDF => self % matShelf(matIdx) % nucCDF

  end subroutine getMatNucCDF

  !!
  !! Subroutine to attach pointer to a material's macroscopic XS set
  !!
  subroutine getMatMacroXS(self,macroXS,E,matIdx)
    class(byNucNoMT), intent(inout)        :: self
    type(xsMacroSet),pointer,intent(inout) :: macroXS
    real(defReal),intent(in)               :: E
    integer(shortInt),intent(in)           :: matIdx

    call self % matShelf(matIdx) % setEnergy(E)
    macroXS => self % matShelf(matIdx) % XS

  end subroutine



end module byNucNoMT_class
