module byNucNoMT_class

  use numPrecision
  use genericProcedures,       only : fatalError
  use byNucNoMT_Data_class ,   only : byNucNoMT_Data
  use nuclideMemoryNoMT_class, only : nuclideMemoryNoMT
  use aceNoMT_class,           only : aceNoMT

  use xsMainCDF_class,         only : xsMainCDF
  use xsMainSet_class,         only : xsMainSet
  implicit none
  private

  type, public :: byNucNoMT
    private
    ! Material Data Handle State Variabels
    real(defReal)            ::  E
    integer(shortInt)        ::  material
    integer(shortInt)        ::  isotope
    integer(shortInt)        ::  MT

    ! Material Data Handle storage
    type(nuclideMemoryNoMT), dimension(:), pointer        :: nucShelf  => null()
    real(defReal)                                         :: majorant
    type(byNucNoMT_Data),pointer                          :: dataBlock => null()

  contains
    procedure :: readFrom
    procedure :: getMainCDF
    procedure :: getMainXS
  !  procedure :: initDebug
  !  procedure :: init
    ! Get majorant XS
    ! Get material XS
    ! Get pointer to isotope total XS and pointer to isotope index map
    ! Get neutron emission given nuclide index and MT
    ! Sample deflection mu and secondary energy with CF flag given random number generator
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
    numMaterials = size(self % dataBlock % matNames)

    allocate (self % nucShelf (numNuclide))
    ! *** Allocate material shelf *** !



    ! Attach nuclides to the shelf
    do i=1,numNuclide
      nucPtr => self % dataBlock % nucXsData(i)
      call self % nucShelf(i) % init(i,nucPtr)
    end do

  end subroutine readFrom

  !!
  !! Subroutine to attach pointer to CDF for the main reaction channel
  !!
  subroutine getMainCDF(self,cdfPtr,E,nucIdx)
    class(byNucNoMT), intent(inout)         :: self
    class(xsMainCDF),pointer,intent(inout)  :: cdfPtr
    real(defReal)                           :: E
    integer(shortInt)                       :: nucIdx

    ! Set approperiate nuclide shelf to energy
    call self % nucShelf(nucIdx) % setEnergy(E)

    ! Point to interpolated cdf
    cdfPtr => self % nucShelf(nucIdx) % mainCDF

  end subroutine getMainCDF

  !!
  !! Subroutine to attach pointer to the Main xs set of the nuclide
  !!
  subroutine getMainXS(self,xsPtr,E,nucIdx)
    class(byNucNoMT), intent(inout)         :: self
    class(xsMainSet),pointer,intent(inout)  :: xsPtr
    real(defReal)                           :: E
    integer(shortInt)                       :: nucIdx

    ! Set approperiate nuclide shelf to energy
    call self % nucShelf(nucIdx) % setEnergy(E)

    ! Point to interpolated xs set
    xsPtr => self % nucShelf(nucIdx) % xs

  end subroutine getMainXS


end module byNucNoMT_class
