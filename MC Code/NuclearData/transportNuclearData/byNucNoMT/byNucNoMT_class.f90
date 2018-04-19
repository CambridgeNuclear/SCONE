module byNucNoMT_class

  ! Global modules
  use numPrecision
  use genericProcedures,             only : fatalError, linFind, searchError
  use RNG_class,                     only : RNG
  use dictionary_class,              only : dictionary
  use particle_class,                only : particle

  ! Interface
  use perNuclideNuclearDataCE_inter, only : perNuclideNuclearDataCE

  ! Modules specific to byNucNoMT type of XS data storage
  use byNucNoMT_Data_class ,         only : byNucNoMT_Data
  use nuclideMemoryNoMT_class,       only : nuclideMemoryNoMT
  use materialMemoryNoMT_class,      only : materialMemoryNoMT
  use materialDataNoMT_class,        only : materialDataNoMT
  use aceNoMT_class,                 only : aceNoMT

  ! Cross-section packages to interface with Collision Operator
 ! use xsMainCDF_class,               only : xsMainCDF
  use xsMainSet_class,               only : xsMainSet_ptr
  use xsNucMacroSet_class,           only : xsNucMacroSet_ptr
  use xsMacroSet_class,              only : xsMacroSet_ptr

  implicit none
  private

  type, public,extends(perNuclideNuclearDataCE) :: byNucNoMT
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
    procedure :: init
    procedure :: readFrom

    ! Nuclear Data Interface Procedures
    procedure :: getIdx
    procedure :: getName

    ! transportNuclearData Interface Procedures
    procedure :: getTransXS_E
    procedure :: getMajorantXS_E
    procedure :: getTotalMatXS_E
    procedure :: isFissileMat
    procedure :: initFissionSite
    procedure :: setActiveMaterials

    ! perNuclideNuclearDataCE procedures
    procedure :: getMainNucXS
    procedure :: xsOf
    procedure :: invertScattering
    procedure :: sampleMuEout
    procedure :: sampleMu
    procedure :: releaseAt
    procedure :: isInCMframe
    procedure :: isFissileNuc
    procedure :: getMass
    procedure :: getkT
    procedure :: getNucMacroXS

    procedure :: getMatMacroXS

    ! Type specific procedures

  end type byNucNoMT

contains
  subroutine init(self,dict)
    class(byNucNoMT), intent(inout)   :: self
    type(dictionary), intent(inout)   :: dict
    type(aceNoMT),pointer                :: nucPtr
    type(materialDataNoMt),pointer       :: matPtr
    integer(shortInt)                    :: numNuclide, numMaterials
    integer(shortInt)                    :: i
    character(100), parameter            :: Here = 'init (byNUcNoMT_class.f90)'

    ! Check if material data was alrady read. If it was return error becouse prodecures
    ! for cleaning memory from material and xs data are not implemented
    if(associated(self % dataBlock)) call fatalError(Here,'It is forbidden to reinitialise XS data')

    ! Read Material data into a shared "fat" object
    allocate(self % dataBlock)
    call self % dataBlock % init(dict)

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


  end subroutine init


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
  !! Returns material index for given material name
  !! Throws error if material is not found
  !!
  function getIdx(self,matName) result(matIdx)
    class(byNucNoMT), intent(in) :: self
    character(*), intent(in)     :: matName
    integer(shortInt)            :: matIdx
    character(100), parameter    :: Here = 'getIdx (byNucNoMt_class.f90)'
    matIdx = linFind(self % dataBlock % matData(:) % name, matName)
    call searchError(matIdx,Here)

  end function getIdx

  !!
  !! Returns material name for given material index
  !! Throws error if material index does not correspond to valid material
  !!
  function getName(self,matIdx) result(matName)
    class(byNucNoMT), intent(in)   :: self
    integer(shortInt), intent(in)  :: matIdx
    character(nameLen)             :: matName
    character(100),parameter       :: Here = 'getName (byNucNoMT_class.f90)'

    associate (matData => self % dataBlock % matData)

      if ( matIdx < 1 .or. matIdx > size(matData) ) then
        call fatalError(Here,'material index is out of bounds')
      end if

      matName = matData(matIdx) % name

    end associate
  end function getName

  !!
  !! Function to obtain transport XS for material given by index
  !! In this implementation it is equal to total XS
  !!
  function getTransXS_E(self,E,matIdx) result (xs)
    class(byNucNoMT), intent(inout)  :: self
    real(defReal),intent(in)         :: E
    integer(shortInt), intent(in)    :: matIdx
    real(defReal)                    :: xs

    xs = self % matShelf(matIdx) % getTotal(E)

  end function getTransXS_E

  !!
  !! Function to obtain majorant XS for all active materials
  !!
  function getMajorantXS_E(self,E) result (xs)
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

  end function getMajorantXS_E

  !!
  !! Function to obtain total XS for material identified by its index
  !!
  function getTotalMatXS_E(self,E,matIdx) result (xs)
    class(byNucNoMT), intent(inout)  :: self
    real(defReal),intent(in)         :: E
    integer(shortInt), intent(in)    :: matIdx
    real(defReal)                    :: xs

    xs = self % matShelf(matIdx) % getTotal(E)

  end function getTotalMatXS_E

  !!
  !! Returns .true. if material contains fissile nuclides
  !!
  function isFissileMat(self,matIdx) result(isIt)
    class(byNucNoMT), intent(in)   :: self
    integer(shortInt), intent(in)  :: matIdx
    logical(defBool)               :: isIt

    isIt = self % dataBlock % matData(matIdx) % isFissile

  end function isFissileMat


  !!
  !! Function to generate a fission site from a fissile material.
  !! Necassary in initialisation of an eigenvalue calculation:
  !!
  subroutine initFissionSite(self,p)
    class(byNucNoMT), intent(in)   :: self
    class(particle), intent(inout) :: p

  end subroutine initFissionSite


  !!
  !! Set all materials that are present in geometry
  !! Active materials will be included in evaluation of majorant XS
  !!
  !! Provide an array of all active material indexes
  !! Order is not significant (assume there is none)
  !!
  subroutine setActiveMaterials(self,matIdxList)
    class(byNucNoMT), intent(inout)           :: self
    integer(shortInt),dimension(:),intent(in) :: matIdxList
    character(100),parameter                  :: Here= 'setActiveMaterials (byNucNoMT_class.f90)'

    if( allocated(self % activeMaterials)) deallocate( self % activeMaterials)

    ! Check if all active material indices are valid
    if ( any( matIdxList < 1 .or. matIdxList > size( self % matShelf) )) then
      call fatalError(Here,'List of active materials contains out of bounds matIdx ')

    end if

    self % activeMaterials = matIdxList

  end subroutine setActiveMaterials

  !!
  !! Subroutine to attach pointer to the Main xs set of the nuclide
  !!
  subroutine getMainNucXS(self,xsPtr,E,nucIdx)
    class(byNucNoMT), intent(inout)         :: self
    type(xsMainSet_ptr),intent(inout)       :: xsPtr
    real(defReal),intent(in)                :: E
    integer(shortInt), intent(in)           :: nucIdx

    ! Set approperiate nuclide shelf to energy
    call self % nucShelf(nucIdx) % setEnergy(E)

    ! Point to interpolated xs set
    xsPtr = self % nucShelf(nucIdx) % xs

  end subroutine getMainNucXS

  !!
  !! Obtain xs for given energy, nuclide and MT number
  !!
  function xsOf(self,E,nucIdx,MT) result(xs)
    class(byNucNoMT), intent(inout) :: self
    real(defReal),intent(in)        :: E
    integer(shortInt),intent(in)    :: nucIdx
     integer(shortInt),intent(in)   :: MT
    real(defReal)                   :: xs
  end function xsOf

  !!
  !! Invert scattering - for given random number sample MT number of a scattering reaction
  !! All reactions that produce secondary neutrons should be considered
  !!
  function invertScattering(self,E,r) result(MT)
    class(byNucNoMT), intent(inout) :: self
    real(defReal), intent(in)       :: E
    real(defReal), intent(in)       :: r
    integer(shortInt)               :: MT
  end function invertScattering

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
  !! Subroutine to sample mu when it is known that reaction is elastic
  !! For now there should be no performance benefit but better implementation may arrive
  !!
  subroutine sampleMu(self,mu,E_in,rand,MT,nucIdx)
    class(byNucNoMT), intent(in)  :: self
    real(defReal), intent(out)    :: mu
    real(defReal), intent(in)     :: E_in
    class(RNG), intent(inout)     :: rand
    integer(shortInt), intent(in) :: MT
    integer(shortInt), intent(in) :: nucIdx
    real(defReal)                 :: dummyE

    call self % dataBlock % nucXSData(nucIdx) % sampleMuEout( mu, dummyE, E_in, rand, MT)

  end subroutine sampleMu

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
  !! Returns .true. if nuclide under nucIdx is fissile
  !!
  function isFissileNuc(self,nucIdx) result(isIt)
    class(byNucNoMT), intent(in)  :: self
    integer(shortInt), intent(in) :: nucIdx
    logical(defBool)              :: isIt

    isIt = self % dataBlock % nucXSData(nucIdx) % isFissile

  end function isFissileNuc

  !!
  !! Function that returns Mass of the nuclide; A in [Mn - neutron mass ]
  !!
  function getMass(self,nucIdx) result(A)
    class(byNucNoMT),intent(in)  :: self
    integer(shortInt),intent(in) :: nucIdx
    real(defReal)                :: A

    A = self % dataBlock % nucXsData(nucIdx) % atomWeight

  end function getMass



  !!
  !! Function that returns temperature of the nuclide; kT in [MeV]
  !!
  function getkT(self,nucIdx) result(kT)
    class(byNucNoMT),intent(in)  :: self
    integer(shortInt),intent(in) :: nucIdx
    real(defReal)                :: kT

    kT = self % dataBlock % nucXsData(nucIdx) % temp

  end function getkT

  !!
  !! Subroutine to attach pointer to a material's cdf to choose collision nuclide
  !!
  subroutine getNucMacroXS(self,nucMacroXs,E,matIdx)
    class(byNucNoMT),intent(inout)            :: self
    type(xsNucMacroSet_ptr),intent(inout)     :: nucMacroXs
    real(defReal),intent(in)                  :: E
    integer(shortInt),intent(in)              :: matIdx

    call self % matShelf(matIdx) % setTotalToEnergy(E)
    nucMacroXs = self % matShelf(matIdx) % nucCDF

  end subroutine getNucMacroXS

  !!
  !! Subroutine to attach pointer to a material's macroscopic XS set
  !!
  subroutine getMatMacroXS(self,macroXS,E,matIdx)
    class(byNucNoMT), intent(inout)        :: self
    type(xsMacroSet_ptr),intent(inout)     :: macroXS
    real(defReal),intent(in)               :: E
    integer(shortInt),intent(in)           :: matIdx

    call self % matShelf(matIdx) % setEnergy(E)
    macroXS = self % matShelf(matIdx) % XS

  end subroutine getMatMacroXS


end module byNucNoMT_class
