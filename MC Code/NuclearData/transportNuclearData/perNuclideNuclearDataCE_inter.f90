module perNuclideNuclearDataCE_inter

  use numPrecision
  use genericProcedures,          only : fatalError
  use transportNuclearData_inter, only : transportNuclearData
  use particle_class,             only : particle
  use RNG_class,                  only : RNG

  ! XS packages
  use xsMacroSet_class,         only : xsMacroSet_ptr
  use xsMainSet_class,          only : xsMainSet_ptr
  use xsNucMacroSet_class,      only : xsNucMacroSet_ptr

  implicit none
  private

  type, public,extends(transportNuclearData),abstract :: perNuclideNuclearDataCE
    private
  contains
    ! Extend generic procedures to accept calls with energy value
    generic :: getTransXS    => getTransXS_E
    generic :: getMajorantXS => getMajorantXS_E
    generic :: getTotalMatXS => getTotalMatXS_E
    generic :: getMatMacroXS => getMatMacroXS_E

    ! Adapters to transport procedures
    procedure :: getTransXS_p
    procedure :: getMajorantXS_p
    procedure :: getTotalMatXS_p
    procedure :: getMatMacroXS_p

    procedure :: initFissionSite

    ! Per energy access to transport xs data by Energy
    procedure(getTransXS_E),deferred      :: getTransXS_E
    procedure(getMajorantXS_E),deferred   :: getMajorantXS_E
    procedure(getTotalMatXS_E),deferred   :: getTotalMatXS_E
    procedure(getMatMacroXS_E), deferred  :: getMatMacroXS_E
    procedure(initFissionSite_E),deferred :: initFissionSite_E


    ! Procedures to obtain nuclide data
    procedure(getMainNucXS),deferred      :: getMainNucXS
    procedure(xsOf), deferred             :: xsOf
    procedure(invertScattering), deferred :: invertScattering
    procedure(releaseAt),deferred         :: releaseAt
    procedure(sampleMuEout),deferred      :: sampleMuEout
    procedure(sampleMu),deferred          :: sampleMu
    procedure(isInCMframe),deferred       :: isInCMframe
    procedure(isFissileNuc),deferred      :: isFissileNuc
    procedure(getMass),deferred           :: getMass
    procedure(getkT),deferred             :: getkT

    procedure(getNucMacroXS),deferred     :: getNucMacroXS

  end type perNuclideNuclearDataCE

  abstract interface
    !!
    !! Return transport XS of a given material.
    !! Transport XS is used to determine how far neutron flies. For transport-correced
    !! MG calculations it can be diffrent from material total XS. Usually is equal to material
    !! total XS.
    !!
    function getTransXS_E(self,E,matIdx) result (xs)
      import :: shortInt ,&
                defReal  ,&
                perNuclideNuclearDataCE
      class(perNuclideNuclearDataCE), intent(inout) :: self
      real(defReal),intent(in)                      :: E
      integer(shortInt), intent(in)                 :: matIdx
      real(defReal)                                 :: xs
    end function getTransXS_E

    !!
    !! Return majorant XS for all active materials
    !!
    function getMajorantXS_E(self,E) result(xs)
      import :: defReal  ,&
                perNuclideNuclearDataCE
      class(perNuclideNuclearDataCE), intent(inout) :: self
      real(defReal),intent(in)                      :: E
      real(defReal)                                 :: xs
    end function getMajorantXS_E

    !!
    !! Get total XS of a given material
    !!
    function getTotalMatXS_E(self,E,matIdx) result (xs)
      import :: shortInt ,&
                defReal  ,&
                perNuclideNuclearDataCE
      class(perNuclideNuclearDataCE), intent(inout) :: self
      real(defReal),intent(in)                      :: E
      integer(shortInt), intent(in)                 :: matIdx
      real(defReal)                                 :: xs
    end function getTotalMatXS_E

    !!
    !! Get set of material macroscopic xross-sections
    !!
    subroutine getMatMacroXS_E(self,macroXS,E,matIdx)
      import :: xsMacroSet_ptr, &
                defReal, &
                shortInt, &
                perNuclideNuclearDataCE
      class(perNuclideNuclearDataCE), intent(inout) :: self
      type(xsMacroSet_ptr),intent(inout)            :: macroXS
      real(defReal),intent(in)                      :: E
      integer(shortInt),intent(in)                  :: matIdx
    end subroutine getMatMacroXS_E

    !!
    !! Sample fission site knowing that particle is in fissile material and is CE
    !!
    subroutine initFissionSite_E(self,p)
      import :: perNuclideNuclearDataCE, &
                particle
      class(perNuclideNuclearDataCE), intent(in) :: self
      class(particle), intent(inout)             :: p

    end subroutine initFissionSite_E

    !!
    !! Subroutine which attaches a pointer to Main xs set for a given nuclide
    !!
    subroutine getMainNucXS(self,xsPtr,E,nucIdx)
      import :: perNuclideNuclearDataCE, &
                xsMainSet_ptr, &
                defReal, &
                shortInt
      class(perNuclideNuclearDataCE), intent(inout)  :: self
      type(xsMainSet_ptr),intent(inout)              :: xsPtr
      real(defReal), intent(in)                      :: E
      integer(shortInt), intent(in)                  :: nucIdx
    end subroutine getMainNucXS

    !!
    !! Obtain xs for given energy, nuclide and MT number
    !!
    function xsOf(self,E,nucIdx,MT) result(xs)
      import :: perNuclideNuclearDataCE, &
                defReal, &
                shortInt
      class(perNuclideNuclearDataCE), intent(inout) :: self
      real(defReal),intent(in)                      :: E
      integer(shortInt),intent(in)                  :: nucIdx
      integer(shortInt),intent(in)                  :: MT
      real(defReal)                                 :: xs
    end function xsOf

    !!
    !! Invert scattering - for given random number sample MT number of a scattering reaction
    !! All reactions that produce secondary neutrons should be considered
    !!
    function invertScattering(self,E,r,nucIdx) result(MT)
      import :: perNuclideNuclearDataCE, &
                defReal ,&
                shortInt
      class(perNuclideNuclearDataCE), intent(inout) :: self
      real(defReal), intent(in)                     :: E
      real(defReal), intent(in)                     :: r
      integer(shortInt),intent(in)                  :: nucIdx
      integer(shortInt)                             :: MT
    end function invertScattering

    !!
    !! Returns average neutron emission at a given energy
    !!
    function releaseAt(self,E_in,MT,nucIdx) result(nu)
      import :: perNuclideNuclearDataCE, &
                defReal ,&
                shortInt
      class(perNuclideNuclearDataCE), intent(in)  :: self
      real(defReal), intent(in)     :: E_in
      integer(shortInt), intent(in) :: MT
      integer(shortInt), intent(in) :: nucIdx
      real(defReal)                 :: nu
    end function releaseAt

    !!
    !! Sample deflection angle  and emission energy for given MT and nuclide index
    !!
    subroutine sampleMuEout(self,mu,E_out,E_in,rand,MT,nucIdx)
      import :: perNuclideNuclearDataCE, &
                defReal ,&
                shortInt, &
                RNG
      class(perNuclideNuclearDataCE), intent(in)  :: self
      real(defReal), intent(out)    :: mu
      real(defReal), intent(out)    :: E_out
      real(defReal), intent(in)     :: E_in
      class(RNG), intent(inout)     :: rand
      integer(shortInt), intent(in) :: MT
      integer(shortInt), intent(in) :: nucIdx
    end subroutine sampleMuEout

    !!
    !! Sample deflection angle for given MT and nuclide index
    !!
    subroutine sampleMu(self,mu,E_in,rand,MT,nucIdx)
      import :: perNuclideNuclearDataCE, &
                defReal ,&
                shortInt, &
                RNG
      class(perNuclideNuclearDataCE), intent(in)  :: self
      real(defReal), intent(out)                  :: mu
      real(defReal), intent(in)                   :: E_in
      class(RNG), intent(inout)                   :: rand
      integer(shortInt), intent(in)               :: MT
      integer(shortInt), intent(in)               :: nucIdx
    end subroutine sampleMu

    !!
    !! Function which returns .true. if emission data is provided in CM frame for a given MT and
    !! nuclide index.
    !!
    function isInCMframe(self,MT,nucIdx) result(isIt)
      import :: perNuclideNuclearDataCE, &
                defBool ,&
                shortInt
      class(perNuclideNuclearDataCE), intent(in)  :: self
      integer(shortInt), intent(in) :: MT
      integer(shortInt), intent(in) :: nucIdx
      logical(defBool)              :: isIt
    end function isInCMframe

    !!
    !! Returns .true. if nuclide under nucIdx is fissile
    !!
    function isFissileNuc(self,nucIdx) result(isIt)
      import :: perNuclideNuclearDataCE, &
                defBool ,&
                shortInt
      class(perNuclideNuclearDataCE), intent(in)  :: self
      integer(shortInt), intent(in) :: nucIdx
      logical(defBool)              :: isIt
    end function isFissileNuc

    !!
    !! Returns Mass of nuclide in neutron masses
    !!
    function getMass(self,nucIdx) result(A)
      import :: perNuclideNuclearDataCE, &
                defReal ,&
                shortInt
      class(perNuclideNuclearDataCE), intent(in)  :: self
      integer(shortInt),intent(in) :: nucIdx
      real(defReal)                :: A
    end function getMass

    !!
    !! Returns temperature of the nuclide; kT in [MeV]
    !!
    function getkT(self,nucIdx) result(kT)
      import :: perNuclideNuclearDataCE, &
                defReal ,&
                shortInt
      class(perNuclideNuclearDataCE), intent(in)  :: self
      integer(shortInt),intent(in) :: nucIdx
      real(defReal)                :: kT
    end function getKT

    !!
    !! Get set of material nuclide macroscopic cross-sections
    !!
    subroutine getNucMacroXS(self,nucMacroXS,E,matIdx)
      import :: xsNucMacroSet_ptr, &
                defReal, &
                shortInt, &
                perNuclideNuclearDataCE
      class(perNuclideNuclearDataCE), intent(inout) :: self
      type(xsNucMacroSet_ptr),intent(inout)         :: nucMacroXS
      real(defReal),intent(in)                      :: E
      integer(shortInt),intent(in)                  :: matIdx
    end subroutine getNucMacroXS

  end interface

contains

  !!
  !! getTransXS adapter to translate call with particle to a call with energy value
  !! Returns error if multigroup neutron is provided
  !!
  function getTransXS_p(self,p,matIdx) result (xs)
    class(perNuclideNuclearDataCE), intent(inout) :: self
    class(particle), intent(in)                   :: p
    integer(shortInt), intent(in)                 :: matIdx
    real(defReal)                                 :: xs
    character(100), parameter            :: Here='getTransXS_p (perNuclideNuclearDataCE_inter.f90)'

    if (p % isMG) then
      call fatalError(Here,'Multigroup neutron given to CE nuclear Data')
    end if

    xs = self % getTransXS_E(p % E, matIdx)

  end function getTransXS_p

  !!
  !! getMajorantXS adapter to translate call with particle to a call with energy value
  !! Returns error if multigroup neutron is provided
  !!
  function getMajorantXS_p(self,p) result (xs)
    class(perNuclideNuclearDataCE), intent(inout) :: self
    class(particle), intent(in)                   :: p
    real(defReal)                                 :: xs
    character(100), parameter         :: Here='getMajorantXS_p (perNuclideNuclearDataCE_inter.f90)'

    if (p % isMG) then
      call fatalError(Here,'Multigroup neutron given to CE nuclear Data')
    end if

    xs = self % getMajorantXS_E(p % E)

  end function getMajorantXS_p

  !!
  !! getTotalMatXS adapter to translate call with particle to a call with energy value
  !! Returns error if multigroup neutron is provided
  !!
  function getTotalMatXS_p(self,p,matIdx) result (xs)
    class(perNuclideNuclearDataCE), intent(inout) :: self
    class(particle), intent(in)                   :: p
    integer(shortInt), intent(in)                 :: matIdx
    real(defReal)                                 :: xs
    character(100), parameter         :: Here='getTotalMatXS_p (perNuclideNuclearDataCE_inter.f90)'

    if (p % isMG) then
      call fatalError(Here,'Multigroup neutron given to CE nuclear Data')
    end if

    xs = self % getTotalMatXS_E(p % E, matIdx)

  end function getTotalMatXS_p

  !!
  !! getMatMacroXS adapter to translate call with particle to a call with energy value
  !! Returns error if multigroup neutron is provided
  !!
  subroutine getMatMacroXS_p(self,macroXS,p,matIdx)
    class(perNuclideNuclearDataCE), intent(inout) :: self
    type(xsMacroSet_ptr),intent(inout)            :: macroXS
    class(particle), intent(in)                   :: p
    integer(shortInt),intent(in)                  :: matIdx
    character(100), parameter         :: Here='getMatMacroXS_p (perNuclideNuclearDataCE_inter.f90)'

    if (p % isMG) then
      call fatalError(Here,'Multigroup neutron given to CE nuclear Data')
    end if

    call self % getMatMacroXS_E(macroXS, p % E, matIdx)

  end subroutine getMatMacroXS_p

  !!
  !! Function to generate a fission site from a fissile material.
  !! Necassary in initialisation of an eigenvalue calculation:
  !! Is it to avoid reapeting code to determine if material is fissile and CE
  !!
  subroutine initFissionSite(self,p)
    class(perNuclideNuclearDataCE), intent(in)   :: self
    class(particle), intent(inout)               :: p
    integer(shortInt)                            :: matIdx
    character(100), parameter      :: Here = 'initFissionSite (perNuclideNuclearDataCE_inter.f90)'

    matIdx = p % matIdx()

    ! Determine if material is fissile
    if ( .not.self % isFissileMat(matIdx) ) then
      p % isDead = .true.
      return
    end if

    ! Determine if particle is CE and throw error if it is MG
    if ( p % isMG ) then
      call fatalError(Here,'Trying to create fission site for MG particle with CE data')
    end if

    ! Call specific instance implementation of sampling
    call self % initFissionSite_E(p)

  end subroutine initFissionSite


end module perNuclideNuclearDataCE_inter
