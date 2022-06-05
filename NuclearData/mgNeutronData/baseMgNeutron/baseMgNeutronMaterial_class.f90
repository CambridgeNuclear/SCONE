module baseMgNeutronMaterial_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError, numToChar
  use RNG_class,         only : RNG
  use dictionary_class,  only : dictionary
  use dictDeck_class,    only : dictDeck

  ! Nuclear Data Interfaces
  use materialHandle_inter,    only : materialHandle
  use mgNeutronMaterial_inter, only : mgNeutronMaterial, kill_super => kill
  use neutronXSPackages_class, only : neutronMacroXSs

  ! Reaction objects
  use reactionMG_inter,        only : reactionMG
  use fissionMG_class,         only : fissionMG
  use multiScatterMG_class,    only : multiScatterMG
  use multiScatterP1MG_class,  only : multiScatterP1MG

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: baseMgNeutronMaterial_TptrCast
  public :: baseMgNeutronMaterial_CptrCast

  ! Public data location parameters
  ! Use them if accessing data entries directly
  integer(shortInt), parameter, public :: TOTAL_XS      = 1
  integer(shortInt), parameter, public :: IESCATTER_XS  = 2
  integer(shortInt), parameter, public :: CAPTURE_XS    = 3
  integer(shortInt), parameter, public :: FISSION_XS    = 4
  integer(shortInt), parameter, public :: NU_FISSION    = 5

  !!
  !! Basic type of MG material data
  !!
  !! Stores MG data in a table.
  !! Fission is treated as a seperate reaction
  !! All other scattering reactions are lumped into single multiplicative scattering,
  !! which is stored as INELASTIC scatering in macroXSs package! After all it is inelastic in
  !! the sense that outgoing group can change. Diffrent types of multiplicative scattering can be
  !! build. See doc of "init" procedure for details.
  !!
  !! Public members:
  !!   data -> Rank 2 array with all XSs data
  !!
  !! Interface:
  !!   materialHandle interface
  !!   mgNeutronMaterial interface
  !!   init -> initialise Basic MG Material from dictionary and config keyword
  !!   nGroups -> returns number of energy groups
  !!
  !! Note:
  !!   Order of "data" array is: data(XS_type, Group #)
  !!   Dictionary with data must contain following entries:
  !!     -> numberOfGroups
  !!     -> capture [nGx1]
  !!     -> scatteringMultiplicity [nGxnG]
  !!     -> P0 [nGxnG]
  !!   Optional entries:
  !!     -> fission [nGx1]
  !!     -> nu [nGx1]
  !!     -> chi [nGx1]
  !!     -> P# [nGxnG]
  !!
  type, public, extends(mgNeutronMaterial) :: baseMgNeutronMaterial
    real(defReal),dimension(:,:), allocatable :: data
    class(multiScatterMG), allocatable        :: scatter
    type(fissionMG), allocatable              :: fission
    integer(shortInt)                         :: nG

  contains
    ! Superclass procedures
    procedure :: kill
    procedure :: getMacroXSs_byG
    procedure :: getTotalXS
    procedure :: getNuFissionXS
    procedure :: getFissionXS
    procedure :: getChi
    procedure :: getScatterXS

    ! Local procedures
    procedure :: init
    procedure :: nGroups
    procedure :: getTotalPtr
    procedure :: getNuFissionPtr
    procedure :: getChiPtr
    procedure :: getScatterPtr

  end type baseMgNeutronMaterial

contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(baseMgNeutronMaterial), intent(inout) :: self

    ! Call superclass procedure
    call kill_super(self)

    ! Kill local content
    if(allocated(self % data))    deallocate(self % data)
    if(allocated(self % scatter)) deallocate(self % scatter)
    if(allocated(self % fission)) deallocate(self % fission)

  end subroutine kill

  !!
  !! Load Macroscopic XSs into the provided package for a given group index G
  !!
  !! See mgNeutronMaterial documentation for more details
  !!
  subroutine getMacroXSs_byG(self, xss, G, rand)
    class(baseMgNeutronMaterial), intent(in) :: self
    type(neutronMacroXSs), intent(out)       :: xss
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    character(100), parameter :: Here = ' getMacroXSs (baseMgNeutronMaterial_class.f90)'

    ! Verify bounds
    if(G < 1 .or. self % nGroups() < G) then
      call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                           ' Data has only: ' // numToChar(self % nGroups()))
    end if

    ! Get XSs
    xss % total            = self % data(TOTAL_XS, G)
    xss % elasticScatter   = ZERO
    xss % inelasticScatter = self % data(IESCATTER_XS, G)
    xss % capture          = self % data(CAPTURE_XS, G)

    if(self % isFissile()) then
      xss % fission        = self % data(FISSION_XS, G)
      xss % nuFission      = self % data(NU_FISSION, G)
    else
      xss % fission        = ZERO
      xss % nuFission      = ZERO
    end if

  end subroutine getMacroXSs_byG

  !!
  !! Return Total XSs for energy group G
  !!
  !! See mgNeutronMaterial documentationfor details
  !!
  function getTotalXS(self, G, rand) result(xs)
    class(baseMgNeutronMaterial), intent(in) :: self
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    real(defReal)                            :: xs
    character(100), parameter :: Here = ' getTotalXS (baseMgNeutronMaterial_class.f90)'

    ! Verify bounds
    if (G < 1 .or. self % nGroups() < G) then
      call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                           ' Data has only: ' // numToChar(self % nGroups()))
      xs = ZERO ! Avoid warning
    end if

    xs = self % data(TOTAL_XS, G)

  end function getTotalXS

  !!
  !! Return NuFission XS for energy group G
  !!
  !! See mgNeutronMaterial documentationfor details
  !!
  function getNuFissionXS(self, G, rand) result(xs)
    class(baseMgNeutronMaterial), intent(in) :: self
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    real(defReal)                            :: xs
    character(100), parameter :: Here = ' getNuFissionXS (baseMgNeutronMaterial_class.f90)'

    ! Verify bounds
    if (self % isFissile()) then
      if(G < 1 .or. self % nGroups() < G) then
        call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                             ' Data has only: ' // numToChar(self % nGroups()))
        xs = ZERO ! Avoid warning
      end if
      xs = self % data(NU_FISSION, G)
    else
      xs = ZERO
    end if

  end function getNuFissionXS

  !!
  !! Return Fission XS for energy group G
  !!
  !! See mgNeutronMaterial documentationfor details
  !!
  function getFissionXS(self, G, rand) result(xs)
    class(baseMgNeutronMaterial), intent(in) :: self
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    real(defReal)                            :: xs
    character(100), parameter :: Here = ' getFissionXS (baseMgNeutronMaterial_class.f90)'

    ! Verify bounds
    if (self % isFissile()) then
      if(G < 1 .or. self % nGroups() < G) then
        call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                             ' Data has only: ' // numToChar(self % nGroups()))
        xs = ZERO ! Avoid warning
      end if
      xs = self % data(FISSION_XS, G)
    else
      xs = ZERO
    end if

  end function getFissionXS

  !!
  !! Return chi for energy group G
  !!
  !! See mgNeutronMaterial documentationfor details
  !!
  function getChi(self, G, rand) result(chi)
    class(baseMgNeutronMaterial), intent(in) :: self
    integer(shortInt), intent(in)            :: G
    class(RNG), intent(inout)                :: rand
    real(defReal)                            :: chi
    character(100), parameter :: Here = ' getChi (baseMgNeutronMaterial_class.f90)'

    if (self % isFissile()) then
      ! Verify bounds
      if(G < 1 .or. self % nGroups() < G) then
        call fatalError(Here,'Invalid group number: '//numToChar(G)// &
                             ' Data has only: ' // numToChar(self % nGroups()))
        chi = ZERO ! Avoid warning
      end if
    
      chi = self % fission % data(G,2)
    else
      chi = ZERO
    end if

  end function getChi

  !!
  !! Return scatter XS for incoming energy group Gin and outgoing group Gout
  !!
  !! See mgNeutronMaterial documentationfor details
  !!
  function getScatterXS(self, Gin, Gout, rand) result(xs)
    class(baseMgNeutronMaterial), intent(in) :: self
    integer(shortInt), intent(in)            :: Gin
    integer(shortInt), intent(in)            :: Gout
    class(RNG), intent(inout)                :: rand
    real(defReal)                            :: xs
    character(100), parameter :: Here = ' getScatterXS (baseMgNeutronMaterial_class.f90)'

    ! Verify bounds
    if(Gin < 1 .or. self % nGroups() < Gin .or. Gout < 1 .or. self % nGroups() < Gout) then
      call fatalError(Here,'Invalid group numbers: '//numToChar(Gin)//' and '//numToChar(Gout) &
                           //' Data has only: ' // numToChar(self % nGroups()))
      xs = ZERO ! Avoid warning
    end if
    xs = self % scatter % P0(Gout,Gin) 
  
  end function getScatterXS


  !!
  !! Initialise Base MG Neutron Material fromdictionary
  !!
  !! Args:
  !!   dict       [in] -> Input dictionary with all required XSs
  !!   scatterKey [in] -> String with keyword to choose approperiate multiplicative scatering
  !!                        type
  !! Errors:
  !!   FatalError if scatteKey is invalid
  !!   FatalError if data in dictionary is invalid (inconsistant # of groups;
  !!     -ve entries in P0 XSs)
  !!
  !! Note:
  !!   Some time in the future scattering MG reaction objects will have factory. For now
  !!   the factory is hardcoded into this procedure. Not the best solution but is fine at this
  !!   stage. The following scatterKey are supported:
  !!     -> P0
  !!     -> P1
  !!
  subroutine init(self, dict, scatterKey)
    class(baseMgNeutronMaterial), intent(inout) :: self
    class(dictionary),target, intent(in)        :: dict
    character(nameLen), intent(in)              :: scatterKey
    integer(shortInt)                           :: nG, N, i
    real(defReal), dimension(:), allocatable    :: temp
    type(dictDeck)                              :: deck
    character(100), parameter :: Here = 'init (baseMgNeutronMaterial_class.f90)'


    ! Read number of groups
    call dict % get(nG, 'numberOfGroups')
    if(nG < 1) call fatalError(Here,'Number of groups is invalid' // numToChar(nG))
    self % nG = nG

    ! Set fissile flag
    call self % set(fissile = dict % isPresent('fission'))

    ! Build scattering reaction
    ! Prepare input deck
    deck % dict => dict

    ! Choose Scattering type
    select case(scatterKey)
      case ('P0')
        allocate( multiScatterMG :: self % scatter)

      case ('P1')
        allocate( multiScatterP1MG :: self % scatter)

      case default
        call fatalError(Here,'scatterKey: '//trim(scatterKey)//'is wrong. Must be P0 or P1')

    end select

    ! Initialise
    call self % scatter % init(deck, macroAllScatter)

    ! Deal with fission
    if(self % isFissile()) allocate(self % fission)
    if(self % isFissile()) call self % fission % init(deck, macroFission)

    ! Allocate space for data
    if(self % isFissile()) then
      N = 5
    else
      N = 3
    end if

    allocate(self % data(N, nG))

    ! Load cross sections
    call dict % get(temp, 'capture')
    if(size(temp) /= nG) then
      call fatalError(Here,'Capture XSs have wong size. Must be: ' &
                          // numToChar(nG)//' is '//numToChar(size(temp)))
    end if
    self % data(CAPTURE_XS,:) = temp

    ! Extract values of scattering XS
    if(size(self % scatter % scatterXSs) /= nG) then
      call fatalError(Here, 'Somthing went wrong. Inconsistant # of groups in material and reaction&
                            &. Clearly programming error.')
    end if
    self % data(IESCATTER_XS,:) = self % scatter % scatterXSs

    ! Load Fission-data
    if( self % isFissile()) then
      ! Load Fission
      call dict % get(temp, 'fission')
      if(size(temp) /= nG) then
        call fatalError(Here,'Fission XSs have wong size. Must be: ' &
                            // numToChar(nG)//' is '//numToChar(size(temp)))
      end if
      self % data(FISSION_XS,:) = temp

      ! Calculate nuFission
      call dict % get(temp, 'nu')
      if(size(temp) /= nG) then
        call fatalError(Here,'Nu vector has wong size. Must be: ' &
                            // numToChar(nG)//' is '//numToChar(size(temp)))
      end if
      self % data(NU_FISSION,:) = temp * self % data(FISSION_XS,:)
    end if

    ! Calculate total XS
    do i =1,nG
      self % data(TOTAL_XS, i) = self % data(IESCATTER_XS, i) + self % data(CAPTURE_XS, i)
      if(self % isFissile()) then
        self % data(TOTAL_XS, i) = self % data(TOTAL_XS, i) + self % data(FISSION_XS, i)
      end if
    end do
  end subroutine init

  !!
  !! Return number of energy groups
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  pure function nGroups(self) result(nG)
    class(baseMgNeutronMaterial), intent(in) :: self
    integer(shortInt)                        :: nG

    if(allocated(self % data)) then
      nG = self % nG
    else
      nG = 0
    end if

  end function nGroups
  
  !!
  !! Return pointer to Total XSs 
  !!
  function getTotalPtr(self) result(xs)
    class(baseMgNeutronMaterial), intent(in), target :: self
    real(defReal), dimension(:), pointer             :: xs

    xs => self % data(TOTAL_XS, :)

  end function getTotalPtr
  
  !!
  !! Return pointer to NuFission XSs 
  !!
  function getNuFissionPtr(self) result(xs)
    class(baseMgNeutronMaterial), intent(in), target :: self
    real(defReal), dimension(:), pointer             :: xs

    if (self % isFissile()) then
      xs => self % data(NU_FISSION, :)
    else
      xs => null()
    end if

  end function getNuFissionPtr
  
  !!
  !! Return pointer to Chis 
  !!
  function getChiPtr(self) result(chi)
    class(baseMgNeutronMaterial), intent(in), target :: self
    real(defReal), dimension(:), pointer             :: chi

    if (self % isFissile()) then
      chi => self % fission % data(:,2)
    else
      chi => null()
    end if

  end function getChiPtr
  
  !!
  !! Return pointer to scatter XSs 
  !!
  function getScatterPtr(self) result(xs)
    class(baseMgNeutronMaterial), intent(in), target :: self
    real(defReal), dimension(:,:), pointer           :: xs

    xs => self % scatter % P0(:, :)

  end function getScatterPtr
  
  !!
  !! Cast materialHandle pointer to baseMgNeutronMaterial type pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null if source is not of baseMgNeutronMaterial type
  !!   Target points to source if source is baseMgNeutronMaterialtype
  !!
  pure function baseMgNeutronMaterial_TptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    type(baseMgNeutronMaterial), pointer           :: ptr

    select type(source)
      type is(baseMgNeutronMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgNeutronMaterial_TptrCast

  !!
  !! Cast materialHandle pointer to baseMgNeutronMaterial class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null if source is not of baseMgNeutronMaterial class
  !!   Target points to source if source is baseMgNeutronMaterial class
  !!
  pure function baseMgNeutronMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(baseMgNeutronMaterial), pointer          :: ptr

    select type(source)
      class is(baseMgNeutronMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function baseMgNeutronMaterial_CptrCast


end module baseMgNeutronMaterial_class
