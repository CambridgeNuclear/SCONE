module ceNeutronMaterial_class

  use numPrecision
  use genericProcedures, only : fatalError, binarySearch, numToChar
  use RNG_class,         only : RNG
  use particle_class,    only : particle

  ! Nuclear Data Handles
  use materialHandle_inter,    only : materialHandle
  use neutronMaterial_inter,   only : neutronMaterial
  use neutronXsPackages_class, only : neutronMacroXSs

  ! CE Neutron Interfaces
  use ceNeutronDatabase_inter, only : ceNeutronDatabase

  ! Cache
  use ceNeutronCache_mod,      only : materialCache, nuclideCache

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public ceNeutronMaterial_CptrCast
  public ceNeutronMaterial_TptrCast

  !!
  !! An abstract class that represent all CE Neutron Material Data
  !!
  !! Exist mainly in order to decouple caching logic from the database implementation
  !! so there is no need to repeat it in every database type. Thus it will be easier to
  !! mantain and optimise.
  !!
  !! Note that a material without any composition is not allowed.
  !! Makes no assumption about the range of nucIdx. Allows for -ve values
  !!
  !! Interface:
  !!   materialHandle Interface
  !!   getMacroXSs -> return package of macroscopic XSs directly from Energy and RNG
  !!   set         -> Set data related to material by keyword association
  !!   setComposition -> Set composition of material from densities and nucIdxs
  !!   sampleNuclide  -> sample collision nuclide
  !!   sampleFission  -> sample collision nuclide given that fission reaction has happened
  !!   sampleScatter  -> sample collision nuclide given that Scattering has happened
  !!   sampleScatterWithFission -> sample collision nuclide given ther scatter or fission has
  !!     happened
  !!
  type, public, extends(neutronMaterial) :: ceNeutronMaterial
    integer(shortInt)                            :: matIdx = 0
    class(ceNeutronDatabase), pointer            :: data => null()
    real(defReal), dimension(:), allocatable     :: dens
    integer(shortInt), dimension(:), allocatable :: nuclides
    logical(defBool)                             :: fissile =.false.
    real(defReal), dimension(:), allocatable     :: unionGrid
    real(defReal), dimension(:), allocatable     :: totalXS
    integer(shortInt),dimension(:,:),allocatable :: gridIdx

  contains
    ! Superclass procedures
    procedure :: kill
    generic   :: getMacroXSs => getMacroXSs_byE
    procedure :: getMacroXSs_byP

    ! Local procedures
    procedure, non_overridable :: set
    procedure, non_overridable :: setComposition
    procedure, non_overridable :: getMacroXSs_byE
    procedure                  :: isFissile
    procedure, non_overridable :: sampleNuclide
    procedure, non_overridable :: sampleFission
    procedure, non_overridable :: sampleScatter
    procedure, non_overridable :: sampleScatterWithFission
    procedure                  :: search

  end type ceNeutronMaterial

contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(ceNeutronMaterial), intent(inout) :: self

    self % matIdx  = 0
    self % data    => null()
    if(allocated(self % dens))     deallocate(self % dens )
    if(allocated(self % nuclides)) deallocate (self % nuclides)
    self % fissile = .false.

  end subroutine kill

  !!
  !! Return Macroscopic XSs for the material given particle
  !!
  !! See neutronMaterial_inter for details
  !!
  subroutine getMacroXSs_byP(self, xss, p)
    class(ceNeutronMaterial), intent(in) :: self
    type(neutronMacroXSs), intent(out)   :: xss
    class(particle), intent(in)          :: p
    character(100), parameter :: Here = 'getMacroXSs_byP (ceNeutronMaterial_class.f90)'

    if(.not.p % isMG) then
      call self % getMacroXSs(xss, p % E, p % pRNG)

    else
      call fatalError(Here,'MG neutron given to CE data')

    end if
  end subroutine getMacroXSs_byP

  !!
  !! Set composition of the material in terms of nucIdx and atomic density
  !!
  !! Use this procedure ONLY during build. NEVER during transport.
  !! IT IS NOT THREAD SAFE!
  !!
  !! Args:
  !!   dens    [in] -> array of atomic densities [1/barn/cm] of nuclides
  !!   nucIdxs [in] -> correpsonding array with nucIdxs
  !!
  !! Errors:
  !!   FatalError if arrays have diffrent size
  !!   FatalError if dens contains -ve values
  !!   FatalError if dens has size of 0 -> no composition
  !!
  subroutine setComposition(self, dens, nucIdxs)
    class(ceNeutronMaterial), intent(inout)     :: self
    real(defReal), dimension(:), intent(in)     :: dens
    integer(shortInt), dimension(:), intent(in) :: nucIdxs
    character(100), parameter :: Here = 'setComposition (ceNeutronMaterial_class.f90)'

    ! Check input
    if(size(dens) /= size(nucIdxs)) call fatalError(Here,'Diffrent sizes of density and nuclide vector')
    if(any(dens < ZERO)) call fatalError(Here,'-ve nuclide densities are present')
    if(size(dens) == 0)  call fatalError(Here,'Empty composition is not allowed')

    ! Clean any current content
    if(allocated(self % dens))     deallocate(self % dens)
    if(allocated(self % nuclides)) deallocate(self % nuclides)

    ! Load values
    self % dens     = dens
    self % nuclides = nucIdxs

  end subroutine setComposition

  !!
  !! Set matIdx, pointer to a database and fissile flag
  !!
  !! All arguments are optional. Use with keyword association e.g.
  !!   call mat % set( matIdx = 7)
  !!
  !! Use this procedure ONLY during build. NEVER during transport.
  !! IT IS NOT THREAD SAFE!
  !!
  !! Args:
  !!   matIdx [in]    -> material index
  !!   database [in]  -> pointer to a database that updates XSs on the ceNeutronCache
  !!   fissile [in] -> flag indicating whether fission data is present
  !!
  subroutine set(self, matIdx, database, fissile)
    class(ceNeutronMaterial), intent(inout)               :: self
    integer(shortInt), intent(in),optional                :: matIdx
    class(ceNeutronDatabase),pointer, optional,intent(in) :: database
    logical(defBool),intent(in), optional                 :: fissile

    if(present(matIdx))    self % matIdx  = matIdx
    if(present(database))  self % data    => database
    if(present(fissile))   self % fissile = fissile

  end subroutine set

  !!
  !! Return Macroscopic XSs for the material
  !!
  !! Args:
  !!   xss [out]    -> Cross section package to store the data
  !!   E [in]       -> Requested energy [MeV]
  !!   rand [inout] -> Random Number Generator
  !!
  !! Errors:
  !!   fatalError if E is out-of-bounds for the stored data
  !!
  subroutine getMacroXSs_byE(self, xss, E, rand)
    class(ceNeutronMaterial), intent(in) :: self
    type(neutronMacroXSs), intent(out)   :: xss
    real(defReal), intent(in)            :: E
    class(RNG), intent(inout)            :: rand

    ! Check Cache and update if needed
    if(materialCache(self % matIdx) % E_tail /= E) then
      call self % data % updateMacroXSs(E, self % matIdx, rand)
    end if

    xss = materialCache(self % matIdx) % xss

  end subroutine getMacroXSs_byE

  !!
  !! Return .true. if nuclide is fissile
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   .TRUE. if fissile, .FALSE. otherwise
  !!
  !! Errors:
  !!   None
  !!
  elemental function isFissile(self) result(isIt)
    class(ceNeutronMaterial), intent(in) :: self
    logical(defBool)                     :: isIt

    isIt = self % fissile

  end function isFissile

  !!
  !! Sample collision nuclide at energy E
  !!
  !! This function randomly determines the exact nuclide for a collision
  !! It uses nuclide total XSs to determine nuclide
  !!
  !! Args:
  !!   E [in]       -> incident energy [MeV]
  !!   rand [inout] -> random number generator
  !!
  !! Result:
  !!   nucIdx of the sampled nuclide for collision
  !!
  !! Errors:
  !!   fatalError if sampling fails for some reason (E.G. random number > 1)
  !!   fatalError if E is out-of-bounds of the present data
  !!
  function sampleNuclide(self, E, rand) result(nucIdx)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: E
    class(RNG), intent(inout)            :: rand
    integer(shortInt)                    :: nucIdx
    real(defReal)                        :: xs
    integer(shortInt)                    :: i
    character(100), parameter :: Here = 'sampleNuclide (ceNeutronMaterial_class.f90)'

    ! Get total material XS
    if(E /= materialCache(self % matIdx) % E_tot) then
      call self % data % updateTotalMatXS(E, self % matIdx, rand)
    end if

    xs = materialCache(self % matIdx) % xss % total * rand % get()

    ! Loop over all nuclides
    do i=1,size(self % nuclides)
      nucIdx = self % nuclides(i)
      if(E /= nuclideCache(nucIdx) % E_tot) call self % data % updateTotalNucXS(E, nucIdx, rand)
      xs = xs - nuclideCache(nucIdx) % xss % total * self % dens(i)
      if(xs < ZERO) return
    end do

    ! Print error message as the inversion failed
    call fatalError(Here,'Nuclide sampling loop failed to terminate')

  end function sampleNuclide

  !!
  !! Sample fission nuclide given that a fission neutron was produced
  !!
  !! Basicly samples from P(nucIdx| fission neutron produced in material)
  !! Usefull when generating fission sites
  !!
  !! As such it uses nu*sigma_f
  !!
  !! For a non-fissile material return nucIdx <= 0 !
  !!
  !! Args:
  !!   E [in]       -> incident energy [MeV]
  !!   rand [inout] -> random number generator
  !!
  !! Result:
  !!   nucIdx of the sampled nuclide for collision.
  !!
  !! Errors:
  !!   fatalError if sampling fails for some reason (E.G. random number > 1)
  !!   fatalError if E is out-of-bounds of the present data
  !!   Returns nucIdx <= if material is not fissile
  !!
  function sampleFission(self, E, rand) result(nucIdx)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: E
    class(RNG), intent(inout)            :: rand
    integer(shortInt)                    :: nucIdx
    real(defReal)                        :: xs
    integer(shortInt)                    :: i
    character(100), parameter :: Here = 'sampleFission (ceNeutronMaterial_class.f90)'

    ! Short-cut for nonFissile material
    if(.not.self % fissile) then
      nucIdx = 0
      return
    end if

    ! Calculate material macroscopic nuFission
    if(E /= materialCache(self % matIdx) % E_tail) then
      call self % data % updateMacroXSs(E, self % matIdx, rand)
    end if

    xs = materialCache(self % matIdx) % xss % nuFission * rand % get()

    ! Loop over all nuclides
    do i=1,size(self % nuclides)
      nucIdx = self % nuclides(i)
      if(E /= nuclideCache(nucIdx) % E_tail) call self % data % updateMicroXSs(E, nucIdx, rand)
      xs = xs - nuclideCache(nucIdx) % xss % nuFission * self % dens(i)
      if(xs < ZERO) return
    end do

    ! Print error message as the inversion failed
    call fatalError(Here,'Nuclide sampling loop failed to terminate')

  end function sampleFission

  !!
  !! Sample collision nuclide given that any scattering has happened
  !! Treats fission as a capture!
  !!
  !! For a pure-absorbing material return nucIdx <= 0 !
  !!
  !! Args:
  !!   E [in]       -> incident energy [MeV]
  !!   rand [inout] -> random number generator
  !!
  !! Result:
  !!   nucIdx of the sampled nuclide for collision.
  !!
  !! Errors:
  !!   fatalError if sampling fails for some reason (E.G. random number > 1)
  !!   fatalError if E is out-of-bounds of the present data
  !!   Returns nucIdx <= if material is a pure-absorber (with fission as absorbtion)
  !!
  function sampleScatter(self, E, rand) result(nucIdx)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: E
    class(RNG), intent(inout)            :: rand
    integer(shortInt)                    :: nucIdx
    real(defReal)                        :: xs
    integer(shortInt)                    :: i
    character(100), parameter :: Here = 'sampleScatter (ceNeutronMaterial_class.f90)'

    ! Calculate material macroscopic cross section of all scattering
    if(E /= materialCache(self % matIdx) % E_tail) then
      call self % data % updateMacroXSs(E, self % matIdx, rand)
    end if

    xs = rand % get() * (materialCache(self % matIdx) % xss % elasticScatter + &
                         materialCache(self % matIdx) % xss % inelasticScatter)

    ! Loop over all nuclides
    do i=1,size(self % nuclides)
      nucIdx = self % nuclides(i)
      if(E /= nuclideCache(nucIdx) % E_tail) call self % data % updateMicroXSs(E, nucIdx, rand)
      xs = xs - (nuclideCache(nucIdx) % xss % elasticScatter + &
                 nuclideCache(nucIdx) % xss % inelasticScatter ) * self % dens(i)
      if(xs < ZERO) return
    end do

    ! Print error message as the inversion failed
    call fatalError(Here,'Nuclide sampling loop failed to terminate')

  end function sampleScatter

  !!
  !! Sample collision nuclide given that any scattering or fission has happened
  !! Treats fission as a scattering!
  !!
  !! For a pure-capture material return nucIdx <= 0 !
  !!
  !! Args:
  !!   E [in]       -> incident energy [MeV]
  !!   rand [inout] -> random number generator
  !!
  !! Result:
  !!   nucIdx of the sampled nuclide for collision.
  !!
  !! Errors:
  !!   fatalError if sampling fails for some reason (E.G. random number > 1)
  !!   fatalError if E is out-of-bounds of the present data
  !!   Returns nucIdx <= if material is a pure-capture (with fission as scattering)
  !!
  function sampleScatterWithFission(self, E, rand) result(nucIdx)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: E
    class(RNG), intent(inout)            :: rand
    integer(shortInt)                    :: nucIdx
    real(defReal)                        :: xs
    integer(shortInt)                    :: i
    character(100), parameter :: Here = 'sampleScatterWithFission (ceNeutronMaterial_class.f90)'

    ! Calculate material macroscopic cross section of all scattering
    if(E /= materialCache(self % matIdx) % E_tail) then
      call self % data % updateMacroXSs(E, self % matIdx, rand)
    end if

    xs = rand % get() * (materialCache(self % matIdx) % xss % elasticScatter + &
                         materialCache(self % matIdx) % xss % inelasticScatter + &
                         materialCache(self % matIdx) % xss % fission)

    ! Loop over all nuclides
    do i=1,size(self % nuclides)
      nucIdx = self % nuclides(i)
      if(E /= nuclideCache(nucIdx) % E_tail) call self % data % updateMicroXSs(E, nucIdx, rand)
      xs = xs - (nuclideCache(nucIdx) % xss % elasticScatter + &
                 nuclideCache(nucIdx) % xss % inelasticScatter + &
                 nuclideCache(nucIdx) % xss % fission ) * self % dens(i)
      if(xs < ZERO) return
    end do

    ! Print error message as the inversion failed
    call fatalError(Here,'Nuclide sampling loop failed to terminate')

  end function sampleScatterWithFission

  !!
  !! Cast materialHandle pointer to ceNeutronMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of ceNeutronMaterial
  !!   Pointer to source if source is ceNeutronMaterial class
  !!
  pure function ceNeutronMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(ceNeutronMaterial), pointer          :: ptr

    select type(source)
      class is(ceNeutronMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function ceNeutronMaterial_CptrCast

  !!
  !! Cast materialHandle pointer to ceNeutronMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of ceNeutronMaterial
  !!   Pointer to source if source is ceNeutronMaterial class
  !!
  pure function ceNeutronMaterial_TptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    type(ceNeutronMaterial), pointer            :: ptr

    select type(source)
      type is(ceNeutronMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function ceNeutronMaterial_TptrCast

  !!
  !! Search energy for grid and interpolation factor for energy E
  !! Called only by aceNeutronDatabaseUni and UniIdx
  !!
  !! Interpolation factor definition:
  !!   f = (E - E_low) / (E_top - E_low)
  !!   E = E_top * f + E_low * (1-f)
  !!
  !! Args:
  !!   idx [out] -> index of the bottom bin for energy E
  !!   f   [out] -> value of the interpolation factor for energy E
  !!   E   [in]  -> Energy to search for [MeV]
  !!
  !! Errors:
  !!   If energy E is beyond range terminate with fatalError
  !!
  subroutine search(self, idx, f, E)
    class(ceNeutronMaterial), intent(in) :: self
    integer(shortInt), intent(out)       :: idx
    real(defReal), intent(out)           :: f
    real(defReal), intent(in)            :: E
    character(100), parameter :: Here = 'search (ceNeutronMaterialUni_class.f90)'

    idx = binarySearch(self % unionGrid, E)
    if(idx <= 0) then
      call fatalError(Here,'Failed to find energy: '//numToChar(E)//&
                           ' for material '//numToChar(self % matIdx))
    end if

    associate(E_top => self % unionGrid(idx + 1), E_low  => self % unionGrid(idx))
      f = (E - E_low) / (E_top - E_low)
    end associate

  end subroutine search


end module ceNeutronMaterial_class
