module ceNeutronMaterial_class

  use numPrecision
  use universalVariables
  use genericProcedures,  only : fatalError, numToChar
  use RNG_class,          only : RNG
  use particle_class,     only : particle

  ! Nuclear Data Handles
  use materialHandle_inter,    only : materialHandle
  use neutronMaterial_inter,   only : neutronMaterial
  use neutronXsPackages_class, only : neutronMacroXSs

  ! CE Neutron Interfaces
  use ceNeutronDatabase_inter, only : ceNeutronDatabase
  use ceNeutronNuclide_inter,  only : ceNeutronNuclide, ceNeutronNuclide_CptrCast

  ! Cache
  use ceNeutronCache_mod,      only : materialCache, nuclideCache

  ! Scattering procedures
  use scatteringKernels_func,  only : relativeEnergy_constXS, dopplerCorrectionFactor

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
    real(defReal)                                :: kT     = ZERO
    class(ceNeutronDatabase), pointer            :: data => null()
    real(defReal), dimension(:), allocatable     :: dens
    integer(shortInt), dimension(:), allocatable :: nuclides
    logical(defBool)                             :: fissile = .false.
    logical(defBool)                             :: hasTMS  = .false.
    real(defReal)                                :: eUpperSab = ZERO
    real(defReal)                                :: eLowerURR = ZERO

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
    procedure                  :: useTMS
    procedure, non_overridable :: sampleNuclide
    procedure, non_overridable :: sampleFission
    procedure, non_overridable :: sampleScatter
    procedure, non_overridable :: sampleScatterWithFission

  end type ceNeutronMaterial

contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(ceNeutronMaterial), intent(inout) :: self

    self % matIdx  = 0
    self % kT      = ZERO
    self % data    => null()
    if (allocated(self % dens))     deallocate(self % dens)
    if (allocated(self % nuclides)) deallocate (self % nuclides)
    self % fissile = .false.
    self % hasTMS  = .false.
    self % eUpperSab = ZERO
    self % eLowerURR = ZERO

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

    if (.not.p % isMG) then
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
  !!   FatalError if arrays have different size
  !!   FatalError if dens contains -ve values
  !!   FatalError if dens has size of 0 -> no composition
  !!
  subroutine setComposition(self, dens, nucIdxs)
    class(ceNeutronMaterial), intent(inout)     :: self
    real(defReal), dimension(:), intent(in)     :: dens
    integer(shortInt), dimension(:), intent(in) :: nucIdxs
    character(100), parameter :: Here = 'setComposition (ceNeutronMaterial_class.f90)'

    ! Check input
    if (size(dens) /= size(nucIdxs)) call fatalError(Here,'Different sizes of density and nuclide vector')
    if (any(dens < ZERO)) call fatalError(Here,'-ve nuclide densities are present')
    if (size(dens) == 0)  call fatalError(Here,'Empty composition is not allowed')

    ! Clean any current content
    if (allocated(self % dens))     deallocate(self % dens)
    if (allocated(self % nuclides)) deallocate(self % nuclides)

    ! Load values
    self % dens     = dens
    self % nuclides = nucIdxs

  end subroutine setComposition

  !!
  !! Set matIdx, pointer to a database and fissile flag
  !!
  !! All arguments are optional. Use with keyword association e.g.
  !!   call mat % set(matIdx = 7)
  !!
  !! Use this procedure ONLY during build. NEVER during transport.
  !! IT IS NOT THREAD SAFE!
  !!
  !! NOTE: eUpperSab and eLowerURR are fed by the aceNeutronDatabase, and they are the
  !! strictest (respecitvely highest and lowest) energy limits among all nuclides
  !! in the material composition.
  !!
  !! Args:
  !!   matIdx [in]    -> material index
  !!   database [in]  -> pointer to a database that updates XSs on the ceNeutronCache
  !!   fissile [in]   -> flag indicating whether fission data is present
  !!   hasTMS [in]    -> flag indicating whether TMS is on
  !!   temp [in]      -> TMS material temperature
  !!   eUpperSab [in] -> upper energy of S(a,b) range in the material
  !!   eLowerURR [in] -> lower energy of ures range in the material
  !!
  subroutine set(self, matIdx, database, fissile, hasTMS, temp, eUpperSab, eLowerURR)
    class(ceNeutronMaterial), intent(inout)                 :: self
    integer(shortInt), intent(in), optional                 :: matIdx
    class(ceNeutronDatabase), pointer, optional, intent(in) :: database
    logical(defBool), intent(in), optional                  :: fissile
    logical(defBool), intent(in), optional                  :: hasTMS
    real(defReal), intent(in), optional                     :: temp
    real(defReal), intent(in), optional                     :: eUpperSab
    real(defReal), intent(in), optional                     :: eLowerURR
    character(100), parameter :: Here = 'set (ceNeutronMaterial_class.f90)'

    if (present(database))  self % data    => database
    if (present(fissile))   self % fissile = fissile
    if (present(matIdx))    self % matIdx  = matIdx
    if (present(hasTMS))    self % hasTMS  = hasTMS
    if (present(temp))      self % kT = (kBoltzmann * temp) / joulesPerMeV

    if (present(eUpperSab)) then
      if (eUpperSab < ZERO) call fatalError (Here, 'Upper Sab energy limit of material '&
                                             &//numToChar(matIdx)//' is negative.')
      self % eUpperSab  = eUpperSab
    end if

    if (present(eLowerURR)) then
      if (eLowerURR < ZERO) call fatalError (Here, 'Lower URR energy limit of material '&
                                             &//numToChar(matIdx)//' is negative.')
      self % eLowerURR  = eLowerURR
    end if

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
    if (materialCache(self % matIdx) % E_tail /= E .or. materialCache(self % matIdx) % E_tot /= E) then
      call self % data % updateMacroXSs(E, self % matIdx, rand)
    end if

    xss = materialCache(self % matIdx) % xss

  end subroutine getMacroXSs_byE

  !!
  !! Return .true. if material is fissile
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   .true. if fissile, .false. otherwise
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
  !! Return .true. if TMS is on in the material and the provided energy is not
  !! within the ures or S(a,b) range of any nuclide in the material.
  !!
  !! Args:
  !!   E [in] -> test energy
  !!
  !! Result:
  !!   .true. if conditions for using TMS are satisfied, .false. otherwise
  !!
  !! Errors:
  !!   None
  !!
  elemental function useTMS(self, E) result(shouldIt)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: E
    logical(defBool)                     :: shouldIt

    shouldIt = self % hasTMS .and. (E > self % eUpperSab) .and. (E < self % eLowerURR)

  end function useTMS

  !!
  !! Sample collision nuclide at energy E
  !!
  !! This function randomly determines the exact nuclide for a collision
  !! It uses nuclide total XSs to determine nuclide
  !!
  !! Args:
  !!   E [in]       -> incident energy [MeV]
  !!   rand [inout] -> random number generator
  !!   nucIdx [out] -> sampled nuclide index
  !!   eOut [out]   -> relative energy between neutron and target (may be /= E in case of TMS)
  !!
  !! Errors:
  !!   fatalError if sampling fails for some reason (E.G. random number > 1)
  !!   fatalError if E is out-of-bounds of the present data
  !!
  subroutine sampleNuclide(self, E, rand, nucIdx, eOut)
    class(ceNeutronMaterial), intent(in) :: self
    real(defReal), intent(in)            :: E
    class(RNG), intent(inout)            :: rand
    integer(shortInt), intent(out)       :: nucIdx
    real(defReal), intent(out)           :: eOut
    class(ceNeutronNuclide), pointer     :: nuc
    integer(shortInt)                    :: i
    real(defReal)                        :: P_acc, eMin, eMax, A, eRel, &
                                            trackMatXS, totNucXS, dens
    character(100), parameter :: Here = 'sampleNuclide (ceNeutronMaterial_class.f90)'

    ! Get material tracking XS
    if (E /= materialCache(self % matIdx) % E_track) then
      call self % data % updateTrackMatXS(E, self % matIdx, rand)
    end if

    trackMatXS = materialCache(self % matIdx) % trackXS * rand % get()

    ! Loop over nuclides
    do i = 1,size(self % nuclides)

      nucIdx = self % nuclides(i)
      dens = self % dens(i)

      associate (nucCache => nuclideCache(nucIdx))

        ! Retrieve nuclide XS from cache
        if (self % useTMS(E)) then

          ! If the material is using TMS, the nuclide temperature majorant is needed
          ! The check for the right values stored in cache happens inside the subroutine
          call self % data % updateTotalTempNucXS(E, self % kT, nucIdx)

          totNucXS = nucCache % tempMajXS * nucCache % doppCorr

        else

          ! Update nuclide cache if needed
          if (E /= nucCache % E_tot) call self % data % updateTotalNucXS(E, nucIdx, rand)
          totNucXS = nucCache % xss % total

        end if

        trackMatXS = trackMatXS - totNucXS * dens

        ! Nuclide temporarily accepted: check TMS condition
        if (trackMatXS < ZERO) then

          ! Save energy to be used to sample reaction
          eOut = E

          if (self % useTMS(E)) then

            ! If the material is using TMS, retrieve nuclide and nuclide information
            nuc => ceNeutronNuclide_CptrCast(self % data % getNuclide(nucIdx))
            if (.not. associated(nuc)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

            A = nuc % getMass()

            ! Sample relative energy
            eRel = relativeEnergy_constXS(E, A, nucCache % deltakT, rand)

            ! Call through system minimum and maximum energies
            call self % data % energyBounds(eMin, eMax)

            ! Ensure relative energy is within energy bounds
            if (eRel < eMin) eRel = eMin
            if (eMax < eRel) eRel = eMax

            ! Get relative energy nuclide cross section
            call self % data % updateTotalNucXS(eRel, nucIdx, rand)

            ! Calculate acceptance probability using ratio of relative energy xs to temperature majorant
            P_acc = nucCache % xss % total * nucCache % doppCorr / totNucXS

            ! Accept or reject the sampled nuclide
            if (rand % get() >= P_acc) nucIdx = REJECTED

            ! Overwrite energy to be used to sample reaction
            eOut = eRel

          end if

          ! Exit function, return the sampled nucIdx
          return

        end if

      end associate

    end do

    ! Print error message as the inversion failed
    call fatalError(Here,'Nuclide sampling loop failed to terminate')

  end subroutine sampleNuclide

  !!
  !! Sample fission nuclide given that a fission neutron was produced
  !!
  !! Basically samples from P(nucIdx| fission neutron produced in material)
  !! Useful when generating fission sites
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
    class(ceNeutronNuclide), pointer     :: nuc
    integer(shortInt)                    :: nucIdx, i
    real(defReal)                        :: xs, doppCorr, A, nuckT, deltakT
    character(100), parameter :: Here = 'sampleFission (ceNeutronMaterial_class.f90)'

    ! Short-cut for nonFissile material
    if (.not. self % fissile) then
      nucIdx = 0
      return
    end if

    ! Calculate material macroscopic nuFission
    ! The cache is updated without checking the energy to get the correct results with TMS
    ! The relative energy flag cached is cleaned to make sure cross sections are updated
    materialCache(self % matIdx) % E_rel = ZERO
    call self % data % updateMacroXSs(E, self % matIdx, rand)

    xs = materialCache(self % matIdx) % xss % nuFission * rand % get()

    ! Loop over all nuclides
    do i = 1,size(self % nuclides)

      nucIdx = self % nuclides(i)

      ! If the material uses TMS, the Doppler correction factor is needed
      if (self % useTMS(E)) then

        nuc => ceNeutronNuclide_CptrCast(self % data % getNuclide(nucIdx))
        if (.not. associated(nuc)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

        A     = nuc % getMass()
        nuckT = nuc % getkT()
        deltakT = self % kT - nuckT
        doppCorr = dopplerCorrectionFactor(E, A, deltakT)

      else
        doppCorr = ONE
      end if

      ! The nuclide cache should be at the right energy after updating the material
      ! In the case of TMS where the macro xss are at a relative energy, the nuclide
      ! xss to be used are at the relative energy just sampled
      xs = xs - nuclideCache(nucIdx) % xss % nuFission * self % dens(i) * doppCorr

      if (xs < ZERO) return

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
    class(ceNeutronNuclide), pointer     :: nuc
    integer(shortInt)                    :: nucIdx, i
    real(defReal)                        :: xs, doppCorr, A, nuckT, deltakT
    character(100), parameter :: Here = 'sampleScatter (ceNeutronMaterial_class.f90)'

    ! Calculate material macroscopic cross section of all scattering
    ! The cache is updated without checking the energy to get the correct results with TMS
    ! The relative energy flag cached is cleaned to make sure cross sections are updated
    materialCache(self % matIdx) % E_rel = ZERO
    call self % data % updateMacroXSs(E, self % matIdx, rand)

    xs = rand % get() * (materialCache(self % matIdx) % xss % elasticScatter + &
                         materialCache(self % matIdx) % xss % inelasticScatter)

    ! Loop over all nuclides
    do i = 1,size(self % nuclides)

      nucIdx = self % nuclides(i)

      ! If the material uses TMS, the Doppler correction factor is needed
      if (self % useTMS(E)) then

        ! If the material is using TMS, the Doppler correction factor is needed
        nuc => ceNeutronNuclide_CptrCast(self % data % getNuclide(nucIdx))
        if (.not. associated(nuc)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

        A     = nuc % getMass()
        nuckT = nuc % getkT()
        deltakT = self % kT - nuckT
        doppCorr = dopplerCorrectionFactor(E, A, deltakT)

      else
        doppCorr = ONE
      end if

      ! The nuclide cache should be at the right energy after updating the material
      ! In the case of TMS where the macro xss are at a relative energy, the nuclide
      ! xss to be used are at the relative energy just sampled
      xs = xs - (nuclideCache(nucIdx) % xss % elasticScatter + &
                 nuclideCache(nucIdx) % xss % inelasticScatter ) * self % dens(i) * doppCorr

      if (xs < ZERO) return

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
    class(ceNeutronNuclide), pointer     :: nuc
    integer(shortInt)                    :: nucIdx, i
    real(defReal)                        :: xs, doppCorr, A, nuckT, deltakT
    character(100), parameter :: Here = 'sampleScatterWithFission (ceNeutronMaterial_class.f90)'

    ! Calculate material macroscopic cross section of all scattering and fission
    ! The cache is updated without checking the energy to get the correct results with TMS
    ! The relative energy flag cached is cleaned to make sure cross sections are updated
    materialCache(self % matIdx) % E_rel = ZERO
    call self % data % updateMacroXSs(E, self % matIdx, rand)

    xs = rand % get() * (materialCache(self % matIdx) % xss % elasticScatter + &
                         materialCache(self % matIdx) % xss % inelasticScatter + &
                         materialCache(self % matIdx) % xss % fission)

    ! Loop over all nuclides
    do i = 1,size(self % nuclides)

      nucIdx = self % nuclides(i)

      ! If the material uses TMS, the Doppler correction factor is needed
      if (self % useTMS(E)) then

        ! If the material is using TMS, the Doppler correction factor is needed
        nuc => ceNeutronNuclide_CptrCast(self % data % getNuclide(nucIdx))
        if (.not. associated(nuc)) call fatalError(Here, 'Failed to retrieve CE Neutron Nuclide')

        A     = nuc % getMass()
        nuckT = nuc % getkT()
        deltakT = self % kT - nuckT
        doppCorr = dopplerCorrectionFactor(E, A, deltakT)

      else
        doppCorr = ONE
      end if

      ! The nuclide cache should be at the right energy after updating the material
      ! In the case of TMS where the macro xss are at a relative energy, the nuclide
      ! xss to be used are at the relative energy just sampled
      xs = xs - (nuclideCache(nucIdx) % xss % elasticScatter + &
                 nuclideCache(nucIdx) % xss % inelasticScatter + &
                 nuclideCache(nucIdx) % xss % fission) * self % dens(i) * doppCorr

      if (xs < ZERO) return

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


end module ceNeutronMaterial_class
