module mgXsClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use universalVariables
  use genericProcedures,          only : fatalError
  use dictionary_class,           only : dictionary
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use outputFile_class,           only : outputFile

  ! Nuclear Data interface
  use nuclearDatabase_inter,      only : nuclearDatabase
  use neutronXSPackages_class,    only : neutronMacroXSs
  use neutronMaterial_inter,      only : neutronMaterial,neutronMaterial_CptrCast

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Tally Interfaces
  use scoreMemory_class,          only : scoreMemory
  use tallyClerk_inter,           only : tallyClerk, kill_super => kill

  implicit none
  private

  !! Size of clerk memory
  integer(shortInt), parameter  :: ARRAY_SCORE_SIZE   = 7 ,&  ! Size of data to store as 1D arrays
                                   MATRIX_SCORE_SMALL = 3 ,&  ! Size of data to store as 2D arrays when scoring up to P1
                                   MATRIX_SCORE_FULL  = 9     ! Size of data to store as 2D arrays when scoring up to P7

  !! Locations of different bins wrt memory address of the clerk
  integer(longInt), parameter   :: FLUX_idx      = 1 ,&  ! Flux
                                   SCATT_idx     = 2 ,&  ! Scattering macroscopic reaction rate
                                   CAPT_idx      = 3 ,&  ! Capture macroscopic reaction rate
                                   FISS_idx      = 4 ,&  ! Fission macroscopic reaction rate
                                   NUBAR_idx     = 5 ,&  ! NuBar
                                   CHI_idx       = 6 ,&  ! Fission neutron spectrum
                                   SCATT_EV_idx  = 7     ! Analog: number of scattering events

  !!
  !! Multi-group macroscopic cross section calculation
  !!
  !! It prints out:
  !! capture xs, fission xs, transport xs, nu, chi, the P0 and P1 scattering matrices,
  !! the P0 scattering production matrix. On request, also the P2 -> P7 scattering matrices.
  !!
  !! NOTE:
  !! - the cross sections are tallied with a collision estimator;
  !! - the scattering matrices and chi are calculated with an analog estimator: their
  !!   probability is calculated directly by scoring events
  !! - the transport cross section is calculated with both the out-scatter and
  !!   flux-limited approximation
  !!
  !! Private Members:
  !!   energyMap -> tally map for energy group structure
  !!   spacemap  -> tally map for material or spatial bins
  !!   energyN   -> number of energy groups
  !!   matN      -> number of materials or spatial bins
  !!   width     -> number of memory bins needed
  !!   PN        -> flag to indicate whether matrices above P1 are needed
  !!
  !! Interface
  !!   tallyClerk Interface
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myMgXsClerk {
  !!   type mgXsClerk;
  !!   # energyMap { energyMap definition } #
  !!   # spaceMap  { <other tallyMap definition> } #
  !!   # PN 1; #
  !! }
  !!
  !!  NOTE: spaceMap can be any type of spatial tallyMap, including a multiMap
  !!  NOTE: PN is a flag that indicates whether scattering matrices from P2 to P7
  !!        have to be calculated
  !!
  type, public, extends(tallyClerk) :: mgXsClerk
    private
    ! Maps
    class(tallyMap), allocatable :: spaceMap
    class(tallyMap), allocatable :: energyMap

    ! Useful data
    integer(shortInt) :: energyN = 0
    integer(shortInt) :: matN = 0
    integer(shortInt) :: width = 0
    logical(defBool)  :: PN = .false.

    ! Settings
    logical(defBool) :: handleVirtual = .true.

  contains
    ! Procedures used during build
    procedure  :: init
    procedure  :: kill
    procedure  :: validReports
    procedure  :: getSize

    ! File reports -> run-time procedures
    procedure  :: reportInColl
    procedure  :: reportOutColl
    procedure  :: reportSpawn

    ! Output procedures
    procedure  :: print
    procedure  :: processRes
    procedure  :: processPN
    procedure  :: display

  end type mgXsClerk


contains

  !!
  !! Initialise clerk from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(mgXsClerk), intent(inout)   :: self
    class(dictionary), intent(in)     :: dict
    character(nameLen), intent(in)    :: name

    ! Assign name
    call self % setName(name)

    ! Load energy map and bin number
    if (dict % isPresent('energyMap')) then
      call new_tallyMap(self % energyMap, dict % getDictPtr('energyMap'))
      self % energyN = self % energyMap % bins(0)
    else
      self % energyN = 1
    end if

    ! Load space/material map and bin number
    if (dict % isPresent('spaceMap')) then
      call new_tallyMap(self % spaceMap, dict % getDictPtr('spaceMap'))
      self % matN = self % spaceMap % bins(0)
    else
      self % matN = 1
    end if

    ! Get PN flag
    call dict % getOrDefault(self % PN, 'PN', .true.)

    ! Set memory width
    if (self % PN) then
      self % width = ARRAY_SCORE_SIZE + MATRIX_SCORE_FULL * self % energyN
    else
      self % width = ARRAY_SCORE_SIZE + MATRIX_SCORE_SMALL * self % energyN
    end if

    ! Handle virtual collisions
    call dict % getOrDefault(self % handleVirtual,'handleVirtual', .true.)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(mgXsClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill and deallocate maps
    call self % energyMap % kill()

    if (allocated(self % spaceMap)) then
      call self % spaceMap % kill()
    end if

    ! Reset parameters
    self % matN    = 0
    self % energyN = 0
    self % width   = 0
    self % PN      = .false.
    self % handleVirtual = .false.

  end subroutine kill

  !!
  !! Returns array of codes that represent different reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(mgXsClerk),intent(in)                :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [inColl_CODE, outColl_CODE, spawn_CODE]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(mgXsClerk), intent(in) :: self
    integer(shortInt)            :: S

    S = self % width * self % energyN * self % matN

  end function getSize

  !!
  !! Process incoming collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportInColl(self, p, xsData, mem, virtual)
    class(mgXsClerk), intent(inout)       :: self
    class(particle), intent(in)           :: p
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    logical(defBool), intent(in)          :: virtual
    type(particleState)                   :: state
    type(neutronMacroXSs)                 :: xss
    class(neutronMaterial), pointer       :: mat
    real(defReal)                         :: nuFissXS, captXS, fissXS, scattXS, flux
    integer(shortInt)                     :: enIdx, locIdx, binIdx
    integer(longInt)                      :: addr
    character(100), parameter :: Here =' reportInColl (mgXsClerk_class.f90)'

    ! Return if collision is virtual but virtual collision handling is off
    if ((.not. self % handleVirtual) .and. virtual) return

    ! Get current particle state
    state = p

    ! Find bin indexes
    ! Energy
    if (allocated(self % energyMap)) then
      enIdx = self % energyN + 1 - self % energyMap % map(state)
    else
      enIdx = 1
    end if
    ! Space
    if (allocated(self % spaceMap)) then
      locIdx = self % spaceMap % map(state)
    else
      locIdx = 1
    end if

    ! Return if invalid bin index
    if ((enIdx == self % energyN + 1) .or. locIdx == 0) return

    ! Calculate bin address
    binIdx = self % energyN * (locIdx - 1) + enIdx
    addr = self % getMemAddress() + self % width * (binIdx - 1) - 1

    ! Calculate flux with the right cross section according to virtual collision handling
    if (self % handleVirtual) then
      flux = p % w / xsData % getTrackingXS(p, p % matIdx(), TRACKING_XS)
    else
      flux = p % w / xsData % getTotalMatXS(p, p % matIdx())
    end if

    ! Check if the particle is in void. This call might happen when handling virtual collisions.
    ! This is relevant in the case of homogenising materials that include void: the flux
    ! in void will be different than zero, and the zero reaction rates have to be averaged
    if (p % matIdx() /= VOID_MAT) then

      ! Get material pointer
      mat => neutronMaterial_CptrCast(xsData % getMaterial(p % matIdx()))
      if (.not.associated(mat)) then
        call fatalError(Here,'Unrecognised type of material was retrived from nuclearDatabase')
      end if

      ! Retrieve material cross sections
      call mat % getMacroXSs(xss, p)

      ! Calculate reaction rates
      nuFissXS = xss % nuFission * flux
      captXS   = xss % capture * flux
      fissXS   = xss % fission * flux
      scattXS  = (xss % elasticScatter + xss % inelasticScatter) * flux

    else

      ! Reaction rates in void are zero
      nuFissXS = ZERO
      captXS   = ZERO
      fissXS   = ZERO
      scattXS  = ZERO

    end if

    ! Add scores to counters
    call mem % score(flux,     addr + FLUX_idx)
    call mem % score(nuFissXS, addr + NUBAR_idx)
    call mem % score(captXS,   addr + CAPT_idx)
    call mem % score(fissXS,   addr + FISS_idx)
    call mem % score(scattXS,  addr + SCATT_idx)

  end subroutine reportInColl

  !!
  !! Process outgoing collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportOutColl(self, p, MT, muL, xsData, mem)
    class(mgXsClerk), intent(inout)      :: self
    class(particle), intent(in)          :: p
    integer(shortInt), intent(in)        :: MT
    real(defReal), intent(in)            :: muL
    class(nuclearDatabase),intent(inout) :: xsData
    type(scoreMemory), intent(inout)     :: mem
    type(particleState)                  :: preColl, postColl
    real(defReal)                        :: score, prod, mu, mu2, mu3, mu4, mu5
    integer(shortInt)                    :: enIdx, locIdx, binIdx, binEnOut
    integer(longInt)                     :: addr

    ! Get pre and post collision particle state
    preColl  = p % preCollision
    postColl = p

    ! Find multipliers to account for scattering multiplicity
    select case(MT)
      case(N_2N, N_2Na, N_2Nd, N_2Nf, N_2Np, N_2N2a, N_2Nl(1):N_2Nl(16))
        score = 2.0_defReal
      case(N_3N, N_3Na, N_3Nf, N_3Np)
        score = 3.0_defReal
      case(N_4N)
        score = 4.0_defReal
      case default
        score = ONE
    end select

    ! Score in case of scattering events
    select case(MT)
      case ( N_N_ELASTIC, N_N_INELASTIC, N_N_ThermINEL, N_Nl(1):N_Nl(40), N_Ncont, &
             N_2N, N_2Na, N_2Nd, N_2Nf, N_2Np, N_2N2a, N_2Nl(1):N_2Nl(16),  &
             N_3N, N_3Na, N_3Nf, N_3Np, N_4N, N_Na, N_Np, N_Nd, N_Nt)

        ! Find bin indexes
        ! Energy
        if (allocated(self % energyMap)) then
          enIdx = self % energyN + 1 - self % energyMap % map(preColl)
        else
          enIdx = 1
        end if
        ! Space
        if (allocated(self % spaceMap)) then
          locIdx = self % spaceMap % map(preColl)
        else
          locIdx = 1
        end if

        ! Return if invalid bin index
        if ((enIdx == self % energyN + 1) .or. locIdx == 0) return

        ! Calculate bin address
        binIdx = self % energyN * (locIdx - 1) + enIdx
        addr = self % getMemAddress() + self % width * (binIdx - 1) - 1

        ! Score a scattering event from group g
        call mem % score(preColl % wgt, addr + SCATT_EV_idx)

        ! Get bin of outgoing energy
        if (allocated(self % energyMap)) then
          binEnOut = self % energyN + 1 - self % energyMap % map(postColl)
        else
          binEnOut = 1
        end if

        ! Return if invalid bin index
        if (binEnOut == self % energyN + 1) return

        ! Score scattering event from group g to g'
        call mem % score(preColl % wgt, addr + SCATT_EV_idx + binEnOut)

        ! Score outgoing scattering angle for P1 matrix
        mu = muL * preColl % wgt
        call mem % score(mu, addr + SCATT_EV_idx + self % energyN + binEnOut)

        ! Score multiplicity matrix
        prod = score * preColl % wgt
        call mem % score(prod, addr + SCATT_EV_idx + 2*self % energyN + binEnOut)

        ! Higher order scattering matrices
        if (self % PN) then

          ! Pre-computing powers
          mu2 = muL * muL
          mu3 = mu2 * muL
          mu4 = mu3 * muL
          mu5 = mu4 * muL

          ! Score outgoing scattering angle for P2 matrix
          mu = HALF * (3.0_defReal*mu2 - ONE) * preColl % wgt
          call mem % score(mu, addr + SCATT_EV_idx + 3*self % energyN + binEnOut)

          ! Score outgoing scattering angle for P3 matrix
          mu = HALF * (5.0_defReal*mu3 - 3.0_defReal*muL) * preColl % wgt
          call mem % score(mu, addr + SCATT_EV_idx + 4*self % energyN + binEnOut)

          ! Score outgoing scattering angle for P4 matrix
          mu = (35.0_defReal*mu4 - 30.0_defReal*mu2 + 3.0_defReal)/8.0_defReal * preColl % wgt
          call mem % score(mu, addr + SCATT_EV_idx + 5*self % energyN + binEnOut)

          ! Score outgoing scattering angle for P5 matrix
          mu = (63.0_defReal*mu5 - 70.0_defReal*mu3 + 15.0_defReal*muL)/8.0_defReal * preColl % wgt
          call mem % score(mu, addr + SCATT_EV_idx + 6*self % energyN + binEnOut)

          ! Score outgoing scattering angle for P6 matrix
          mu = (231.0_defReal*mu5*muL - 315.0_defReal*mu4 + 105.0_defReal*mu2 &
                - 5.0_defReal)/16.0_defReal * preColl % wgt
          call mem % score(mu, addr + SCATT_EV_idx + 7*self % energyN + binEnOut)

          ! Score outgoing scattering angle for P7 matrix
          mu = (429.0_defReal*mu5*mu2 - 693.0_defReal*mu5 + 315.0_defReal*mu3 &
                - 35.0_defReal*muL)/16.0_defReal * preColl % wgt
          call mem % score(mu, addr + SCATT_EV_idx + 8*self % energyN + binEnOut)

        end if

      case default
        ! Do nothing

    end select

  end subroutine reportOutColl

  !!
  !! Process fission report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportSpawn(self, MT, pOld, pNew, xsData, mem)
    class(mgXsClerk), intent(inout)       :: self
    integer(shortInt), intent(in)         :: MT
    class(particle), intent(in)           :: pOld
    class(particleState), intent(in)      :: pNew
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    integer(longInt)                      :: addr, binIdx, enIdx, locIdx

    if (MT == N_FISSION) then

      ! Find bin indexes
      ! Energy
      if (allocated(self % energyMap)) then
        enIdx = self % energyN + 1 - self % energyMap % map(pNew)
      else
        enIdx = 1
      end if
      ! Space
      if (allocated(self % spaceMap)) then
        locIdx = self % spaceMap % map(pNew)
      else
        locIdx = 1
      end if

      ! Return if invalid bin index
      if ((enIdx == self % energyN + 1) .or. locIdx == 0) return

      ! Calculate bin address
      binIdx = self % energyN * (locIdx - 1) + enIdx
      addr = self % getMemAddress() + self % width * (binIdx - 1) - 1

      ! Score energy group of fission neutron
      call mem % score(ONE,  addr + CHI_idx)

    end if

  end subroutine reportSpawn

  !!
  !! Final processing to calculate the multi-group cross sections, fission data and
  !! P0, P1 and production scattering matrices
  !!
  !! Args:
  !!   mem [in] -> memory
  !!   capt_res [out]    -> capture MG xss and uncertainties
  !!   fiss_res [out]    -> fission MG xss and uncertainties
  !!   transFL_res [out] -> transport MG xss and uncertainties with flux limited approximation
  !!   transOS_res [out] -> transport MG xss and uncertainties with out-scatter approximation
  !!   nu_res   [out]    -> MG nuBar and uncertainties
  !!   chi_res  [out]    -> MG fission spectrum and uncertainties
  !!   P0_res   [out]    -> P0 scattering matrix and uncertainties
  !!   P1_res   [out]    -> P1 scattering matrix and uncertainties
  !!   prod_res [out]    -> P0 scattering production matrix and uncertainties
  !!
  !! Errors:
  !!   none
  !!
  pure subroutine processRes(self, mem, capt_res, fiss_res, transFL_res, transOS_res, &
                             nu_res, chi_res, P0_res, P1_res, prod_res)
    class(mgXsClerk), intent(in)    :: self
    type(scoreMemory), intent(in)   :: mem
    real(defReal), dimension(:,:), allocatable, intent(out) :: capt_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: fiss_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: transFL_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: transOS_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: nu_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: chi_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: P0_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: P1_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: prod_res
    real(defReal), dimension(:,:), allocatable  :: delta, deltaStd
    real(defReal), dimension(:), allocatable    :: tot, totStd, fluxG, fluxGstd
    integer(longInt)  :: addr
    integer(shortInt) :: N, M, i, j, k, g1, gEnd, idx
    real(defReal)     :: capt, fiss, scatt, nu, chi, P0, P1, prod, sumChi, flux, scattProb, &
                         captStd, fissStd, scattStd, nuStd, chiStd, P0std, P1std, prodStd,  &
                         fluxStd, scattProbStd, scattXS, scattXSstd

    ! Get number of bins
    N = self % energyN
    M = self % matN

    ! Allocate arrays for MG xss
    allocate( capt_res(2,N*M), fiss_res(2,N*M), transFL_res(2,N*M), transOS_res(2,N*M), &
              nu_res(2,N*M), chi_res(2,N*M), P0_res(2,N*N*M), P1_res(2,N*N*M),          &
              prod_res(2,N*N*M), tot(N*M), fluxG(N*M), delta(M, N), totStd(N*M),        &
              fluxGstd(N*M), deltaStd(M, N) )

    ! Initialise values
    sumChi = 0    ! to normalise chi
    k      = 1    ! to calculate transport xss
    delta  = ZERO
    deltaStd = ZERO

    ! Loop over all energies and spatial bins
    do i = 1, N*M

      ! Retrieve reaction rates and flux from memory
      addr = self % getMemAddress() + self % width * (i - 1) - 1
      call mem % getResult(flux,  fluxStd,  addr + FLUX_idx)
      call mem % getResult(fiss,  fissStd,  addr + FISS_idx)
      call mem % getResult(capt,  captStd,  addr + CAPT_idx)
      call mem % getResult(scatt, scattStd, addr + SCATT_idx)
      call mem % getResult(nu,    nuStd,    addr + NUBAR_idx)
      call mem % getResult(chi,   chiStd,   addr + CHI_idx)
      call mem % getResult(scattProb, scattProbStd, addr + SCATT_EV_idx)

      ! Calculate MG constants, being careful to avoid division by zero
      ! If flux is zero all cross sections must be set to zero
      if (flux == ZERO) then
        capt_res(1:2,i) = ZERO
        scattXS         = ZERO
        scattXSstd      = ZERO
      else
        capt_res(1,i) = capt/flux
        capt_res(2,i) = capt_res(1,i) * sqrt((captStd/capt)**2 + (fluxStd/flux)**2)
        scattXS       = scatt/flux
        scattXSstd    = scattXS * sqrt((scattStd/scatt)**2 + (fluxStd/flux)**2)
      end if

      ! Calculate fission production term and uncertainties
      if (fiss == ZERO) then
        fiss_res(1:2,i) = ZERO
        nu_res(1:2,i)   = ZERO
      else
        fiss_res(1,i) = fiss/flux
        fiss_res(2,i) = fiss_res(1,i) * sqrt((fissStd/fiss)**2 + (fluxStd/flux)**2)
        nu_res(1,i)   = nu/fiss
        nu_res(2,i)   = nu_res(1,i) * sqrt((nuStd/nu)**2 + (fissStd/fiss)**2)
      end if

      ! Store fission spectrum
      chi_res(1,i) = chi
      chi_res(2,i) = chiStd
      sumChi = sumChi + chi
      ! If this is the last energy group for a material, normalise the spectrum
      if (mod(i,N) == 0) then
        g1   = i+1-N
        gEnd = i
        if (sumChi /= ZERO) chi_res(1:2, g1:gEnd) = chi_res(1:2, g1:gEnd)/sumChi
        sumChi = 0
      end if

      ! Store total cross section and flux for this energy group
      tot(i)    = capt_res(1,i) + fiss_res(1,i) + scattXS
      totStd(i) = sqrt(capt_res(2,i)**2 + fiss_res(2,i)**2 + scattXSstd**2)
      fluxG(i)  = flux
      fluxGstd(i) = fluxStd

      ! Loop over outgoing energies to calculate the scattering matrices
      do j = 1, N

        ! Retrieve results from memory
        call mem % getResult(P0, P0std, addr + SCATT_EV_idx + j)
        call mem % getResult(P1, P1std, addr + SCATT_EV_idx + N + j)
        call mem % getResult(prod, prodStd, addr + SCATT_EV_idx + 2*N + j)

        ! Calculate scattering matrices and uncertainties
        idx = N*(i-1) + j
        if (P0 == ZERO) then
          P0_res(1:2, idx) = ZERO
          P1_res(1:2, idx) = ZERO
          prod_res(1, idx) = ONE
          prod_res(2, idx) = ZERO
        else
          P0_res(1, idx) = P0/scattProb*scattXS
          P0_res(2, idx) = P0_res(1,idx) * sqrt((P0std/P0)**2 + (scattProbStd/scattProb)**2 + (scattXSstd/scattXS)**2)
          P1_res(1, idx) = P1/scattProb*scattXS
          P1_res(2, idx) = P1_res(1,idx) * sqrt((P1std/P1)**2 + (scattProbStd/scattProb)**2 + (scattXSstd/scattXS)**2)
          prod_res(1, idx) = prod/P0
          prod_res(2, idx) = prod_res(1,idx) * sqrt((prodStd/prod)**2 + (P0std/P0)**2)
          ! Accumulate quantities over all outgoing groups for transport cross section generation
          delta(k,j)    = delta(k,j) + P1/scattProb*scatt
          deltaStd(k,j) = sqrt(deltaStd(k,j)**2 + ((P1std/P1)**2 + (scattProbStd/scattProb)**2 + &
                          (scattStd/scatt)**2) * (P1/scattProb*scatt)**2)
        end if

      end do

      ! If this is the last energy group for a material increment index
      if (mod(i, N) == 0) k = k + 1

      ! Calculate out-scatter transport cross section
      g1   = N*(i-1) + 1
      gEnd = N*(i-1) + N
      transOS_res(1,i) = tot(i) - sum(P1_res(1, g1:gEnd))
      transOS_res(2,i) = sqrt(totStd(i)**2 + sum(P1_res(2, g1:gEnd)))

    end do

    ! Calculate flux-limited transport cross section
    do i = 1, M
      do j = 1, N
        idx = N*(i-1) + j
        if (fluxG(idx) == ZERO) then
          transFL_res(1:2,idx) = ZERO
        else
          transFL_res(1,idx) = tot(idx) - delta(i,j)/fluxG(idx)
          transFL_res(2,idx) = sqrt(totStd(idx)**2 + ((deltaStd(i,j)/delta(i,j))**2 + &
                                    (fluxGstd(idx)/fluxG(idx))**2) * (delta(i,j)/fluxG(idx))**2)
        end if
      end do
    end do

  end subroutine processRes

  !!
  !! Final processing to calculate the high order scattering matrices, P2 to P7
  !!
  !! Args:
  !!   mem [in] -> memory
  !!   P2_res [out] -> P2 scattering matrix and uncertainties
  !!   P3_res [out] -> P3 scattering matrix and uncertainties
  !!   P4_res [out] -> P4 scattering matrix and uncertainties
  !!   P5_res [out] -> P5 scattering matrix and uncertainties
  !!   P6_res [out] -> P6 scattering matrix and uncertainties
  !!   P7_res [out] -> P7 scattering matrix and uncertainties
  !!
  !! Errors:
  !!   none
  !!
  pure subroutine processPN(self, mem, P2_res, P3_res, P4_res, P5_res, P6_res, P7_res)
    class(mgXsClerk), intent(in)     :: self
    type(scoreMemory), intent(in)    :: mem
    real(defReal), dimension(:,:), allocatable, intent(out) :: P2_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: P3_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: P4_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: P5_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: P6_res
    real(defReal), dimension(:,:), allocatable, intent(out) :: P7_res
    integer(longInt)  :: addr
    integer(shortInt) :: N, M, i, j, idx
    real(defReal)     :: scatt, scattXS, scattProb, flux, scattStd, scattXSstd, &
                         scattProbStd, fluxStd, P0, P2, P3, P4, P5, P6, P7,     &
                         P0std, P2std, P3std, P4std, P5std, P6std, P7std

    ! Get number of bins
    N = self % energyN
    M = self % matN

    ! Allocate arrays for scattering matrices
    allocate( P2_res(2,N*N*M), P3_res(2,N*N*M), P4_res(2,N*N*M), &
              P5_res(2,N*N*M), P6_res(2,N*N*M), P7_res(2,N*N*M) )

    ! Loop over all energies and spatial bins
    do i = 1, N*M

      ! Retrieve reaction rates and flux from memory
      addr = self % getMemAddress() + self % width * (i - 1) - 1
      call mem % getResult(scatt, scattStd, addr + SCATT_idx)
      call mem % getResult(flux,  fluxStd,  addr + FLUX_idx)
      call mem % getResult(scattProb, scattProbStd, addr + SCATT_EV_idx)

      ! Calculate MG constants, being careful to avoid division by zero
      ! If flux is zero all cross sections must be set to zero
      if (flux == ZERO) then
        scattXS    = ZERO
        scattXSstd = ZERO
      else
        scattXS    = scatt/flux
        scattXSstd = scattXS * sqrt((scattStd/scatt)**2 + (fluxStd/flux)**2)
      end if

      ! Loop over outgoing energies to calculate the scattering matrices
      do j = 1, N

        ! Retrieve results from memory
        call mem % getResult(P0, P0std, addr + SCATT_EV_idx + j)
        call mem % getResult(P2, P2std, addr + SCATT_EV_idx + 3*N + j)
        call mem % getResult(P3, P3std, addr + SCATT_EV_idx + 4*N + j)
        call mem % getResult(P4, P4std, addr + SCATT_EV_idx + 5*N + j)
        call mem % getResult(P5, P5std, addr + SCATT_EV_idx + 6*N + j)
        call mem % getResult(P6, P6std, addr + SCATT_EV_idx + 7*N + j)
        call mem % getResult(P7, P7std, addr + SCATT_EV_idx + 8*N + j)

        ! Calculate scattering matrices and uncertainties
        idx = N*(i-1) + j
        if (P0 == ZERO) then
          P2_res(1:2,idx) = ZERO
          P3_res(1:2,idx) = ZERO
          P4_res(1:2,idx) = ZERO
          P5_res(1:2,idx) = ZERO
          P6_res(1:2,idx) = ZERO
          P7_res(1:2,idx) = ZERO
        else
          P2_res(1,idx) = P2/scattProb*scattXS
          P2_res(2,idx) = P2_res(1,idx) * sqrt((P2std/P2)**2 + (scattProbStd/scattProb)**2 + (scattXSstd/scattXS)**2)
          P3_res(1,idx) = P3/scattProb*scattXS
          P3_res(2,idx) = P3_res(1,idx) * sqrt((P3std/P3)**2 + (scattProbStd/scattProb)**2 + (scattXSstd/scattXS)**2)
          P4_res(1,idx) = P4/scattProb*scattXS
          P4_res(2,idx) = P4_res(1,idx) * sqrt((P4std/P4)**2 + (scattProbStd/scattProb)**2 + (scattXSstd/scattXS)**2)
          P5_res(1,idx) = P5/scattProb*scattXS
          P5_res(2,idx) = P5_res(1,idx) * sqrt((P5std/P5)**2 + (scattProbStd/scattProb)**2 + (scattXSstd/scattXS)**2)
          P6_res(1,idx) = P6/scattProb*scattXS
          P6_res(2,idx) = P6_res(1,idx) * sqrt((P6std/P6)**2 + (scattProbStd/scattProb)**2 + (scattXSstd/scattXS)**2)
          P7_res(1,idx) = P7/scattProb*scattXS
          P7_res(2,idx) = P7_res(1,idx) * sqrt((P7std/P7)**2 + (scattProbStd/scattProb)**2 + (scattXSstd/scattXS)**2)
        end if

      end do

    end do

  end subroutine processPN

  !!
  !! Write contents of the clerk to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(mgXsClerk), intent(in)               :: self
    class(outputFile), intent(inout)           :: outFile
    type(scoreMemory), intent(in)              :: mem
    integer(shortInt),dimension(:),allocatable :: resArrayShape
    real(defReal), dimension(:,:), allocatable :: fiss, capt, transFL, transOS, &
                                                  nu, chi, P0, P1, P2, P3, P4,  &
                                                  P5, P6, P7, prod
    character(nameLen)                         :: name
    integer(shortInt)                          :: i

    ! Begin block
    call outFile % startBlock(self % getName())

    ! Allocate space for resultShape array
    if (allocated(self % spaceMap)) then
      allocate(resArrayShape(self % spaceMap % dimensions() + 1))
    else
      allocate(resArrayShape(1))
    end if

    ! Print energy map information
    if (allocated(self % energyMap)) call self % energyMap % print(outFile)
    resArrayShape(1) = self % energyN

    ! If a space map print map information
    if (allocated(self % spaceMap)) then
      call self % spaceMap % print(outFile)
      resArrayShape(2:(self % spaceMap % dimensions() + 1)) = self % spaceMap % binArrayShape()
    end if

    ! Process and get results
    call self % processRes(mem, capt, fiss, transFL, transOS, nu, chi, P0, P1, prod)

    ! Print results
    name = 'capture'
    call outFile % startArray(name, resArrayShape)
    do i=1,product(resArrayShape)
      call outFile % addResult(capt(1,i),capt(2,i))
    end do
    call outFile % endArray()

    name = 'fission'
    call outFile % startArray(name, resArrayShape)
    do i=1,product(resArrayShape)
      call outFile % addResult(fiss(1,i),fiss(2,i))
    end do
    call outFile % endArray()

    name = 'transportFluxLimited'
    call outFile % startArray(name, resArrayShape)
    do i=1,product(resArrayShape)
      call outFile % addResult(transFL(1,i),transFL(2,i))
    end do
    call outFile % endArray()

    name = 'transportOutScatter'
    call outFile % startArray(name, resArrayShape)
    do i=1,product(resArrayShape)
      call outFile % addResult(transOS(1,i),transOS(2,i))
    end do
    call outFile % endArray()

    name = 'nu'
    call outFile % startArray(name, resArrayShape)
    do i=1,product(resArrayShape)
      call outFile % addResult(nu(1,i),nu(2,i))
    end do
    call outFile % endArray()

    name = 'chi'
    call outFile % startArray(name, resArrayShape)
    do i=1,product(resArrayShape)
      call outFile % addResult(chi(1,i),chi(2,i))
    end do
    call outFile % endArray()

    ! Modify the shape of the result array to print matrices
    resArrayShape(1) = resArrayShape(1) * self % energyN
    name = 'P0'
    call outFile % startArray(name, resArrayShape)
    do i=1,product(resArrayShape)
      call outFile % addResult(P0(1,i),P0(2,i))
    end do
    call outFile % endArray()

    name = 'P1'
    call outFile % startArray(name, resArrayShape)
    do i=1,product(resArrayShape)
      call outFile % addResult(P1(1,i),P1(2,i))
    end do
    call outFile % endArray()

    name = 'prod'
    call outFile % startArray(name, resArrayShape)
    do i=1,product(resArrayShape)
      call outFile % addResult(prod(1,i),prod(2,i))
    end do
    call outFile % endArray()

    ! Deallocate to limit memory consumption when writing to the output file
    deallocate(capt, fiss, transFL, transOS, nu, chi, P0, P1, prod)

    ! If high order scattering is requested, print the other matrices
    if (self % PN) then

      call self % processPN(mem, P2, P3, P4, P5, P6, P7)

      name = 'P2'
      call outFile % startArray(name, resArrayShape)
      do i = 1,product(resArrayShape)
        call outFile % addResult(P2(1,i),P2(2,i))
      end do
      call outFile % endArray()

      name = 'P3'
      call outFile % startArray(name, resArrayShape)
      do i = 1,product(resArrayShape)
        call outFile % addResult(P3(1,i),P3(2,i))
      end do
      call outFile % endArray()

      name = 'P4'
      call outFile % startArray(name, resArrayShape)
      do i = 1,product(resArrayShape)
        call outFile % addResult(P4(1,i),P4(2,i))
      end do
      call outFile % endArray()

      name = 'P5'
      call outFile % startArray(name, resArrayShape)
      do i = 1,product(resArrayShape)
        call outFile % addResult(P5(1,i),P5(2,i))
      end do
      call outFile % endArray()

      name = 'P6'
      call outFile % startArray(name, resArrayShape)
      do i = 1,product(resArrayShape)
        call outFile % addResult(P6(1,i),P6(2,i))
      end do
      call outFile % endArray()

      name = 'P7'
      call outFile % startArray(name, resArrayShape)
      do i = 1,product(resArrayShape)
        call outFile % addResult(P7(1,i),P7(2,i))
      end do
      call outFile % endArray()

      ! Deallocate to limit memory consumption when writing to the output file
      deallocate(P2, P3, P4, P5, P6, P7)

    end if

    call outFile % endBlock()

  end subroutine print

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(mgXsClerk), intent(in)  :: self
    type(scoreMemory), intent(in) :: mem

    print *, 'mgXsClerk does not support display yet'

  end subroutine display

end module mgXsClerk_class
