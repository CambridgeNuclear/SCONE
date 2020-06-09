module fissionSource_class

  use numPrecision
  use endfConstants
  use universalVariables,      only : OUTSIDE_MAT, VOID_MAT
  use genericProcedures,       only : fatalError, rotateVector
  use dictionary_class,        only : dictionary
  use RNG_class,               only : RNG

  use particle_class,          only : particleState, P_NEUTRON
  use source_inter,            only : source, kill_super => kill

  use geometry_inter,          only : geometry
  use neutronMaterial_inter,   only : neutronMaterial, neutronMaterial_CptrCast
  use ceNeutronMaterial_class, only : ceNeutronMaterial, ceNeutronMaterial_CptrCast
  use fissionCE_class,         only : fissionCE, fissionCE_TptrCast
  use fissionMG_class,         only : fissionMG, fissionMG_TptrCast
  use nuclearDataReg_mod,      only : ndReg_getNeutronCE => getNeutronCE, &
                                      ndReg_getNeutronMG => getNeutronMG
  use nuclearDatabase_inter,   only : nuclearDatabase
  use ceNeutronDatabase_inter, only : ceNeutronDatabase
  use mgNeutronDatabase_inter, only : mgNeutronDatabase

  implicit none
  private

  !!
  !! Neutron Source from distributed fission sites
  !!
  !! Places fission sites uniformly in regions with fissile material.
  !! Spectrum of the fission neutron is such as if it fission was caused by incdent
  !! neutron with CE energy E or MG with group G.
  !! Angular distribution is isotropic.
  !!
  !! Private members:
  !!   isMG   -> is the source multi-group? (default = .false.)
  !!   bottom -> Bottom corner (x_min, y_min, z_min)
  !!   top    -> Top corner (x_max, y_max, z_max)
  !!   E      -> Fission site energy [MeV] (default = 1.0E-6)
  !!   G      -> Fission site Group (default = 1)
  !!
  !! Interface:
  !!   source_inter Interface
  !!
  !! Sample Dictionary Input:
  !!   fission {
  !!     type fissionSource;
  !!     #data MG; #
  !!     #E 15.0;  #
  !!     #G 7;     #
  !!   }
  !!
  type, public,extends(source) :: fissionSource
    private
    logical(defBool)            :: isMG   = .false.
    real(defReal), dimension(3) :: bottom = ZERO
    real(defReal), dimension(3) :: top    = ZERO
    real(defReal)               :: E      = ZERO
    integer(shortInt)           :: G      = 0
  contains
    procedure :: init
    procedure :: sampleParticle
    procedure :: kill
  end type fissionSource

contains

  !!
  !! Initialise Fission Source
  !!
  !! See source_inter for details
  !!
  subroutine init(self, dict, geom)
    class(fissionSource), intent(inout)      :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    character(nameLen)                       :: type
    real(defReal), dimension(6)              :: bounds
    character(100), parameter :: Here = 'init (fissionSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    ! Select Energy Type
    call dict % getOrDefault(type, 'data', 'ce')
    select case(type)
      case('ce')
        self % isMG = .false.

      case('mg')
        self % isMG = .true.

      case default
        call fatalError(Here, 'Invalid source data type specified: must be ce or mg')
    end select

    ! Select Required fission group/energy
    call dict % getOrDefault(self % E, 'E', 1.0E-6_defReal)
    call dict % getOrDefault(self % G, 'G', 1)

    ! Set bounding region
    bounds = self % geom % bounds()
    self % bottom = bounds([1, 3, 5])
    self % top    = bounds([2, 4, 6])

  end subroutine init

  !!
  !! Sample particle's phase space co-ordinates
  !!
  !! See source_inter for details
  !!
  function sampleParticle(self, rand) result(p)
    class(fissionSource), intent(inout)  :: self
    class(RNG), intent(inout)            :: rand
    type(particleState)                  :: p
    class(nuclearDatabase), pointer      :: nucData
    class(neutronMaterial), pointer      :: mat
    class(ceNeutronMaterial), pointer    :: matCE
    type(fissionCE), pointer             :: fissCE
    type(fissionMG), pointer             :: fissMG
    real(defReal), dimension(3)          :: r, rand3
    real(defReal)                        :: mu, phi, E_out, E_up, E_down
    integer(shortInt)                    :: matIdx, uniqueID, nucIdx, i, G_out
    character(100), parameter :: Here = 'sampleParticle (fissionSource_class.f90)'

    ! Get pointer to approperiate nuclear database
    if (self % isMG) then
      nucData => ndReg_getNeutronMG()
    else
      nucData => ndReg_getNeutronCE()
    end if
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

    i = 0
    rejection : do
      ! Protect against infinate loop
      i = i +1
      if ( i > 200) then
        call fatalError(Here, 'Infinate loop in sampling of fission sites. Please check that&
                              & defined volume contains fissile material.')
      end if

      ! Sample Position
      rand3(1) = rand % get()
      rand3(2) = rand % get()
      rand3(3) = rand % get()
      r = (self % top - self % bottom) * rand3 + self % bottom

      ! Find material under position
      call self % geom % whatIsAt(r, matIdx, uniqueID)

      ! Reject if there is no material
      if (matIdx == VOID_MAT .or. matIdx == OUTSIDE_MAT) cycle rejection

      mat => neutronMaterial_CptrCast(nucData % getMaterial(matIdx))
      if (.not.associated(mat)) call fatalError(Here, "Nuclear data did not return neutron material.")

      ! Resample position if material is not fissile
      if (.not.mat % isFissile()) cycle

      ! Assign basic phase-space coordinates
      p % matIdx   = matIdx
      p % uniqueID = uniqueID
      p % wgt      = ONE
      p % time     = ZERO
      p % type     = P_NEUTRON
      p % r        = r

      ! Set Energy
      select type (nucData)
        class is (ceNeutronDatabase)
          ! Get energy bounds
          call nucData% energyBounds(E_down, E_up)

          ! Get material
          matCE => ceNeutronMaterial_CptrCast(nucData % getMaterial(matIdx))
          if (.not.associated(mat)) then
            call fatalError(Here, 'Failed to get ceNeutronMaterial from ceDatabase.')
          end if

          ! Get Nuclide
          nucIdx = matCE % sampleFission(self % E, rand)

          ! Get reaction object
          fissCE => fissionCE_TptrCast(nucData % getReaction(N_FISSION, nucIdx))
          if(.not.associated(fissCE)) then
            call fatalError(Here, "Failed to get CE Fission Reaction Object")
          end if

          ! Get mu, phi, E_out
          call fissCE % sampleOut(mu, phi, E_out, self % E, rand)

          ! Put data into particle
          p % E = E_out
          p % isMG = .false.
          p % dir  = rotateVector([ONE, ZERO, ZERO], mu, phi)

          ! Apply upper energy cut-off
          if (p % E > E_up) p % E = E_up

        class is (mgNeutronDatabase)
          ! Get reaction object
          fissMG => fissionMG_TptrCast(nucData % getReaction(macroFission, matIdx))
          if(.not.associated(fissMG)) then
            call fatalError(Here, "Failed to get MG Fission Reaction Object")
          end if

          ! Get mu, phi, G_out
          call fissMG % sampleOut(mu, phi, G_out, self % G, rand)

          ! Set outgoing state
          p % G = G_out
          p % isMG = .true.
          p % dir = rotateVector([ONE, ZERO, ZERO], mu, phi)

        class default
          call fatalError(Here, "Uncrecognised type of nuclearDatabase")

      end select
      ! Exit the loop
      exit rejection

    end do rejection

  end function sampleParticle

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(fissionSource), intent(inout) :: self

    call kill_super(self)

    self % isMG   = .false.
    self % bottom = ZERO
    self % top    = ZERO
    self % E      = ZERO
    self % G      = 0

  end subroutine kill

end module fissionSource_class
