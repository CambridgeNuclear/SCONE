module fissionSource_class

  use numPrecision
  use endfConstants
  use universalVariables,      only : OUTSIDE_MAT, VOID_MAT, UNDEF_MAT, OVERLAP_MAT
  use genericProcedures,       only : rotateVector, numToChar
  use errors_mod,              only : fatalError
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
  !! Spectrum of the fission neutron is such as if it fission was caused by
  !! incident neutron with CE energy E or MG with group G.
  !! Angular distribution is isotropic.
  !!
  !! Private members:
  !!   isMG     -> is the source multi-group? (default = .false.)
  !!   bottom   -> Bottom corner (x_min, y_min, z_min)
  !!   top      -> Top corner (x_max, y_max, z_max)
  !!   E        -> Fission site energy [MeV] (default = 1.0E-6)
  !!   G        -> Fission site Group (default = 1)
  !!   attempts -> Maximum number of attempts to find a fissile material
  !! Interface:
  !!   source_inter Interface
  !!
  !! Sample Dictionary Input:
  !!   fission {
  !!     type fissionSource;
  !!     #data MG; #
  !!     #E 15.0;  #
  !!     #G 7;     #
  !!     #attempts 100000; # (default: 10000)
  !!     #top (1.0 1.0 1.0);    #
  !!     #bottom (0.0 0.0 0.0); #
  !!   }
  !!
  type, public,extends(source) :: fissionSource
    private
    logical(defBool)            :: isMG   = .false.
    real(defReal), dimension(3) :: bottom = ZERO
    real(defReal), dimension(3) :: top    = ZERO
    real(defReal)               :: E      = ZERO
    integer(shortInt)           :: G      = 0
    integer(shortInt)           :: attempts = 10000
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
    real(defReal), allocatable, dimension(:) :: temp
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

    ! Load the number of allowed attempts
    call dict % getOrDefault(self % attempts, 'attempts', 10000)

    if (self % attempts < 1) then
      call fatalError(Here, 'Number of attempts must be greater than 0. Is: ' // numToChar(self % attempts))
    end if

    ! Set bounding region
    if (dict % isPresent('top') .or. dict % isPresent('bottom')) then
      ! Get top and bottom from dictionary
      call dict % get(temp, 'top')
      if (size(temp) /= 3) then
        call fatalError(Here, "Top point must have 3 coordinates.")
      end if
      self % top = temp

      call dict % get(temp, 'bottom')
      if (size(temp) /= 3) then
        call fatalError(Here, "Bottom point must have 3 coordinates.")
      end if
      self % bottom = temp

      ! Check that top and bottom are valid
      if (any(self % top < self % bottom)) then
        call fatalError(Here, "Top point must have all coordinates greater than bottom point.")
      end if

    else ! Get bounds from geometry
      bounds = self % geom % bounds()
      self % bottom = bounds(1:3)
      self % top    = bounds(4:6)
    end if


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

    ! Get pointer to appropriate nuclear database
    if (self % isMG) then
      nucData => ndReg_getNeutronMG()
    else
      nucData => ndReg_getNeutronCE()
    end if
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

    i = 0
    rejection : do
      ! Protect against infinite loop
      i = i +1
      if (i > self % attempts) then
        call fatalError(Here, "Failed to find a fissile material in: "// numToChar(self % attempts) // " attempts.&
                              & Increase the number of maximum attempts or verify that fissile materials are present.")
      end if

      ! Sample Position
      rand3(1) = rand % get()
      rand3(2) = rand % get()
      rand3(3) = rand % get()
      r = (self % top - self % bottom) * rand3 + self % bottom

      ! Find material under position
      call self % geom % whatIsAt(matIdx, uniqueID, r)

      ! Reject if there is no material
      if (matIdx == VOID_MAT .or. matIdx == OUTSIDE_MAT) cycle rejection

      ! Terminate if there is an error in the geometry
      if (matIdx == UNDEF_MAT) then
        print *, r
        call fatalError(Here, 'Particle position was sampled in an undefined material')
      elseif (matIdx == OVERLAP_MAT) then
        print *, r
        call fatalError(Here, 'Particle position was sampled in an overlapping cell region')
      end if

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
          call fatalError(Here, "Unrecognised type of nuclearDatabase")

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
