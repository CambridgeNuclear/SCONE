module generalSource_class

  use numPrecision
  use universalVariables,      only: OUTSIDE_MAT, NOT_FOUND
  use genericProcedures,       only: fatalError
  use dictionary_class,        only: dictionary
  use RNG_class,               only: RNG  

  use particle_class,          only: particleState, P_NEUTRON, P_PHOTON
  use source_inter,            only: source
  
  use geometry_inter,          only: geometry
  use materialMenu_mod,        only: matIdx
  use nuclearDatabase_inter
  use fissionCE_class,         only: fissionCE, fissionCE_TptrCast
  use fissionMG_class,         only: fissionMG, fissionMG_TptrCast
  use ceNeutronDatabase_inter, only: ceNeutronDatabase, ceNeutronDatabase_CptrCast
  use reactionHandle_inter,    only: reactionHandle
  use uncorrelatedReactionCE_inter, only: uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use reactionMG_inter,        only: reactionMG_CptrCast, reactionMG
  implicit none
  private

  !!
  !! Class describing a general particle source
  !!
  !! Can sample particles from a chosen material, cell, or point
  !! Providing a bounding box can increase sampling efficiency
  !! Particle energy can be fixed or sampled from a particular reaction
  !! Particle direction is isotropic or mono-directional
  !! Particles have a single specified type
  !!
  !! Private members:
  !!   r            -> optional source position
  !!   dir          -> optional source direction
  !!   E            -> optional source energy
  !!   G            -> optional source energy group
  !!   particleType -> source particle type
  !!   isMG         -> is the source multi-group?
  !!   isIsotropic  -> is the source isotropic?
  !!   reaction     -> optional reaction to sample particle from
  !!   isotope      -> optional isotope to fully specify the reaction source
  !!   cell         -> optional unique cell ID in which source resides
  !!   mat          -> optional material index in which source resides
  !!   bb           -> optional bounding box to improve sampling efficiency.
  !!                   If not provided as an input, use bounds from geometry.
  !!                   bb is an array defined: [x1 x2 y1 y2 z1 z2]
  !!
  !! Interface:
  !!   init              -> initialise point source
  !!   sampleType        -> set particle type
  !!   samplePosition    -> set particle position
  !!   sampleEnergy      -> set particle energy
  !!   sampleEnergyAngle -> sample particle angle
  !!   kill              -> terminate source
  !!
  type, public,extends(source) :: generalSource
    private
    real(defReal),dimension(:), allocatable  :: r
    real(defReal),dimension(:), allocatable  :: dir
    real(defReal)                            :: E
    integer(shortInt)                        :: G
    integer(shortInt)                        :: particleType
    logical(defBool)                         :: isMG = .false.
    logical(defBool)                         :: isIsotropic
    integer(shortInt)                        :: reaction
    integer(shortInt)                        :: isotope
    integer(shortInt)                        :: cell = -1
    integer(shortInt)                        :: mat = -1
    real(defReal), dimension(:), allocatable :: bb
    class(uncorrelatedReactionCE), pointer   :: reacCE => null()
    class(reactionMG), pointer               :: reacMG => null()
  contains
    procedure :: init
    procedure :: sampleType
    procedure :: samplePosition
    procedure :: sampleEnergy
    procedure :: sampleEnergyAngle
    procedure :: kill
  end type generalSource

contains

  !!
  !! Initialise general source
  !!
  !! Read dictionary to obtain source information and provide
  !! geometry to allow basic check or position sampling
  !!
  !! Args:
  !!   dict [in] -> dict containing point source information
  !!   geom [in] -> pointer to geometry, for checking that point source is inside geometry
  !!                or for sampling position if mat/cell are specified
  !!
  !! Result:
  !!   An initialised general source
  !!
  !! Errors:
  !!   - error if an unrecognised particle type is provided
  !!   - error if either direction or position have more than 3 components
  !!   - error if both CE and MG is specified
  !!   - error if position/cell/mat are specified with each other
  !!   - error if neither energy type is specified
  !!   - error if energy and reaction are specified
  !!
  subroutine init(self, dict, geom)
    class(generalSource), intent(inout)      :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    character(30)                            :: type
    character(:), allocatable                :: matName
    integer(shortInt)                        :: matID, uniqueID, cellID, MT, iso
    logical(defBool)                         :: isCE, isMG, hasDir, hasPos,&
                                                hasMat, hasCell, hasMT, hasIso
    character(100)                           :: dataType
    character(100), parameter :: Here = 'init (generalSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    ! Establish source particle type distribution

    ! Identify which particle is used in the source
    ! Presently limited to neutron and photon
    call dict % getOrDefault(type, 'particle' ,'neutron')
    select case(type)
      case('neutron')
        self % particleType = P_NEUTRON

      case('photon')
        self % particleType = P_PHOTON

      case default
        call fatalError(Here, 'Unrecognised particle type')
     
    end select

    ! Establish source spatial distribution

    ! Check if point source and check it's inside geometry
    hasPos = dict % isPresent('r')
    if (hasPos) then
      call dict % get(self % r, 'r') 
      if (size(self % r) /= 3) then
        call fatalError(Here, 'Source position must have three components')
      end if

      call self % geom % whatIsAt(self % r, matID, uniqueID)
      if (matID == OUTSIDE_MAT) then
        call fatalError(Here, 'Source has been placed outside geometry')
      end if
    end if

    ! Check if source is distributed in a material
    hasMat = dict % isPresent('mat')
    if (hasPos .and. hasMat) then
      call fatalError(Here, 'Source cannot be both a point source and material source')
    end if

    ! Identify material idx and check ensure it is present in the geometry
    if (hasMat) then

      call dict % get(matName, 'mat')
      self % mat = matIdx(matName)

      if (self % mat == NOT_FOUND) then
        call fatalError(Here, 'Specified material not found')
      end if

    end if

    ! Check if source is distributed in a cell
    hasCell = dict % isPresent('cell')
    if (hasCell .and. hasPos) then
      call fatalError(Here, 'Source cannot be both a point source and cell source')
    end if
    if (hasCell .and. hasMat) then
      call fatalError(Here, 'Source cannot be both a material source and cell source')
    end if

    ! Note: can't check that cellID exists in the problem
    !       need to wait for sampling to fail to establish this!
    if (hasCell) then
      call dict % get(cellID, 'cell')
      self % cell = cellID
    end if

    ! Check for bounding box to improve sampling efficiency when using
    ! cell or material sampling
    ! Otherwise, use bb from geometry
    if (hasCell .or. hasMat) then
      hasBb = dict % isPresent('bb')

      ! Take user-provided bounding box
      if (hasBb) then
        call dict % get(self % bb, 'bb')

        if (size(self % bb /= 6)) then
          call fatalError(Here, 'Bounding box must have 6 entries: [x1 x2 y1 y2 z1 z2]')
        end if

      ! Use geometry's bounding box
      else
        allocate(self % bb(6))
        self % bb = self % geom % bounds()

      end if

    end if

    ! Establish source direction distribution

    ! Get beam direction and normalise - otherwise, assume isotropic
    hasDir = dict % isPresent('dir')
    if (hasDir) then
      
      call dict % get(self % dir, 'dir')
      if (size(self % dir) /= 3) then
        call fatalError(Here, 'Source direction must have three components')
      end if

      self % dir = self % dir / norm2(self % dir)
      self % isIsotropic = .false.

    else
      self % isIsotropic = .true.
    end if

    ! Establish source energy distribution

    hasMT = dict % isPresent('mt')
    if hasMT then

      ! Must inform of energy type - assume CE otherwise
      call dict % getOrDefault(dataType,'data','ce')
      select case(dataType)

        case('ce')
          self % isMG = .false.

        case('mg')
          self % isMG = .true.
        
        case default
          call fatalError(Here, 'Invalid source data type specified: must be ce or mg')

      end select

      ! May need to add some extra reactions here - only included the obvious ones!
      call dict % get(MT, 'mt')
      select case(MT)

        ! Microscopic fission
        case(18)

          if (self % particleType /= P_NEUTRON) then
            call fatalError(Here, 'MT 18 source has only been implemented for neutrons')
          end if

          if (self % isMG) then
            call fatalError(Here, 'Cannot sample a microscopic reaction rate with MG data')
          else

            hasIso = dict % isPresent('iso')
            if hasIso then
              call dict % get(iso,'iso')

            else
              call fatalError(Here, &
              'Must provide an isotope (e.g., 92235) for fixed source '//&
              'microscopic reaction sampling')
            end if
     
   !         neutronData => ceNeutronDatabase_CptrCast(self % nucData)
   !         call neutronData % energyBounds(E_down, E_up)
            
   !         ! Get reaction object
   !         self % reacCE => fissionCE_TptrCast(self % nucData % getReaction(N_FISSION, nucIdx))
   !         if(.not.associated(self % reacCE)) call fatalError(Here, "Failed to get CE Fission Reaction Object")

          end if

        ! Macroscopic reactions
        case(:-1)

          if (self % particleType /= P_NEUTRON) then
            call fatalError(Here, 'Macroscopic reaction source has only been implemented for neutrons')
          end if



        case default
          call fatalError(Here, 'The provided MT number is unsupported')
     
      end select

    else

      ! Source is mono-energetic
      ! Get particle energy/group
      isCE = dict % isPresent('E')
      isMG = dict % isPresent('G')
      if (isCE .and. isMG) then
        call fatalError(Here, 'Source may be either continuous energy or MG, not both')
      elseif (isCE) then
        call dict % get(self % E, 'E')
      elseif (isMG)
        call dict % get(self % G, 'G')
        self % isMG = .true.
      else
        call fatalError(Here, 'Must specify source energy, either CE or MG')
      end if

    end if
 
  end subroutine init

  !!
  !! Provide particle type
  !!
  !! Particle type is fixed on initialisation, this routine simply passes the particle type
  !!
  !! Inputs:
  !!   p [inout] -> particle to be given a type
  !!   rand [in] -> pointer to random number generator
  !!
  !! Result:
  !!   Particle is provided with a type
  !!
  subroutine sampleType(self, p, rand)
    class(generalSource), intent(inout) :: self
    class(particleState), intent(inout) :: p
    class(RNG), pointer, intent(in)     :: rand 
    
    p % type = self % particleType

  end subroutine sampleType

  !!
  !! Provide particle position
  !!
  !! Particle location is either taken from the fixed value specified on
  !! initialisation, or sampled from a bounding box until either the specified
  !! material or cell of the source is located
  !!
  !! Inputs:
  !!   p [inout] -> particle to be given a position
  !!   rand [in] -> pointer to random number generator
  !!
  !! Result:
  !!   Particle is provided with a position
  !!
  !! Error:
  !!   Will return an error when sampling position if material/cell isn't
  !!   located after maxSample attempts. This can occur either if the bounding
  !!   box is very large compared to the source material/cell, or because the cell
  !!   ID to be sampled from doesn't exist in the geometry (this is tricky to check
  !!   on initialisation)
  !!   In the former case, this can be remedied by providing a better bounding box
  !!   In the latter case, a correct cell uniqueID should be specified in the input
  !!
  subroutine samplePosition(self, p, rand)
    class(generalSource), intent(inout) :: self
    class(particleState), intent(inout) :: p
    class(RNG), pointer, intent(in)     :: rand
    real(defReal), dimension(3)         :: r
    real(defReal)                       :: dx, dy, dz, r1, r2, r3
    integer(shortInt)                   :: i, matIdx, cellID
    integer(shortInt), parameter        :: maxSample = 10000
    character(100), parameter :: Here = 'samplePosition (generalSource_class.f90)'

    i = 0
    if ((self % mat > OUTSIDE) .or. (self % cell > 0)) then
      
      dx = self % bb(2) - self % bb(1)
      dy = self % bb(4) - self % bb(3)
      dz = self % bb(6) - self % bb(5)

      ! Sample until a valid point is found or until
      ! maxSample is exceeded
      do while i < maxSample

        r1 = rand % get()
        r2 = rand % get()
        r3 = rand % get()

        r = [dx * r1, dy * r2, dz * r3]

        call self % geom % whatIsAt(r, matIdx, cellID)
        
        ! Exit loop if the correct material/cell is located
        ! Provide matIdx to accelerate macroscopic reaction sampling
        if ((self % mat == matIdx) .or. (self % cell == cellID)) then
          p % matIdx = matIdx
          exit
        end if

        i = i + 1
      end do

      if (i == maxSample) then
        call fatalError(Here, 'Low sampling efficiency: this is due either to a '//&
          'large bounding box or because a cellID was specified that is not in the geometry'
      end if

      p % r = r

    else
      p % r = self % r
    end if

  end subroutine samplePosition
 
  !!
  !! Provide angle or sample if isotropic
  !!
  !! Particle direction is either fixed or isotropic
  !! This subroutine either passes the direction or samples from the unit sphere
  !!
  !! Inputs:
  !!   p [inout] -> particle to be given a direction
  !!   rand [in] -> pointer to random number generator
  !!
  !! Result:
  !!   Particle is provided with a direction
  !!
  subroutine sampleEnergyAngle(self, p, rand)
    class(generalSource), intent(inout) :: self
    class(particleState), intent(inout) :: p
    class(RNG), pointer, intent(in)     :: rand 
    real(defReal)                       :: r, phi, theta
    
    if (self % isIsotropic) then
      r = rand % get()
      phi = TWO * PI * r
      r = rand % get()
      theta = acos(1 - TWO * r)
      p % dir = [cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)]
    else
      p % dir = self % dir
    end if

  end subroutine sampleEnergyAngle


  !!
  !! Provide particle energy
  !!
  !! Particle energy is fixed - subroutine simply provides it
  !! Applies to either CE or MG particles
  !!
  !! Inputs:
  !!   p [inout] -> particle to be given an energy
  !!   rand [in] -> pointer to random number generator
  !!
  !! Result:
  !!   Particle is provided with an energy
  !!
  subroutine sampleEnergy(self, p, rand)
    class(generalSource), intent(inout) :: self
    class(particleState), intent(inout) :: p
    class(RNG), pointer, intent(in)     :: rand 
    character(100), parameter :: Here = 'sampleEnergy (generalSource_class.f90)'

    ! Source is a macroscopic reaction
    if (MACROSCOPIC CONDITION) then
            ! MOSTLY LAZY COPY PASTING
      mat => neutronMaterial_CptrCast(self % nucData % getMaterial(matIdx))
      if(.not.associated(mat)) call fatalError(Here, "Nuclear data is not for neutrons")

      ! Fatal error if material isn't fissile
      if( .not.mat % isFissile()) then
        call fatalError(Here, "Macroscopic source material is not fissile")
      end if

      ! Generate new fission site
      select type(mat)
        class is(ceNeutronMaterial)

          neutronData => ceNeutronDatabase_CptrCast(self % nucData)
          call neutronData % energyBounds(E_down, E_up)
          
          ! Select nuclide
          nucIdx = mat % sampleFission(1.0E-6_defReal, rand)

          ! Get reaction object
          fissCE => fissionCE_TptrCast(self % nucData % getReaction(N_FISSION, nucIdx))
          if(.not.associated(fissCE)) call fatalError(Here, "Failed to get CE Fission Reaction Object")

          ! Get mu, phi, E_out
          call fissCE % sampleOut(mu, phi, E_out, 1.0E-6_defReal, rand)

          ! Put Data into particle State
          call p% point(rotateVector([ONE, ZERO, ZERO], mu, phi))
          p % E = E_out
          p % isMG = .false.

          ! Prevent particle energy above upper bound being sampled
          if (p% E > E_up) p % E = E_up

        class is(mgNeutronMaterial)
          ! Get reaction object
          fissMG => fissionMG_TptrCast(self % nucData % getReaction(macroFission, matIdx))
          if(.not.associated(fissMG)) call fatalError(Here, "Failed to get MG Fission Reaction Object")

          ! Get mu, phi, E_out
          call fissMG % sampleOut(mu, phi, G_out, 1, self % pRNG)

          ! Put Data into particle State
          call neutron % teleport(r)


    ! Source is a CE reaction
    if (associated(self % reacCE)) then


    ! Source is an MG reaction
    else if (associated(self % reacMG)) then


    ! Source is monoenergetic
    else
      if (self % isMG) then
        p % G = self % G
        p % isMG = .true.
      else
        p % E = self % E
        p % isMG = .false.
      end if
    end if

  end subroutine sampleEnergy

  !!
  !! Terminate general source
  !!
  !! Cleans up pointers and allocatables
  !!
  subroutine kill(self)
    class(generalSource), intent(inout) :: self

    self % geom => null()
    if allocated(self % bb) deallocate(self % bb)
    if allocated(self % r) deallocate(self % r)
    if allocated(self % dir) deallocate(self % dir)
    if associated(self % reacCE) self % reacCE => null()
    if associated(self % reacMG) self % reacMG => null()

  end subroutine kill

end module generalSource_class
