module IMCSource_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures,       only : fatalError, rotateVector
  use dictionary_class,        only : dictionary
  use RNG_class,               only : RNG

  use particle_class,          only : particle, particleState, P_PHOTON
  use particleDungeon_class,   only : particleDungeon
  use source_inter,            only : source, kill_super => kill

  use geometry_inter,          only : geometry
  use IMCMaterial_inter,       only : IMCMaterial, IMCMaterial_CptrCast
  use nuclearDataReg_mod,      only : ndReg_getIMCMG => getIMCMG
  use nuclearDatabase_inter,   only : nuclearDatabase
  use mgIMCDatabase_inter,     only : mgIMCDatabase
  use materialMenu_mod,        only : MMnMat => nMat

  implicit none
  private

  !!
  !! IMC Source for uniform generation of photons within material regions
  !!
  !! Angular distribution is isotropic.
  !!
  !! Private members:
  !!   isMG    -> is the source multi-group? (default = .true.)
  !!   bottom  -> Bottom corner (x_min, y_min, z_min)
  !!   top     -> Top corner (x_max, y_max, z_max)
  !!   G       -> Group (default = 1)
  !!   matPops -> Array to store the number of particles sampled in each material for
  !!              normalisation of weight
  !!
  !! Interface:
  !!   source_inter Interface
  !!
  !! SAMPLE INPUT:
  !!   imcSource { type IMCSource; }
  !!
  type, public,extends(source) :: imcSource
    private
    logical(defBool)                             :: isMG   = .true.
    real(defReal), dimension(3)                  :: bottom = ZERO
    real(defReal), dimension(3)                  :: top    = ZERO
    integer(shortInt)                            :: G      = 0
    integer(shortInt), dimension(:), allocatable :: matPops
  contains
    procedure :: init
    procedure :: append
    procedure :: sampleParticle
    procedure :: kill
  end type imcSource

contains

  !!
  !! Initialise IMC Source
  !!
  !! See source_inter for details
  !!
  subroutine init(self, dict, geom)
    class(imcSource), intent(inout)          :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    real(defReal), dimension(6)              :: bounds
    integer(shortInt)                        :: i, n
    character(100), parameter :: Here = 'init (imcSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    call dict % getOrDefault(self % G, 'G', 1)

    ! Set bounding region
    bounds = self % geom % bounds()
    self % bottom = bounds(1:3)
    self % top    = bounds(4:6)

    ! Initialise array to store numbers of particles
    n = MMnMat()
    allocate(self % matPops(n))
    do i=1, n
      self % matPops(i) = 0
    end do

  end subroutine init

  !!
  !! Generate n particles to add to a particleDungeon without overriding
  !! particles already present. More complex than superclass 'append' subroutine,
  !! needed for multiregion functionality.
  !!
  !! The number of particles sampled in each matIdx is recorded and used to normalise
  !! each particle weight, so that the total energy emitted in each region is as
  !! required
  !!
  !! Args:
  !!   dungeon [inout] -> particle dungeon to be added to
  !!   n [in]          -> number of particles to place in dungeon
  !!   rand [inout]    -> particle RNG object
  !!
  !! Result:
  !!   A dungeon populated with n particles sampled from the source, plus particles
  !!   already present in dungeon
  !!
  subroutine append(self, dungeon, N, rand)
    class(imcSource), intent(inout)      :: self
    type(particleDungeon), intent(inout) :: dungeon
    integer(shortInt), intent(in)        :: N
    class(RNG), intent(inout)            :: rand
    type(particleDungeon)                :: tempDungeon
    type(particle)                       :: p
    integer(shortInt)                    :: i
    real(defReal)                        :: normFactor
    character(100), parameter            :: Here = "append (IMCSource_class.f90)"

    ! Reset particle population counters
    do i = 1, size( self % matPops )
      self % matPops(i) = 0
    end do

    ! Set temporary dungeon size
    call tempDungeon % setSize(n)

    ! Generate n particles to populate temporary dungeon
    do i = 1, n
      call tempDungeon % replace(self % sampleParticle(rand), i)
    end do

    ! Call error if any region contains no generated particles (due to small regions and/or
    !   not enough particles used), needed for now as otherwise will lead to energy imbalance
    !   as mat energy will be reduced by emittedRad but no particles will be carrying it
    ! Note that matPops is set to 1 in sample_particle if region is of 0 temperature to avoid
    !   this error for such a case
    if ( minval(self % matPops) == 0 ) then
      call fatalError(Here, "Not all regions emitted particles, use more particles")
    end if

    ! Loop through again and add to input dungeon, normalising energies based on material
    do i = 1, n

      call tempDungeon % release(p)

      ! Place inside geometry to set matIdx, for some reason resets when released from dungeon
      call self % geom % placeCoord( p % coords )

      ! Normalise
      normFactor = self % matPops( p % coords % matIdx )
      p % w = p % w / normFactor

      ! Add to input dungeon
      call dungeon % detain(p)

    end do
   
  end subroutine append

  !!
  !! Sample particle's phase space co-ordinates
  !!
  !! See source_inter for details
  !!
  function sampleParticle(self, rand) result(p)
    class(imcSource), intent(inout)      :: self
    class(RNG), intent(inout)            :: rand
    type(particleState)                  :: p
    class(nuclearDatabase), pointer      :: nucData
    class(IMCMaterial), pointer          :: mat
    real(defReal), dimension(3)          :: r, rand3, dir
    ! Here, i is a float to allow more precise control of loop
    real(defReal)                        :: mu, phi, i
    integer(shortInt)                    :: matIdx, uniqueID
    character(100), parameter :: Here = 'sampleParticle (imcSource_class.f90)'

    ! Get pointer to appropriate nuclear database
    nucData => ndReg_getIMCMG()
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

    ! Position is sampled by taking a random point from within geometry bounding box
    ! If in valid material, position is accepted
    i = 0
    rejection : do
      ! Protect against infinite loop
      i = i + 1
      if ( i > 200) then
        call fatalError(Here, '200 particles in a row sampled in void or outside material.&
                             & Check that geometry is as intended')
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

      ! Point to material
      mat => IMCMaterial_CptrCast(nucData % getMaterial(matIdx))
      if (.not.associated(mat)) call fatalError(Here, "Nuclear data did not return IMC material.")

      ! Sample Direction - chosen uniformly inside unit sphere
      mu = 2 * rand % get() - 1
      phi = rand % get() * 2*pi
      dir(1) = mu
      dir(2) = sqrt(1-mu**2) * cos(phi)
      dir(3) = sqrt(1-mu**2) * sin(phi)

      ! Assign basic phase-space coordinates
      p % matIdx   = matIdx
      p % uniqueID = uniqueID
      p % time     = ZERO
      p % type     = P_PHOTON
      p % r        = r
      p % dir      = dir
      p % G        = self % G
      p % isMG     = .true.

      ! Set weight to be equal to total emitted radiation from material
      ! This weight is then normalised later - see appendIMC (source_inter.f90)
      ! There may be more intuitive ways to do this, but works well for now
      p % wgt = mat % getEmittedRad()

      ! Don't sample particles from areas of 0 temperature
      if( p % wgt == 0 ) then
        self % matPops(matIdx) = 1   ! Set to 1 to avoid error in append subroutine
        i = i - 0.9                  ! To allow more attempts if large regions with 0 temp
        cycle rejection
      end if

      ! Increase counter of number of particles in material in order to normalise later
      self % matPops(matIdx) = self % matPops(matIdx) + 1

      ! Exit the loop
      exit rejection

    end do rejection

  end function sampleParticle

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(imcSource), intent(inout) :: self

    call kill_super(self)

    self % isMG   = .true.
    self % bottom = ZERO
    self % top    = ZERO
    self % G      = 0

  end subroutine kill

end module IMCSource_class
