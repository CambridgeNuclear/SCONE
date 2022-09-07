module ISMCSource_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures,       only : fatalError, rotateVector
  use dictionary_class,        only : dictionary
  use RNG_class,               only : RNG

  use particle_class,          only : particleState, P_MATERIAL
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
  !! ISMC Source for uniform generation of photons within material regions
  !!
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
  type, public,extends(source) :: ismcSource
    private
    logical(defBool)            :: isMG   = .true.
    real(defReal), dimension(3) :: bottom = ZERO
    real(defReal), dimension(3) :: top    = ZERO
    real(defReal)               :: E      = ZERO
    integer(shortInt)           :: G      = 0
    integer(shortInt)           :: N      = 10
    real(defReal)               :: boundingVol = ZERO
  contains
    procedure :: init
    procedure :: sampleParticle
    procedure :: kill
  end type ismcSource

contains

  !!
  !! Initialise IMC Source
  !!
  !! See source_inter for details
  !!
  subroutine init(self, dict, geom)
    class(ismcSource), intent(inout)          :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    character(nameLen)                       :: type
    real(defReal), dimension(6)              :: bounds
    real(defReal), dimension(3)              :: boundSize
    integer(shortInt)                        :: i, n
    character(100), parameter :: Here = 'init (ismcSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    call dict % getOrDefault(self % G, 'G', 1)
    call dict % getOrDefault(self % N, 'N', 10)

    ! Set bounding region
    bounds = self % geom % bounds()
    self % bottom = bounds(1:3)
    self % top    = bounds(4:6)

    ! Calculate volume of bounding region
    boundSize = self % top - self % bottom
    self % boundingVol = boundSize(1) * boundSize(2) * boundSize(3)

!    ! Initialise array to store numbers of particles
!    n = MMnMat()
!    allocate( self % matPops(n) )
!    do i=1, n
!      self % matPops(i) = 0
!    end do

  end subroutine init

  !!
  !! Sample particle's phase space co-ordinates
  !!
  !! See source_inter for details
  !!
  function sampleParticle(self, rand) result(p)
    class(ismcSource), intent(inout)      :: self
    class(RNG), intent(inout)            :: rand
    type(particleState)                  :: p
    class(nuclearDatabase), pointer      :: nucData
    class(IMCMaterial), pointer          :: mat
    real(defReal), dimension(3)          :: r, rand3, dir
    ! Here, i is a float to allow more precise control of loop
    real(defReal)                        :: mu, phi, i
    integer(shortInt)                    :: matIdx, uniqueID, nucIdx
    character(100), parameter :: Here = 'sampleParticle (ismcSource_class.f90)'

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
      p % type     = P_MATERIAL
      p % r        = r
      p % dir      = dir
      p % G        = self % G
      p % isMG     = .true.

      p % wgt = mat % getEnergyDens() * self % boundingVol / self % N

!      ! Don't sample particles from areas of 0 temperature
!      if( p % wgt == 0 ) then
!        self % matPops(matIdx) = 1   ! Set to 1 to avoid error in appendIMC (source_inter.f90)
!        i = i - 0.9                  ! To allow more attempts if large regions with 0 temp
!        cycle rejection
!      end if

!      ! Increase counter of number of particles in material in order to normalise later
!      self % matPops(matIdx) = self % matPops(matIdx) + 1

      ! Exit the loop
      exit rejection

    end do rejection

  end function sampleParticle

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(ismcSource), intent(inout) :: self

    !call kill_super(self)

    self % isMG   = .true.
    self % bottom = ZERO
    self % top    = ZERO
    self % E      = ZERO
    self % G      = 0
    self % N      = 10

  end subroutine kill

end module ISMCSource_class
