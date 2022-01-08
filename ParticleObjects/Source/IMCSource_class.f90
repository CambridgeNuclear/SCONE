module IMCSource_class

  use numPrecision
  use endfConstants
  use universalVariables
  use genericProcedures,       only : fatalError, rotateVector
  use dictionary_class,        only : dictionary
  use RNG_class,               only : RNG

  use particle_class,          only : particleState, P_PHOTON
  use source_inter,            only : source, kill_super => kill

  use geometry_inter,          only : geometry
  use IMCMaterial_inter,       only : IMCMaterial, IMCMaterial_CptrCast
  use nuclearDataReg_mod,      only : ndReg_getIMCMG => getIMCMG
  use nuclearDatabase_inter,   only : nuclearDatabase
  use mgIMCDatabase_inter,     only : mgIMCDatabase

  implicit none
  private

  !!
  !! IMC Source from distributed fission sites
  !!
  !! Places fission sites uniformly in regions with fissile material.
  !! Spectrum of the fission IMC is such as if it fission was caused by incdent
  !! IMC with CE energy E or MG with group G.
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
  !!     type imcSource;
  !!     #data MG; #
  !!     #E 15.0;  #
  !!     #G 7;     #
  !!   }
  !!
  type, public,extends(source) :: imcSource
    private
    logical(defBool)            :: isMG   = .true.
    real(defReal), dimension(3) :: bottom = ZERO
    real(defReal), dimension(3) :: top    = ZERO
    real(defReal)               :: E      = ZERO
    integer(shortInt)           :: G      = 0
  contains
    procedure :: init
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
    character(nameLen)                       :: type
    real(defReal), dimension(6)              :: bounds
    character(100), parameter :: Here = 'init (imcSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    ! Set bounding region
    bounds = self % geom % bounds()
    self % bottom = bounds(1:3)
    self % top    = bounds(4:6)

  end subroutine init

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
    real(defReal)                        :: mu, phi
    integer(shortInt)                    :: matIdx, uniqueID, nucIdx, i
    character(100), parameter :: Here = 'sampleParticle (imcSource_class.f90)'

    ! Get pointer to appropriate nuclear database
    nucData => ndReg_getIMCMG()
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

    i = 0
    rejection : do
      ! Protect against infinite loop
      i = i +1
      if ( i > 200) then
        call fatalError(Here, 'Infinite loop in sampling of fission sites. Please check that&
                              & defined volume contains fissile material.')
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

      mat => IMCMaterial_CptrCast(nucData % getMaterial(matIdx))        ! Currently will only work as intended with 1 cell
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
      p % wgt      = ONE
      p % time     = ZERO
      p % type     = P_PHOTON
      p % r        = r
      p % dir      = dir

      ! Set Energy
      p % E = mat % getEmittedRad() / 5 ! Currently fixed at 5 particles for simplicity
 
      ! Set Time
      p % time = rand % get() * timeStepSize ! + Start of time step time


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
    self % E      = ZERO
    self % G      = 0

  end subroutine kill

end module IMCSource_class
