module materialSource_class

  use numPrecision
  use universalVariables,      only : OUTSIDE_MAT, VOID_MAT
  use genericProcedures,       only : fatalError, rotateVector
  use dictionary_class,        only : dictionary
  use RNG_class,               only : RNG
  use charMap_class,           only : charMap

  use particle_class,          only : particleState, P_NEUTRON
  use source_inter,            only : source, kill_super => kill

  use geometry_inter,          only : geometry
  use neutronMaterial_inter,   only : neutronMaterial, neutronMaterial_CptrCast
  use nuclearDataReg_mod,      only : ndReg_getNeutronCE => getNeutronCE, &
                                      ndReg_getNeutronMG => getNeutronMG
  use nuclearDatabase_inter,   only : nuclearDatabase
  use ceNeutronDatabase_inter, only : ceNeutronDatabase
  use mgNeutronDatabase_inter, only : mgNeutronDatabase

  implicit none
  private

  !!
  !! Neutron Source from a given material
  !!
  !! Places neutrons uniformly in regions with given material.
  !! Angular distribution is isotropic.
  !!
  !! Can be fed a bounding box to increase sampling efficiency
  !!
  !! Private members:
  !!   isMG   -> is the source multi-group? (default = .false.)
  !!   bottom -> Bottom corner (x_min, y_min, z_min)
  !!   top    -> Top corner (x_max, y_max, z_max)
  !!   E      -> Source site energy [MeV] (default = 1.0E-6)
  !!   G      -> Source site group (default = 1)
  !!
  !! Interface:
  !!   source_inter Interface
  !!
  !! Sample Dictionary Input:
  !!   matsource {
  !!     type materialSource;
  !!     #data mg; #
  !!     #E 15.0;  #
  !!     #G 7;     #
  !!     #bb (-x -y -z +x +y +z); #
  !!   }
  !!
  type, public,extends(source) :: materialSource
    private
    logical(defBool)            :: isMG   = .false.
    real(defReal), dimension(3) :: bottom = ZERO
    real(defReal), dimension(3) :: top    = ZERO
    real(defReal)               :: E      = ZERO
    integer(shortInt)           :: G      = 0
    integer(shortInt)           :: matIdx = -1
  contains
    procedure :: init
    procedure :: sampleParticle
    procedure :: kill
  end type materialSource

contains

  !!
  !! Initialise Material Source
  !!
  !! See source_inter for details
  !!
  subroutine init(self, dict, geom)
    class(materialSource), intent(inout)     :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    character(nameLen)                       :: type
    character(nameLen)                       :: matName 
    class(nuclearDatabase), pointer          :: nucData
    class(charMap), pointer                  :: matMap
    real(defReal), dimension(6)              :: bounds
    real(defReal), dimension(:), allocatable :: tempArray
    character(100), parameter :: Here = 'init (materialSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    ! Select energy type
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

    ! Determine which material to sample source in
    call dict % get(matName, 'mat')
    if (self % isMG) then
      nucData => ndReg_getNeutronMG()
    else
      nucData => ndReg_getNeutronCE()
    end if
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')
    matMap => nucData % matNamesMap()
    self % matIdx = matMap % get(matName)

    ! Set bounding region
    if (dict % isPresent('bb')) then
      call dict % get(tempArray, 'bb')
      if (size(tempArray) /= 6) call fatalError(Here,&
               'Bounding box must have 6 entries')
      self % bottom = tempArray(1:3)
      self % top    = tempArray(4:6)

    else
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
    class(materialSource), intent(inout) :: self
    class(RNG), intent(inout)            :: rand
    type(particleState)                  :: p
    class(nuclearDatabase), pointer      :: nucData
    class(neutronMaterial), pointer      :: mat
    real(defReal), dimension(3)          :: r, rand3
    real(defReal)                        :: mu, phi
    integer(shortInt)                    :: matIdx, uniqueID, i
    character(100), parameter :: Here = 'sampleParticle (materialSource_class.f90)'

    ! Get pointer to approperiate nuclear database
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
      if ( i > 200) then
        call fatalError(Here, 'Infinite loop in sampling source. Please check that'//&
                              ' defined volume contains source material.')
      end if

      ! Sample Position
      rand3(1) = rand % get()
      rand3(2) = rand % get()
      rand3(3) = rand % get()
      r = (self % top - self % bottom) * rand3 + self % bottom

      ! Find material under position
      call self % geom % whatIsAt(matIdx, uniqueID, r)

      ! Reject if there is no material
      if (matIdx == OUTSIDE_MAT) cycle rejection

      mat => neutronMaterial_CptrCast(nucData % getMaterial(matIdx))
      if (.not.associated(mat)) call fatalError(Here, "Nuclear data did not return neutron material.")

      ! Resample position if material is not the specified material
      if (.not. (matIdx == self % matIdx)) cycle rejection

      ! Assign basic phase-space coordinates
      p % matIdx   = matIdx
      p % uniqueID = uniqueID
      p % wgt      = ONE
      p % time     = ZERO
      p % type     = P_NEUTRON
      p % r        = r

      mu = TWO * rand % get() - ONE
      phi = TWO_PI * rand % get()

      ! Set Energy
      select type (nucData)
        class is (ceNeutronDatabase)

          p % E = self % E
          p % isMG = .false.
          p % dir  = rotateVector([ONE, ZERO, ZERO], mu, phi)

        class is (mgNeutronDatabase)

          p % G = self % G
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
    class(materialSource), intent(inout) :: self

    call kill_super(self)

    self % isMG   = .false.
    self % bottom = ZERO
    self % top    = ZERO
    self % E      = ZERO
    self % G      = 0

  end subroutine kill

end module materialSource_class
