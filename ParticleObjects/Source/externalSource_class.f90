module externalSource_class

  use numPrecision
  use, intrinsic :: iso_fortran_env, only : iostat_end
  use errors_mod,                    only : fatalError
  use particle_class,                only : particleState, P_NEUTRON
  use source_inter,                  only : source, kill_super => kill
  use geometry_inter,                only : geometry
  use RNG_class,                     only : RNG
  use dictionary_class,              only : dictionary
  use dictParser_func,               only : fileToDict
  use genericProcedures,             only : numToChar, readArray

  use geometry_inter,          only : geometry
  use neutronMaterial_inter,   only : neutronMaterial, neutronMaterial_CptrCast
  use nuclearDataReg_mod,      only : ndReg_getNeutronCE => getNeutronCE, &
                                      ndReg_getNeutronMG => getNeutronMG
  use nuclearDatabase_inter,   only : nuclearDatabase
  use ceNeutronDatabase_inter, only : ceNeutronDatabase
  use mgNeutronDatabase_inter, only : mgNeutronDatabase
  use materialMenu_mod,        only : mm_matIdx => matIdx

  implicit none
  private

  type, public, extends(source) :: externalSource
    private
        integer(shortInt)                            :: numberNeutrons
        real(defReal), dimension(:,:), allocatable   :: r, dir
        real(defReal), dimension(:), allocatable     :: E 
        integer(shortInt), dimension(:), allocatable :: G
        logical(defBool)                             :: isMG = .false.
        logical(defBool)                             :: readBinary = .false. ! Flag to read binary file
    ! Add any data members specific to externalSource here

  contains
    procedure :: init
    procedure :: sampleParticle
    procedure :: kill
  end type externalSource

contains
  !!
  !! Initialise external source from dictionary and list of neutrons
  !! Does not check if source file is consistent with specified input format
  !!
  !! See source_inter for interface details
  !!
  !! Errors if:
  !!   - path to source file not specified
  !!   - source data type inconsistent with nuclear database
  !!   - failed to read source file
  !!
  subroutine init(self, dict, geom)
    class(externalSource), intent(inout)     :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    class(nuclearDatabase), pointer          :: nucData
    character(pathLen)                       :: path
    character(nameLen)                       :: energy   
    integer(shortInt)                        :: i, id
    logical(defBool)                         :: EOF
    real(defReal)                            :: dummy(9)
    character(100),parameter :: Here = 'init (externalSource_class.f90)'

    ! Provide geometry info to source
    self % geom => geom

    ! Select energy type and allocate energy array
    call dict % getOrDefault(energy, 'data', 'ce')
    select case(energy)
      case('ce')
        self % isMG = .false.
      case('mg')
        self % isMG = .true.
      case default
        call fatalError(Here, 'Invalid source data type specified: must be ce or mg')
    end select
    
    ! Get pointer to appropriate nuclear database
    if (self % isMG) then
      nucData => ndReg_getNeutronMG()
    else
      nucData => ndReg_getNeutronCE()
    end if
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

    ! Select path to external source file
    if (dict % isPresent('path')) then
        call dict % get(path, 'path')
    else
        call fatalError(Here, 'path must be specified in the dictionary for externalSource')
    end if

    call dict % getOrDefault(self % readBinary, 'binary', .false.)

    ! Check consistency of energy type with nuclear database
    select type (nucData)
      class is (ceNeutronDatabase)
        if (self % isMG) then
          call fatalError(Here, "Inconsistent external source: ce database with mg data")
        end if
      class is (mgNeutronDatabase)
        if (.not. self % isMG) then
          call fatalError(Here, "Inconsistent external source: mg database with ce data")
        end if
      class default
        call fatalError(Here, "Unrecognised type of nuclearDatabase")
    end select

    ! Read source data from file, binary or ASCII
    self % numberNeutrons = 0
    id = 10
    ! ASCII file reading
    if (.not. self % readBinary) then
      open(unit=id, file=trim(path), status='old', action='read')
      print *, 'Reading external source from ASCII file: ', trim(path)
    else if (self % readBinary) then
      ! Binary file reading
      open(unit=id, file=trim(path), status='old', access='stream', form='unformatted', action='read')
      print *, 'Reading external source from binary file: ', trim(path)
    else 
      call fatalError(Here, 'Invalid external source file reading mode specified')
    end if
      ! Read number of neutrons from start of ASCII file
    do
      call readArray(id, self % readBinary, dummy, EOF)
      if (EOF) exit
      self % numberNeutrons = self % numberNeutrons + 1
    end do
    ! Reset to start of file
    rewind(id)
    ! Allocate position and direction arrays
    allocate(self % r(3, self % numberNeutrons))
    allocate(self % dir(3, self % numberNeutrons))

    ! Allocate energy or group arrays
    if (self % isMG) then 
      allocate(self % G(self % numberNeutrons))
    else
      allocate(self % E(self % numberNeutrons))
    end if

    ! Read and store neutron data from source file
    ! Value for BroodID is ignored
    print *, self % numberNeutrons, ' neutrons read from external source file'
    do i = 1, self % numberNeutrons
      call readArray(id, self % readBinary, dummy, EOF)
      if (EOF) exit
      self % r(:, i) = dummy(1:3)
      self % dir(:, i) = dummy(4:6)
      if (self % isMG) then 
        self % G(i) = int(dummy(8))
      else
        self % E(i) = dummy(7)
      end if
    end do

    close(id)
  end subroutine init

  !!
  !! Sample the particle's phase space co-ordinates from the external source
  !!
  !! See source_inter for details
  !!
  function sampleParticle(self, rand) result(p)
    class(externalSource), intent(inout)     :: self
    class(RNG), intent(inout)                :: rand
    type(particleState)                      :: p
    class(nuclearDatabase), pointer          :: nucData
    class(neutronMaterial), pointer          :: mat
    integer(shortInt)                        :: matIdx, uniqueID, idx
    character(100),parameter :: Here = 'sampleParticle (externalSource_class.f90)'

    ! Get pointer to appropriate nuclear database
    if (self % isMG) then
      nucData => ndReg_getNeutronMG()
    else
      nucData => ndReg_getNeutronCE()
    end if
    if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

    ! Sample index of source neutron to be used
    idx = int(rand % get() * real(self % numberNeutrons, defReal)) + 1

    ! Sample neutron from external source data, check we have not exceeded number of neutrons
    if (idx <= self % numberNeutrons) then
      ! Find material at neutron position
      call self % geom % whatIsAt(matIdx, uniqueID, self % r(:, idx))
      
      ! Get pointer to neutron material, return fatal error if not found
      mat => neutronMaterial_CptrCast(nucData % getMaterial(matIdx))
      if (.not.associated(mat)) call fatalError(Here, "Nuclear data did not return neutron material.")
      
      ! Assign basic phase space coordinates
      p % r = self % r(:, idx)
      p % dir = self % dir(:, idx)
      p % type = P_NEUTRON
      p % time = ZERO
      p % wgt = ONE
      p % uniqueID = uniqueID
      p % matIdx = matIdx

      ! Set energy
      if (self % isMG) then
        p % G = self % G(idx)
        p % isMG = .true.
      else
        p % E = self % E(idx)
        p % isMG = .false.
      end if
    else
        call fatalError(Here, 'Requested neutron is not in the external source file')
    end if
    
  end function sampleParticle

  elemental subroutine kill(self)
    class(externalSource), intent(inout) :: self
    character(100),parameter :: Here = 'kill (externalSource_class.f90)'

    call kill_super(self)

    self % numberNeutrons = 0
    if (allocated(self % r)) deallocate(self % r)
    if (allocated(self % dir)) deallocate(self % dir)
    if (allocated(self % E)) deallocate(self % E)
    if (allocated(self % G)) deallocate(self % G)

  end subroutine kill

end module externalSource_class