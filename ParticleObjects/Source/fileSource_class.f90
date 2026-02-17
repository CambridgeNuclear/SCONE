module fileSource_class

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
  use universalVariables,            only : OUTSIDE_MAT, UNDEF_MAT
  use display_func,                  only : statusMsg

  implicit none
  private

  type, public, extends(source) :: fileSource
    private
        integer(longInt)                             :: numberNeutrons
        real(defReal), dimension(:,:), allocatable   :: r, dir
        real(defReal), dimension(:), allocatable     :: E, w 
        integer(shortInt), dimension(:), allocatable :: G
        logical(defBool)                             :: isMG = .false.
        logical(defBool)                             :: readBinary = .false. ! Flag to read binary file
    ! Add any data members specific to fileSource here

  contains
    procedure :: init
    procedure :: sampleParticle
    procedure :: kill
  end type fileSource

contains
  !!
  !! Initialise file source from dictionary and list of neutrons
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
    class(fileSource), intent(inout)     :: self
    class(dictionary), intent(in)        :: dict
    class(geometry), pointer, intent(in) :: geom
    character(pathLen)                   :: path
    character(nameLen)                   :: energy   
    integer(shortInt)                    :: i, id
    logical(defBool)                     :: EOF
    real(defReal)                        :: dummy(10)
    character(100), parameter :: Here = 'init (fileSource_class.f90)'

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

    ! Select path to file  source file
    if (dict % isPresent('path')) then
        call dict % get(path, 'path')
    else
        call fatalError(Here, 'path must be specified in the dictionary for fileSource')
    end if

    call dict % getOrDefault(self % readBinary, 'binary', .false.)

    ! Read source data from file, binary or ASCII
    self % numberNeutrons = 0
    id = 10
    
    if (.not. self % readBinary) then

      ! ASCII file reading
      open(unit=id, file=trim(path), status='old', action='read')
      call statusMsg('Reading file source from ASCII file: '//trim(path))

    else if (self % readBinary) then

      ! Binary file reading
      open(unit=id, file=trim(path), status='old', access='stream', form='unformatted', action='read')
      call statusMsg('Reading file source from binary file: '//trim(path))

    else 
      call fatalError(Here, 'Invalid file source file reading mode specified')
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
    allocate(self % r(3, self % numberNeutrons), self % dir(3, self % numberNeutrons))
    
    ! Allocate weight array
    allocate(self % w(self % numberNeutrons))

    ! Allocate energy or group arrays
    if (self % isMG) then 
      allocate(self % G(self % numberNeutrons))
    else
      allocate(self % E(self % numberNeutrons))
    end if

    ! Read and store neutron data from source file
    ! Value for BroodID is ignored
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
      self % w(i) = dummy(10)

    end do

    close(id)

  end subroutine init

  !!
  !! Sample the particle's phase space co-ordinates from the file source
  !!
  !! See source_inter for details
  !!
  function sampleParticle(self, rand) result(p)
    class(fileSource), intent(inout) :: self
    class(RNG), intent(inout)        :: rand
    type(particleState)              :: p
    integer(shortInt)                :: matIdx, uniqueID, idx
    character(100), parameter :: Here = 'sampleParticle (fileSource_class.f90)'

    ! Sample index of source neutron to be used
    idx = int(rand % get() * real(self % numberNeutrons, defReal)) + 1

    ! Sample neutron from file source data, check we have not exceeded number of neutrons
    if (idx <= self % numberNeutrons) then
      ! Find material at neutron position
      call self % geom % whatIsAt(matIdx, uniqueID, self % r(:, idx))
      
      ! Check neutron is outside of geometry or in undefined region
      if (matIdx == OUTSIDE_MAT .or. matIdx == UNDEF_MAT) then
          call fatalError(Here, 'Neutron sampled from file source is outside of geometry or in undefined region.')
      endif
      
      ! Assign basic phase space coordinates
      p % r   = self % r(:, idx)
      p % dir = self % dir(:, idx)
      p % type = P_NEUTRON
      p % time = ZERO
      p % wgt  = self % w(idx)
      p % uniqueID = uniqueID
      p % matIdx   = matIdx

      ! Set energy
      if (self % isMG) then
        p % G = self % G(idx)
        p % isMG = .true.
      else
        p % E = self % E(idx)
        p % isMG = .false.
      end if
    else
        call fatalError(Here, 'Requested neutron is not in the file source file')
    end if
    
  end function sampleParticle

  elemental subroutine kill(self)
    class(fileSource), intent(inout) :: self
    character(100), parameter :: Here = 'kill (fileSource_class.f90)'

    call kill_super(self)

    self % numberNeutrons = 0
    if (allocated(self % r)) deallocate(self % r)
    if (allocated(self % dir)) deallocate(self % dir)
    if (allocated(self % E)) deallocate(self % E)
    if (allocated(self % G)) deallocate(self % G)
    if (allocated(self % w)) deallocate(self % w)

  end subroutine kill

end module fileSource_class
