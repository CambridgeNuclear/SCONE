!!
!! Material Menu is a module (Singleton) that contains global definitions of different materials
!!
!! It exists to make it easier for all databases to refer to the same materials by the same
!! name and index. This is necessary to avoid confusion resulting from different materials with the
!! same name or index in different databases.
!!
!! Public Members:
!!   materialDefs -> array of material definitions of type materialItem
!!   nameMap      -> Map that maps material name to matIdx
!!   colourMap    -> Map that maps matIdx to 24bit colour (to use for visualisation)
!!
!! Interface:
!!   init      -> Load material definitions from a dictionary
!!   kill      -> Return to uninitialised state
!!   display   -> Display information about all defined materials to console
!!   getMatPtr -> Return pointer to a detailed material information (materialItem)
!!   nMat      -> Return number of materials
!!   matName   -> Return material Name given Index
!!   matIdx    -> Return material Index given Name
!!
module materialMenu_mod

  use numPrecision
  use universalVariables, only : NOT_FOUND, VOID_MAT, OUTSIDE_MAT, UNDEF_MAT
  use genericProcedures,  only : fatalError, charToInt, numToChar
  use colours_func,       only : rgb24bit
  use intMap_class,       only : intMap
  use charMap_class,      only : charMap
  use dictionary_class,   only : dictionary

  implicit none
  private

  !!
  !! Information about a single nuclide
  !!
  !! Based somewhat on MCNP conventions.
  !! Atomic and Mass number identify clearly a nuclide species
  !! Evaluation number T allows to refer to multiple states/evaluations of the same nuclide species
  !! E.G. at a Different temperature as in MCNP Library.
  !!
  !! Public members:
  !!   Z -> Atomic number
  !!   A -> Mass number
  !!   T -> Evaluation number
  !!   hasSab -> Does the nuclide have S(a,b) data?
  !!   sabMix -> Does the nuclide mix S(a,b) data?
  !!   file_Sab1 -> First (and maybe only) S(a,b) file
  !!   file_Sab2 -> Second S(a,b) file
  !!
  !! Interface:
  !!   init -> build from a string
  !!
  type, public :: nuclideInfo
    integer(shortInt)  :: Z = -1
    integer(shortInt)  :: A = -1
    integer(shortInt)  :: T = -1
    logical(defBool)   :: hasSab = .false.
    logical(defBool)   :: sabMix = .false.
    character(nameLen) :: file_Sab1
    character(nameLen) :: file_Sab2
  contains
    procedure :: init   => init_nuclideInfo
    procedure :: toChar => toChar_nuclideInfo
  end type nuclideInfo


  !!
  !! This is a type which collects all information about a single material definition
  !!
  !! Public Members:
  !!   name      -> name of material
  !!   matIdx    -> material index of the material
  !!   T         -> material temperature [K]
  !!   dens      -> vector of densities [1/barn/cm]
  !!   nuclides  -> associated vector of nuclide types
  !!   extraInfo -> dictionary with extra keywords
  !!
  !! Interface:
  !!   init -> build material item from dictionary
  !!   kill -> return to uninitialised state
  !!
  !! Sample Input Dictionary:
  !!
  !!   matDef {
  !!     temp 273;
  !!     #moder {1001.03 (h-h2o.43);}#
  !!     #tms 1;#
  !!     composition {
  !!       1001.03  5.028E-02;
  !!       8016.03  2.505E-02;
  !!       5010.03  2.0E-005;
  !!     }
  !!     xsFile /home/uberMoffTarkin/XS/mat1.xs;
  !!     #rgb (255 0 0); # // RGB colour to be used in visualisation
  !!   }
  !!
  !! Sample with stochastic mixing:
  !!   matDef {
  !!     temp 300;
  !!     moder {1001.03 (h-h2o.43 h-h2o.53);}
  !!     composition {
  !!       1001.03  5.028E-02;
  !!       8016.03  2.505E-02;
  !!       5010.03  2.0E-005;
  !!     }
  !!   }
  !!
  !! NOTE: the moder dictionary is optional, necessary only if S(a,b) thermal scattering
  !!       data are used. If some nuclides are included in moder but not in composition,
  !!       an error is raised.
  !!       Including two entries in moder will invoke stochastic mixing, i.e.,
  !!       stochastic interpolation between the two data libraries.
  !!
  type, public :: materialItem
    character(nameLen)                         :: name   = ''
    integer(shortInt)                          :: matIdx = 0
    real(defReal)                              :: T      = ZERO
    real(defReal),dimension(:),allocatable     :: dens
    type(nuclideInfo),dimension(:),allocatable :: nuclides
    type(dictionary)                           :: extraInfo
    logical(defBool)                           :: hasTMS = .false.
  contains
    procedure :: init    => init_materialItem
    procedure :: kill    => kill_materialItem
    procedure :: display => display_materialItem
  end type materialItem

  !! Parameters
  integer(shortInt), parameter :: COL_OUTSIDE = int(z'ffffff', shortInt)
  integer(shortInt), parameter :: COL_VOID    = int(z'000000', shortInt)
  integer(shortInt), parameter :: COL_UNDEF   = int(z'00ff00', shortInt)


  !! MODULE COMPONENTS
  type(materialItem),dimension(:),allocatable,target,public :: materialDefs
  type(charMap), target, public                             :: nameMap
  type(intMap), public                                      :: colourMap

  public :: init
  public :: kill
  public :: display
  public :: getMatPtr
  public :: nMat
  public :: matName
  public :: matIdx

contains

  !!
  !! Initialises materialMenu from a dictionary with material definitions
  !!
  !! Args:
  !!   dict [in] -> dictionary with material definitions
  !!
  !! Errors:
  !!   None from Here
  !!
  subroutine init(dict)
    class(dictionary),intent(in)                :: dict
    character(nameLen),dimension(:),allocatable :: matNames
    integer(shortInt)                           :: i
    character(nameLen)                          :: temp

    ! Clean whatever may be already present
    call kill()

    ! Load all material names
    call dict % keys(matNames,'dict')

    ! Allocate space
    allocate(materialDefs(size(matNames)))

    ! Load definitions
    do i=1,size(matNames)
      call materialDefs(i) % init(matNames(i), i, dict % getDictPtr(matNames(i)))
      call nameMap % add(matNames(i), i)
    end do

    ! Add special Material keywords to the dictionary
    temp = 'void'
    call nameMap % add(temp, VOID_MAT)
    temp = 'outside'
    call nameMap % add(temp, OUTSIDE_MAT)

    !! Load colours for the special materials
    call colourMap % add(VOID_MAT, COL_VOID)
    call colourMap % add(OUTSIDE_MAT, COL_OUTSIDE)
    call colourMap % add(UNDEF_MAT, COL_UNDEF)

  end subroutine init


  !!
  !! Returns material Menu to an uninitialised state
  !!
  subroutine kill()
    integer(shortInt) :: i

    call nameMap % kill()
    if(allocated(materialDefs)) then
      do i=1,size(materialDefs)
        call materialDefs(i) % kill()
      end do
      deallocate(materialDefs)
    end if

  end subroutine kill

  !!
  !! Print material definition information to the console
  !!
  !! Args:
  !!   None
  !! Errors:
  !!   None
  !!
  subroutine display()
    integer(shortInt) :: i

    print '(A60)', repeat('<>',30)
    print '(A)', "^^ MATERIAL DEFINITIONS ^^"

    do i = 1,size(materialDefs)
      call materialDefs(i) % display()
      ! Print separation line
      print '(A)', " ><((((*>  +  <*))))><"
    end do

    print '(A60)', repeat('<>',30)

  end subroutine display

  !!
  !! Return Material Name given index
  !!
  !! Args:
  !!   idx [in] -> Material Index
  !!
  !! Result:
  !!   nameLen long character with material name
  !!
  !! Error:
  !!   If idx is -ve or larger then number of defined materials
  !!   Empty string '' is returned as its name
  !!
  function matName(idx) result(name)
    integer(shortInt), intent(in) :: idx
    character(nameLen)            :: name

    if( idx <= 0 .or. nMat() < idx) then
      name = ''

    else
      name = materialDefs(idx) % name
    end if

  end function matName

  !!
  !! Return material index Given Name
  !!
  !! Args:
  !!   name [in] -> material name
  !!
  !! Result:
  !!   matIdx corresponding to name
  !!
  !! Error:
  !!   If name does not correspond to any defined material NOT_FOUND is returned
  !!
  function matIdx(name) result(idx)
    character(*), intent(in) :: name
    integer(shortInt)        :: idx

    idx = nameMap % getOrDefault(name, NOT_FOUND)

  end function matIdx

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! TYPE PROCEDURES
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Initialise material definition from a dictionary and name
  !!
  !! Args:
  !!   name [in] -> character with material name
  !!   idx  [in] -> material index
  !!   dict [in] -> dictionary with material definition
  !!
  !! Errors:
  !!   FatalError if dictionary does not contain valid material definition.
  !!
  subroutine init_materialItem(self, name, idx, dict)
    class(materialItem), intent(inout)            :: self
    character(nameLen), intent(in)                :: name
    integer(shortInt), intent(in)                 :: idx
    class(dictionary), intent(in)                 :: dict
    character(nameLen), dimension(:), allocatable :: keys, moderKeys
    integer(shortInt), dimension(:), allocatable  :: temp
    integer(shortInt)                             :: i
    class(dictionary),pointer                     :: compDict, moderDict
    logical(defBool)                              :: hasSab, foundModer
    character(nameLen), dimension(:), allocatable :: filenames
    character(100), parameter :: Here = 'init_materialItem (materialMenu_mod.f90)'

    ! Return to initial state
    call self % kill()

    ! Load easy components properties
    self % name = name
    self % matIdx = idx

    ! Check TMS flag and read temperature
    call dict % getOrDefault(self % hasTMS, 'tms', .false.)

    if (self % hasTMS .and. .not. dict % isPresent('temp')) then
      call fatalError(Here, 'The material temperature must be specified when TMS is on')
    end if

    call dict % getOrDefault(self % T, 'temp', ZERO)
    if (self % T < ZERO) call fatalError(Here, 'The temperature of material '//numToChar(idx)//&
                                                ' is negative: '//numToChar(self % T))

    ! Get composition dictionary and load composition
    compDict => dict % getDictPtr('composition')
    call compDict % keys(keys)

    ! Allocate space for nuclide information
    allocate(self % nuclides(size(keys)))
    allocate(self % dens(size(keys)))

    ! Check if S(a,b) files are specified.
    if (dict % isPresent('moder')) then
      moderDict => dict % getDictPtr('moder')
      call moderDict % keys(moderKeys)
      hasSab = .true.
    else
      hasSab = .false.
    end if

    ! Load definitions
    foundModer = .false.
    do i =1,size(keys)
      ! Check if S(a,b) is on and required for that nuclide
      if (hasSab .and. moderDict % isPresent(keys(i))) then
        self % nuclides(i) % hasSab = .true.
        foundModer = .true.

        ! Check for stochastic mixing - this will depend on the
        ! size of the array of files produce
        call moderDict % get(filenames, keys(i))
        if (size(filenames) == 2) then
          self % nuclides(i) % file_Sab1 = filenames(1)
          self % nuclides(i) % file_Sab2 = filenames(2)
          self % nuclides(i) % sabMix = .true.
        elseif (size(filenames) == 1) then
          self % nuclides(i) % file_Sab1 = filenames(1)
        else
          print *,filenames
          call fatalError(Here,'Unexpectedly long moder contents. Should be 1 or 2 '//&
                  'entries.')
        end if
      end if

      ! Initialise the nuclides
      call compDict % get(self % dens(i), keys(i))
      call self % nuclides(i) % init(keys(i))
    end do
    ! Make sure if a moderator is provided the nuclide is present
    ! in the composition
    if ((.not. foundModer) .and. hasSab) then
      print *,moderKeys
      call fatalError(Here, 'Nuclides requested for S(alpha,beta) are not present in composition')
    end if

    ! Add colour info if present
    if(dict % isPresent('rgb')) then
      call dict % get(temp, 'rgb')

      if (size(temp) /= 3) then
        call fatalError(Here, "'rgb' keyword must have 3 values")
      end if

      call colourMap % add(idx, rgb24bit(temp(1), temp(2), temp(3)))
    end if

    ! Save dictionary
    self % extraInfo = dict

    ! TODO: Remove composition subdictionary from extraInfo
    !       Or rather do not copy it in the first place

  end subroutine init_materialItem

  !!
  !! Return material Item to uninitialised state
  !!
  subroutine kill_materialItem(self)
    class(materialItem), intent(inout) :: self

    ! Return static components to default
    self % name   = ''
    self % matIdx = 0
    self % T      = ZERO

    ! Deallocate allocatable components
    if(allocated(self % dens)) deallocate(self % dens)
    if(allocated(self % nuclides)) deallocate(self % nuclides)
    call self % extraInfo % kill()

  end subroutine kill_materialItem

  !!
  !! Prints the definition of material to the console
  !! Uses up to 60 columns
  !!
  !! Args:
  !!   None
  !! Errors:
  !!   None
  !!
  subroutine display_materialItem(self)
    class(materialItem), intent(in) :: self
    integer(shortInt)               :: i

    print '(A)', 'Material: '// trim(self % name) //' with index: ' // numToChar(self % matIdx)
    print '(A)', 'Temperature [K]: '//numToChar(self % T)
    print '(A)', 'Nuclide Composition:'
    print '(3A13, A20)', 'Atomic #', 'Mass #', 'Evaluation #', 'Density [1/barn/cm]'

    do i =1,size(self % nuclides)
      print '(3I13, ES20.10)', self % nuclides(i) % Z, self % nuclides(i) % A, self % nuclides(i) % T, &
                           self % dens(i)
    end do


  end subroutine display_materialItem


  !!
  !! Helper function to identify nuclide definition string
  !!
  !! Nuclide definition string has a following format:
  !!   ZZZAAA.TT
  !!
  !! ZZZ -> Up to 3 Digits   [0-9] that specify Atomic Number (minimum 1)
  !! AAA -> EXACTLY 3 Digits [0-9] that specify Mass Number
  !! TT  -> EXACTLY 2 Digits [0-9] that specify Evaluation Number
  !!
  !! Must also be left-adjusted and padded only with spaces.
  !!
  !! Args:
  !!   key [in] -> character string that may or may not contain Nuclide definition
  !!
  !! Result:
  !!   True if key matches the Nuclide Definition format. False otherwise
  !!
  function isNucDefinition(key) result(isIt)
    character(nameLen), intent(in) :: key
    logical(defBool)               :: isIt
    integer(shortInt)              :: dot, za, tt, L

    ! Save trim length of the string
    L = len_trim(key)

    ! Check that length is as expected
    if(L > 9 .or. L < 7) then ! Length cannot fit the format
      isIt = .false.
      return
    end if

    ! Find location of the dot
    dot = scan(key(1:L),'.')

    ! Verify that number of digits before and after dot
    za = verify(key(1:L),'0123456789')
    tt = verify(key(1:L),'0123456789', back = .true.)

    ! Verify that the location of the dot is consistent
    isIt = dot == za .and. dot == tt

  end function isNucDefinition

  !!
  !! Load information into nuclideInfo from string
  !!
  !! Takes ZZZAAA.TT string and converts it into Atomic, Mass and Evaluation number
  !!
  !! Args:
  !!   str [in] -> Input string in ZZZAAA.TT format
  !!
  !! Errors:
  !!   Returns fatal error if it fails to correctly convert string (e.g. string is not ZZZAAA.TT)
  !!
  subroutine init_nuclideInfo(self, str)
    class(nuclideInfo), intent(inout) :: self
    character(nameLen), intent(in)    :: str
    integer(shortInt)                 :: dot
    logical(defBool)                  :: flag
    character(100),parameter :: Here = 'init_nuclideInto (materialMenu_mod.f90)'

    if(.not.isNucDefinition(str)) then
      call fatalError(Here,'Input is not ZZZAAA.TT formated definition: '//trim(str))
    end if

    ! Find location of the dot
    dot = scan(str,'.')

    self % Z = charToInt(str(1:dot-4), error = flag)
    self % A = charToInt(str(dot-3:dot-1), error = flag )
    self % T = charToInt(str(dot+1:len_trim(str)), error = flag)

    if(flag) call fatalError(Here,'Failed to convert: '//trim(str)// ' to nuclide information')

  end subroutine init_nuclideInfo

  !!
  !! Convert nuclide information to the definition character
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Character in format ZZAAA.TT that describes nuclide definition
  !!
  !! Errors:
  !!   None
  !!
  elemental function toChar_nuclideInfo(self) result(str)
    class(nuclideInfo), intent(in) :: self
    character(nameLen)             :: str
    character(3)                   :: ZZ
    character(3)                   :: AAA
    character(2)                   :: TT

    write(ZZ, '(I3)')   self % Z
    write(AAA,'(I3.3)') self % A
    write(TT, '(I2.2)') self % T

    str = trim(adjustl(ZZ)) // AAA // "." // TT

  end function toChar_nuclideInfo

  !!
  !! Get pointer to a material definition under matIdx
  !!
  !! Args:
  !!   idx [in] -> Index of the material
  !!
  !! Result:
  !!   Pointer to a materialItem with the definition
  !!
  !! Errors:
  !!   FatalError if idx does not correspond to any defined material
  !!   FatalError if material definitions were not loaded
  !!
  function getMatPtr(idx) result(ptr)
    integer(shortInt), intent(in) :: idx
    type(materialItem), pointer   :: ptr
    character(100), parameter :: Here = 'getMatPtr (materialMenu_mod.f90)'

    ! Check if materialMenu is initialised
    if(.not.allocated(materialDefs)) then
      call fatalError(Here, "Material definitions were not loaded")
    end if

    ! Verify matIdx
    if( idx <= 0 .or. idx > nMat()) then
      call fatalError(Here,"matIdx: "//numToChar(idx)// &
                           " does not correspond to any defined material")
    end if

    ! Attach pointer
    ptr => materialDefs(idx)

  end function getMatPtr


  !!
  !! Return number of materials
  !!
  !! Args:
  !!   None
  !!
  !! Result:
  !!   Number of defined materials
  !!
  !! Errors:
  !!   Return 0 if materialMenu was not yet loaded
  !!
  function nMat() result(N)
    integer(shortInt) :: N

    if(allocated(materialDefs)) then
      N = size(materialDefs)
    else
      N = 0
    end if

  end function nMat

end module materialMenu_mod
