module materialDataNoMT_class

  use numPrecision
  use genericProcedures, only : fatalError
  use dictionary_class,  only : dictionary

  implicit none
  private

  !!
  !! Type to store all data relevant to material definition in neat little package.
  !! No data encapsulation (private attibute)! Other objects acces data directly.
  !! Should not be modified or used outside byNucNoMT
  !!
  type, public :: materialDataNoMT
    character(matNameLen)                       :: name     !! Material name
    real(defReal)                               :: temp     !! Material temperature [K]
    integer(shortInt)                           :: numNuc   !! Number of nuclides in material
    logical(defBool)                            :: isFissile =.false. !! True is contains fissile nuclides
    integer(shortInt),dimension(:),allocatable  :: nucIdx   !! Index Ids of Nuclides
    real(defReal),dimension(:),allocatable      :: nucDens  !! Nuclide densities 1/barn*cm
    character(ZZidLen),dimension(:),allocatable :: nucNames !! ZZids of nuclides
  contains
    procedure :: init
    procedure :: setNumberOfNuc
    procedure :: printMe
  end type materialDataNoMT

contains

  !!
  !! Initialise from dictionary. Does not read nuclide indexes (they are unavalible)
  !!
  subroutine init(self,dict,name)
    class(materialDataNoMT), intent(inout)      :: self
    class(dictionary), intent(in)               :: dict
    character(*),intent(in)                     :: name
    logical(defBool)                            :: isNotMaterial
    integer(shortInt)                           :: i
    real(defReal)                               :: temp
    character(nameLen),dimension(:),allocatable :: nucKeys
    character(100), parameter                   :: Here ='init (materialDataNoMT_class.f90)'

    ! Verify that provided dictionary contains material
    ! isNotMaterial = ('material' /= adjustl(dict % getChar('type') ))
    ! if (isNotMaterial) then
    !  call fatalError(Here,'Provided Dictionery is not a "material" it is: ' // dict % getChar('type') )
    !
    !end if

    ! Read required data
    self % name = name

    self % temp = dict % getReal('temp')

    ! Load nuclides Keys

    call dict % keys(nucKeys)

    call filterNotNuclides(nucKeys)

    call self % setNumberOfNuc(size(nucKeys))
    self % nucNames = nucKeys

    ! Load nuclide densities
    do i=1,size(nucKeys)
      temp = dict % getReal(nucKeys(i))
      if (temp < 0) call fatalError(Here,'Density of nuclide: '// nucKeys(i) //' is negative')
      self % nucDens(i) = temp

    end do


    contains

      !!
      !! Remove all keywords that are not a valid ZZid from keys
      !!
      subroutine filterNotNuclides(keys)
        character(*),dimension(:),allocatable, intent(inout) :: keys
        character(len(keys)), dimension(:), allocatable      :: tempKeys
        logical(defBool),dimension(:),allocatable            :: mask
        integer(shortInt) :: i

        ! Create mask array
        mask = isZZid(keys)

        tempKeys = pack(keys,mask)

        ! Switch allocation
        deallocate(keys)
        call move_alloc(tempKeys,keys)

      end subroutine filterNotNuclides

      !!
      !! Check is a character is a valid ZZId
      !!
      elemental function isZZid(key) result(isValidZZid)
        character(*), intent(in)   :: key
        character(:),allocatable   :: keyTemp
        integer(shortInt)          :: L, dot_idx, idx
        logical(defBool)           :: isValidZZid

        ! Load into temporary variable with no leading or trailing blanks
        keyTemp = trim(adjustl(key))

        ! Check that length of character is approperiate
        L = len_trim(keyTemp)
        isValidZZid = (L == 8) .or. (L == 9)
        if (.not.isValidZZid) return

        ! Check that position of a dot is approperiate
        dot_idx = index(keyTemp,'.')
        isValidZZid = (dot_idx == 5) .or. (dot_idx == 6)
        if (.not.isValidZZid) return

        ! Check that there are no charcters left of the dot
        idx = verify(keyTemp(1:dot_idx-1),'0123456789')
        isValidZZid = (idx == 0)
        if (.not.isValidZZid) return

        ! Check that dot is followed by 2 numbers
        idx = verify(keyTemp(dot_idx+1:dot_idx+2),'0123456789')
        isValidZZid = (idx == 0)
        if (.not.isValidZZid) return

        ! Check that charcter does not end with number
        idx = verify(keyTemp(L:L),'0123456789')
        isValidZZid = (idx /= 0 )

      end function isZZid

  end subroutine init


  !!
  !! Allocate space for the nuclide data for material
  !! Deallocates any privous data and sets numNuc in the objest
  !!
  subroutine setNumberOfNuc(self,numNuc)
    class(materialDataNoMT), intent(inout) :: self
    integer(shortInt),intent(in)           :: numNuc

    if(allocated( self % nucIdx  )) deallocate (self % nucIdx  )
    if(allocated( self % nucDens )) deallocate (self % nucDens )
    if(allocated( self % nucNames)) deallocate (self % nucNames)

    self % numNuc = numNuc

    allocate(self % nucIdx  (numNuc))
    allocate(self % nucDens (numNuc))
    allocate(self % nucNames(numNuc))

    ! Set default values
    self % nucIdx   = 0
    self % nucDens  = 0.0
    self % nucNames = ''

  end subroutine setNumberOfNuc

  subroutine printMe(self)
    class(materialDataNoMT), intent(in) :: self
    character(100)                      :: format1, format2, line


    format1 = '(100G20.10)'
    format2 = '(A100)'
    line = repeat('<>',100)


    print format2, line
    print format1, self % name, "Temp: ", self % temp, "Contains ", self % numNuc, " nuclides"
    print format1, "isFissile is:", self % isFissile
    print format1, "Nuclide Indexes   :", self % nucIdx
    print format1, "Nuclide Names     :", self % nucNames
    print format1, "Nuclide Densities :", self % nucDens


    print format2, line

  end subroutine printMe

end module materialDataNoMT_class
