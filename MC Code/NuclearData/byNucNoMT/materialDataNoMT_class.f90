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
    integer(shortInt),dimension(:),allocatable  :: nucIdx   !! Index Ids of Nuclides
    real(defReal),dimension(:),allocatable      :: nucDens  !! Nuclide densities 1/barn*cm
    character(ZZidLen),dimension(:),allocatable :: nucNames !! ZZids of nuclides
  contains
    procedure :: init
    procedure :: setNumberOfNuc
    procedure :: print
  end type materialDataNoMT

contains

  !!
  !! Initialise from dictionary. Does not read nuclide indexes (they are unavalible)
  !!
  subroutine init(self,dict,name)
    class(materialDataNoMT), intent(inout) :: self
    class(dictionary), intent(in)          :: dict
    character(*),intent(in)                :: name
    logical(defBool)                       :: isNotMaterial
    character(100), parameter              :: Here ='init (materialDataNoMT_class.f90)'

    ! Verify that provided dictionary contains material
    isNotMaterial = ('material' == adjustl(dict % getChar('type') ))
    if (isNotMaterial) then
      call fatalError(Here,'Provided Dictionery is not a "material" it is:' // dict % getChar('type') )

    end if

    ! Read required data
    self % name = name
   ! self % temp = dict % getReal(

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

  subroutine print(self)
    class(materialDataNoMT), intent(in) :: self
    character(100)                      :: format1, format2, line


    format1 = '(100G20.10)'
    format2 = '(A100)'
    line = repeat('<>',100)


    print format2, line
    print format1, self % name, "Temp: ", self % temp, "Contains ", self % numNuc, " nuclides"

    print format1, "Nuclide Indexes   :", self % nucIdx
    print format1, "Nuclide Names     :", self % nucNames
    print format1, "Nuclide Densities :", self % nucDens


    print format2, line

  end subroutine print

end module materialDataNoMT_class
