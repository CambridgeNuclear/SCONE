module byNucMT_Data_class

  use numPrecision
  use genericProcedures,      only : fatalError, openToRead, removeDuplicates, linFind, &
                                     findDuplicates, arrayConcat, targetNotFound
  use dictionary_class,       only : dictionary

  use aceMT_class,            only : aceMT
  use materialDataNoMT_class, only : materialDataNoMT


  implicit none
  private

  !!
  !! "Fat" data block for XS data. All data of large size is collected here.
  !! In parallel calculations it is supposed to be shared between diffrent threads
  !! so it needs to be strictly read only.
  !!
  type, public :: byNucMT_Data
    !private
    ! Material Data
    type(materialDataNoMT),dimension(:),pointer :: matData     => null()

    ! Isotope Data
    character(nameLen),dimension(:),pointer    :: nucNames    => null()
    type(aceMT),dimension(:), pointer          :: nucXsData   => null()

  contains
    procedure :: init
    procedure :: printMe

    procedure,private :: readMaterials
    procedure,private :: createNuclideList
    procedure,private :: assignNucIndices
    procedure,private :: readNuclides
    procedure,private :: setFissileMaterials

  end type byNucMT_Data

contains


  !!
  !! Subroutine to build byNucMT_Data from dictionary
  !!
  subroutine init(self,matDict)
    class(byNucMT_Data), intent(inout)          :: self
    class(dictionary), intent(in)               :: matDict
    character(pathLen)                          :: nuclideLib

    call self % readMaterials(matDict)

    call self % createNuclideList()

    call self % assignNucIndices()

    ! Read
    call matDict % get(nuclideLib, 'aceLibrary')
    call self % readNuclides(nuclideLib)

    call self % setFissileMaterials

  end subroutine init

  !!
  !! Reads materials data from dictionary.
  !! All data except nuclide indexes is assigned.
  !!
  subroutine readMaterials(self,matDict)
    class(byNucMT_Data), intent(inout)        :: self
    class(dictionary), intent(in)               :: matDict
    type(dictionary)                            :: locDict
    character(nameLen)                          :: matName
    integer(shortInt)                           :: nMat, i
    character(nameLen),dimension(:),allocatable :: matNames

    ! Load material names
    call matDict % keysDict(matNames)

    ! Find number of materials
    nMat = size(matNames )

    ! Allocate material data
    allocate(self % matData(nMat))

    ! Load individual material data
    do i=1,nMat
      matName = matNames(i)
      call matDict % get(locDict, matName)
      call self % matData(i) % init(locDict, matName)

    end do

  end subroutine readMaterials


!  !!
!  !! Check if the dictionary has type 'material'
!  !!
!  function isItMaterial(dict) result(isIt)
!    class(dictionary), intent(in)   :: dict
!    logical(defBool)                :: isIt
!
!    isIt = ('material' == adjustl(dict % getChar('type') ))
!
!  end function isItMaterial

  !!
  !! Remove repetitions from material definitions and create a list of all used nuclide data
  !! cards without any repetitions
  !!
  subroutine createNuclideList(self)
    class(byNucMT_Data),intent(inout)           :: self
    integer(shortInt)                             :: maxNucNames
    character(nameLen),dimension(:),allocatable   :: withRepetition
    character(nameLen),dimension(:),allocatable   :: noRepetition
    integer(shortInt)                             :: i,j

    maxNucNames=sum(self % matData(:) % numNuc)

    ! Crate array of isotope names with repetitions
    allocate(withRepetition(maxNucNames))
    j=1
    do i = 1,size(self % matData)
      ! Load material names for material i
      withRepetition(j:j+self % matData(i) % numNuc) = self % matData(i) % nucNames
      j = j + self % matData(i) % numNuc
    end do

    noRepetition = removeDuplicates(withRepetition)
    allocate(self % nucNames( size(noRepetition) ))
    self % nucNames = noRepetition

  end subroutine createNuclideList


  !!
  !! Assign nuclide indexes to material nuclides on the unified nuclide list.
  !!
  subroutine assignNucIndices(self)
    class(byNucMT_Data),intent(inout)  :: self
    integer(shortInt)                    :: i,j

    do i = 1,size(self % matData)
      do j = 1,self % matData(i) % numNuc
        self % matData(i) % nucIdx(j) = linFind(self % nucNames, self % matData(i) % nucNames(j))
        if (self % matData(i) % nucIdx(j) == -1 ) then
          call fatalError('assignIsoIndices (byNucMT_class.f90)', &
                          'Isotope ' // self % matData(i) % nucNames(j) //' was not found')
        end if
      end do
    end do
  end subroutine assignNucIndices

  !!
  !! Read library of ACE nuclide cards and read nuclide data
  !!
  subroutine readNuclides(self,libraryPath)
    class(byNucMT_Data), intent(inout)         :: self
    character(*),intent(in)                      :: libraryPath
    integer(shortInt),parameter                  :: library=78
    character(99)                                :: readMsg
    character(nameLen),dimension(:),allocatable  :: zzIDs
    integer(shortInt),dimension(:),allocatable   :: startLine
    character(pathLen),dimension(:),allocatable  :: nucPath
    integer(shortInt)                            :: i, j, readStat
    integer(shortInt)                            :: libLen
    character(nameLen)                           :: currentZZId, &
                                                    pastZZId
    call openToRead(library,libraryPath)

    ! Find length of isotope library
    ! There is a discrepancy in behaviour on Linux and Windows-Cygwin system.
    ! On Windows after readeing last line one addoitional read is made in which last line is read
    ! again but with iostat=-1. On lunux iostat = -1 is returned when READING LAST LINE THE FIRST TIME.
    ! Thus if code that runs porperly on linux is reused on windows last line is read twice!. Explicit
    ! check whether that happens is required to obtain correct behaviour on both systems.
    libLen=0
    do
      read(unit = library, fmt=*,iostat=readStat,iomsg = readMsg) currentZZId
      if(readStat == -1 .and. trim(currentZZId)==trim(pastZZId)) exit ! See character equality to indicate double read of last line
      libLen = libLen + 1
      pastZZId = currentZZId
    end do
    rewind(library)

    ! Allocate and read library
    allocate(zzIDs(libLen))
    allocate(startLine(libLen))
    allocate(nucPath(libLen))

    do i=1,libLen
      read(library,"(A10, I12, A100)" ) zzIds(i), startLine(i), nucPath(i)
    end do

    ! Check library for repeted zzIDs
    if (size(zzIds) /= size(removeDuplicates(zzIds))) then
      call fatalError('readIsotopes (byNucMT_Data_class.f90)', &
                      'Duplicate zzIds found in ACE library file: ' //&
                      arrayConcat(findDuplicates(zzIds)))
    end if

    ! Allocate Memory for isotopic data
    allocate(self % nucXsData(size(self % nucNames)))

    ! Read Isotope Data
    do i=1,size(self % nucNames)
      ! **** This Search need to be modernised ****!
      j = linFind(zzIds,self % nucNames(i))
      if (j == targetNotFound) then
        call fatalError('readIsotopes (byNucMT_Data_class.f90)', &
                        'Isotope ' // self % nucNames(i) //' was not found')
      end if

      print *, "Reading : ", zzIds(j), startLine(j), nucPath(j)
      call self % nucXsData(i) % init(nucPath(j),startLine(j))

    end do
    print *, size(self % nucXSData)


    close(library)

  end subroutine readNuclides


!  !!
!  !! Read input file and indentify maximum number of nuclides in the problem
!  !!
!  subroutine readMaxNumNuc(Input,maxNuc)
!    integer(shortInt),intent(in)      :: Input
!    integer(shortInt),intent(out)     :: maxNuc
!    character(99),parameter           :: Here='readMaxNumNuc (byNucMT_Data_class.f90)'
!    character(99)                     :: readMsg
!    integer(shortInt)                 :: numMat, numNuc, readStat, i, j
!    character(3)                      :: dummyChar
!    real(defReal)                     :: dummyReal
!
!    rewind(Input)
!
!    read(unit = Input, fmt=* , iostat = readStat, iomsg = readMsg) numMat
!    if (readStat > 0) call fatalError(Here, readMsg)
!
!    maxNuc = -1
!    do i=1,numMat
!      read(unit = Input, fmt=*, iostat=readStat, iomsg = readMsg) dummyChar, &
!                                                                  dummyReal,    &
!                                                                  numNuc
!      if (readStat > 0) call fatalError(Here, readMsg)
!      if (numNuc <= 0)  call fatalError(Here,'Number of Material Nuclides is 0 or -ve')
!      maxNuc=max(maxNuc,numNuc)
!
!      ! Skip lines to get to next material header
!      do j=1,numNuc
!        read(Input,*) dummyChar, dummyReal
!      end do
!
!    end do
!
!    rewind(Input)
!
!  end subroutine readMaxNumNuc

  !!
  !! Set isFissile flag in materials containing fissile nuclides to .true.
  !!
  subroutine setFissileMaterials(self)
    class(byNucMT_Data), intent(inout) :: self
    integer(shortInt)                    :: i,j,nucIdx

    do i=1,size(self % matData)
      self % matData(i) % isFissile =.false.
      do j=1,size(self % matData(i) % nucIdx)
        nucIdx = self % matData(i) % nucIdx(j)
        self % matData(i) % isFissile = self % matData(i) % isFissile .or. self % nucXSData(nucIdx) % isFissile

      end do
    end do

  end subroutine setFissileMaterials



  !!
  !! Print data about all materials to the console
  !!
  subroutine printMe(self)
    class(byNucMT_Data), intent(in) :: self
    integer(shortInt)                 :: i

    do i =1,size(self % matData)
      call self % matData(i) % printMe()
    end do

  end subroutine printMe
    
end module byNucMT_Data_class
