module ByIsoNoMT_Data_class

  use numPrecision
  use genericProcedures, only : fatalError, openToRead, removeDuplicates, linSearch

  implicit none
  private

  type, public :: byIsoNoMT_Data
    private
    ! Material Data
    character(len=matNameLen),dimension(:),allocatable :: matNames
    integer(kind=shortInt),dimension(:),allocatable    :: matNumIso
    integer(kind=shortInt),dimension(:,:),allocatable  :: matIsoIdx
    character(len=ZZidLen),dimension(:,:),allocatable  :: matIsoNames
    real(kind=defReal),dimension(:,:),allocatable      :: matIsoDens
    real(kind=defReal),dimension(:),allocatable        :: matTemp
    ! Isotope Data
    character(len=zzIdLen),dimension(:),allocatable    :: isoNames
  contains
    procedure :: readFrom
    procedure :: print

    procedure,private :: createMatArrays
    procedure,private :: readMaterials
    procedure,private :: createIsotopeList
    procedure,private :: assignIsoIndices
    procedure,private :: readIsotopes

  end type byIsoNoMT_Data

contains

  subroutine readFrom(self,matInput,isotopeLib)
  !! Reads materials and isotopes from provided material and ACE library input files
    class(byIsoNoMT_Data), intent(inout) :: self
    character(len=*), intent(in)         ::  matInput
    character(len=*), intent(in)         ::  isotopeLib

    call self % readMaterials(matInput)
    call self % createIsotopeList()
    call self % assignIsoIndices()
    call self % readIsotopes(isotopeLib)

  end subroutine readFrom


  subroutine readMaterials(self, inputFile)
    class(byIsoNoMT_Data), intent(inout) :: self
    character(len=*), intent(in)         :: inputFile

    integer(kind=shortInt),parameter     :: input=66
    character(len=99),parameter          :: here='readFrom in ByIsoNoMT_Data_class.f03'
    character(len=99)                    :: readMsg

    character(len=matNameLen)            :: matName
    character(len=ZZidLen)               :: ZZid
    integer(kind=shortInt)               :: numMat, numIso, maxIso
    real(kind=defReal)                   :: temp, numDen
    integer(kind=shortInt)               :: i, j, readStat

    call openToRead(Input,InputFile)

    call readMaxNumIso(Input,maxIso)

    ! Read Number of Materials in input file
    read(unit = Input, fmt=* , iostat = readStat, iomsg = readMsg) numMat
    if (readStat > 0) call fatalError(Here, readMsg)

    call self % createMatArrays(numMat,maxIso)

    ! Read Materials
    do i=1,numMat

      read(unit = Input,    &
           fmt=*,           &
           iostat=readStat, &
           iomsg = readMsg) matName, &
                            temp,    &
                            numIso

      if (readStat > 0) call fatalError(here, readMsg)
      if (numIso <= 0) call fatalError(here,'Number of Material Isotopes is 0 or -ve')

      self % matNames(i)  = matName
      self % matTemp(i)   = temp
      self % matNumIso(i) = numIso

      do j=1,numIso

        read(unit = Input,      &
             fmt=*,             &
             iostat = readStat, &
             iomsg = readMsg)   ZZid ,&
                                numDen

        if (readStat > 0) call fatalError(here, readMsg)
        if (numDen <= 0 ) call fatalError(here,'Numerical Density is -ve')

        self % matIsoNames(i,j) = ZZid
        self % matIsoDens(i,j)  = numDen
      end do
    end do

    ! Give error if there are more enteries in the input
    read (unit = Input,    &
          fmt =*,          &
          iostat = readStat) ZZid, numDen
    if (readStat /= endOfFile) call fatalError(here,'End of file was not reached, check the ' // &
                                                    'number of isotopes in the last material')

    close(Input)

  end subroutine readMaterials

  subroutine createIsotopeList(self)
    class(byIsoNoMT_Data),intent(inout)              :: self
    integer(shortInt)                           :: maxIsoNames
    character(len=zzIdLen),dimension(:),allocatable  :: withRepetition
    integer(kind=shortInt)                           :: i,j

    maxIsoNames=sum(self % matNumIso)

    ! Crate array of isotope names with repetitions
    allocate(withRepetition(maxIsoNames))
    j=1
    do i = 1,size(self % matNames)
      ! Load material names for material i
      withRepetition(j:j+self % matNumIso(i)) = self % matIsoNames(i,1:self % matNumIso(i))
      j = j + self % matNumIso(i)
    end do

    self % isoNames = removeDuplicates(withRepetition)
  end subroutine

  subroutine assignIsoIndices(self)
    class(byIsoNoMT_Data),intent(inout)       :: self
    integer(kind=shortInt)                    :: i,j

    do i = 1,size(self % matNames)
      do j = 1,self % matNumIso(i)
        self % matIsoIdx(i,j) = linSearch(self % isoNames, self % matIsoNames(i,j))
        if (self % matIsoIdx(i,j) == -1 ) then
          call fatalError('assignIsoIndices (byIsoNoMT_class.f90)', &
                          'Isotope ' // self % matIsoNames(i,j) //' was not found')
        end if
      end do
    end do
  end subroutine

  subroutine readIsotopes(self,libraryPath)
    class(byIsoNoMT_Data), intent(inout)             :: self
    character(len=*),intent(in)                      :: libraryPath
    integer(kind=shortInt),parameter                 :: library=78
    character(len=99)                                :: readMsg
    character(len=zzIdLen),dimension(:),allocatable  :: zzIDs
    integer(kind=shortInt),dimension(:),allocatable  :: startLine
    character(len=pathLen),dimension(:),allocatable  :: isoPath
    integer(kind=shortInt)                           :: i, j, readStat
    integer(kind=shortInt)                           :: libLen

    call openToRead(library,libraryPath)

    ! Find length of isotope library
    libLen=0
    do while (readStat /= -1)
      read(unit = library, fmt=*,iostat=readStat,iomsg = readMsg)
      libLen = libLen + 1
    end do
    rewind(library)

    ! Allocate and read library
    allocate(zzIDs(libLen))
    allocate(startLine(libLen))
    allocate(isoPath(libLen))

    do i=1,libLen
      read(library,"(A10 I12 A100)" ) zzIds(i), startLine(i), isoPath(i)
    end do

    ! Read Isotope Data
    do i=1,size(self % isoNames)
      j = linSearch(zzIds,self % isoNames(i))
      if (j == -1) then
        call fatalError('readIsotopes (byIsoNoMT_Data_class.f90)', &
                        'Isotope ' // self % isoNames(i) //' was not found')
      end if
      print *, zzIds(j), startLine(j), isoPath(j) ! Later will call isotopeACE to initialise
    end do


    close(library)

  end subroutine readIsotopes

  subroutine createMatArrays(self,numMat,maxIso)
    class(ByIsoNoMT_Data), intent(inout)  :: self
    integer(kind=shortInt),intent(in)     :: numMat, maxIso

    allocate(self % matNames(numMat))
    allocate(self % matNumIso(numMat))
    allocate(self % matIsoIdx(numMat,maxIso))
    allocate(self % matIsoNames(numMat,maxIso))
    allocate(self % matIsoDens(numMat,maxIso))
    allocate(self % matTemp(numMat))

  end subroutine createMatArrays

  subroutine readMaxNumIso(Input,maxIso)
    !class(ByIsoNoMT_Data), intent(inout)  :: self
    integer(kind=shortInt),intent(in)     :: Input
    integer(kind=shortInt),intent(out)    :: maxIso
    character(len=99),parameter           :: Here='readMaxNumIso in ByIsoNoMT_Data_class.f03'
    character(len=99)                     :: readMsg
    integer(kind=shortInt)                :: numMat, numIso, readStat, i, j
    character(len=3)                      :: dummyChar
    real(kind=defReal)                    :: dummyReal

    rewind(Input)

    read(unit = Input, fmt=* , iostat = readStat, iomsg = readMsg) numMat
    if (readStat > 0) call fatalError(Here, readMsg)

    maxIso = -1
    do i=1,numMat
      read(unit = Input, fmt=*, iostat=readStat, iomsg = readMsg) dummyChar, &
                                                                  dummyReal,    &
                                                                  numIso
      if (readStat > 0) call fatalError(Here, readMsg)
      if (numIso <= 0) call fatalError(Here,'Number of Material Isotopes is 0 or -ve')
      maxIso=max(maxIso,numIso)
      ! Skip lines to get to next material header
      do j=1,numIso
        read(Input,*) dummyChar, dummyReal
      end do
    end do

    rewind(Input)

  end subroutine readMaxNumIso

  subroutine print(self)
    class(byIsoNoMT_Data), intent(in) :: self
    character(len=99)                 :: format
    integer(kind=shortInt)            :: i

    print '(a)', 'Material Names:'
    print '(a)', self % matNames

    print '(a)', 'Isotope Numbers:'
    print '(i5)', self % matNumIso

    print '(a)', 'Materials Temperatures:'
    print '(f10.3)', self % matTemp

    print '(a)', 'Isotope Names:'
    do i=1,size(self % matIsoNames,1)
      print '(9999a)', self % matIsoNames(i,:)
    end do

    print '(a)', 'Isotope Indexes:'
    do i=1,size(self % matIsoNames,1)
      print '(9999I5)', self % matIsoIdx(i,:)
    end do

    print '(a)', 'Isotope Densities:'
    do i=1,size(self % matIsoNames,1)
      print '(9999es15.5)', self % matIsoDens(i,:)
    end do

    print '(a)', 'Isotope Names Array'
    print '(a)', self % isoNames

  end subroutine print
    
end module ByIsoNoMT_Data_class
