module ByIsoNoMT_Data_class

  use numPrecision
  use genericProcedures, only : fatalError, openToRead

  implicit none
  private

  type, public :: ByIsoNoMT_Data
    private
    ! Material Data
    character(len=matNameLen),dimension(:),allocatable :: matNames
    integer(kind=shortInt),dimension(:),allocatable    :: matNumIso
    integer(kind=shortInt),dimension(:,:),allocatable  :: matIsoIdx
    character(len=ZZidLen),dimension(:,:),allocatable  :: matIsoNames
    real(kind=defReal),dimension(:,:),allocatable      :: matIsoDens
    real(kind=defReal),dimension(:),allocatable        :: matTemp
  contains
    procedure :: readFrom
    procedure :: print

    procedure,private :: createMatArrays


  end type ByIsoNoMT_Data

contains

  subroutine readFrom(self, InputFile)
    class(ByIsoNoMT_Data), intent(inout) :: self
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

        read(unit = Input,    &
             fmt=*,           &
             iostat=readStat, &
             iomsg = readMsg) ZZid ,&
                              numDen

        if (readStat > 0) call fatalError(here, readMsg)
        if (numDen <= 0 ) call fatalError(here,'Numerical Density is -ve')

        self % matIsoNames(i,j) = ZZid
        self % matIsoDens(i,j)  = numDen

      end do
    end do

    close(Input)

  end subroutine readFrom

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

  end subroutine print
    
end module ByIsoNoMT_Data_class
