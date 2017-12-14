program test

  use numPrecision
  use ByIsoNoMT_Data_class, only : byIsoNoMT_Data
  use aceNoMT_class

  implicit none
  integer(kind=shortInt) :: i
  integer(kind=shortInt),dimension(99),target :: A
  integer(kind=shortInt),dimension(:),pointer :: B
  INTEGER(SHORTiNT),DIMENSION(:), allocatable :: C
  type(ByIsoNoMT_Data)  :: CEdata
  character(len=pathLen)      :: matInput="./testInput"
  character(len=pathLen)      :: isoInput="/home/mak60/myACE/JEF311.aceXS"
  character(len=99)      :: format
  character(len=99),dimension(2) :: Ach
  character(len=pathLen)    :: acePath = "/home/mak60/myACE/acedata/92238JF311.ace"
  integer(shortInt)         :: firstLine = 1170286
  type(aceNoMT)             :: isotope


  C=[1,2,3,4,5,6,7,8,9,10]
  print *, C(3:5)

  call CEdata % readFrom(matInput,isoInput)
  call CEdata % print()

  call isotope % init(acePath,firstLine)

end program test
