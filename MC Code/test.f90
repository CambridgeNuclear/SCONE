program test

  use numPrecision
  use ByIsoNoMT_Data_class, only : byIsoNoMT_Data

  implicit none
  integer(kind=shortInt) :: i
  integer(kind=shortInt),dimension(2,3) :: A
  integer(kind=shortInt),dimension(:),pointer :: B
  type(ByIsoNoMT_Data)  :: CEdata
  character(len=pathLen)      :: matInput="./testInput"
  character(len=pathLen)      :: isoInput="/home/mak60/myACE/JEF311.aceXS"
  character(len=99)      :: format
  character(len=99),dimension(2) :: Ach


  call CEdata % readFrom(matInput,isoInput)
  call CEdata % print()


end program test
