program test

  use numPrecision
  use ByIsoNoMT_Data_class

  implicit none
  integer(kind=shortInt) :: i
  integer(kind=shortInt),dimension(2,3) :: A
  integer(kind=shortInt),dimension(:),pointer :: B
  type(ByIsoNoMT_Data)  :: CEdata
  character(len=20)      :: Input="./testInput"
  character(len=99)      :: format
  character(len=99),dimension(2) :: Ach


  call CEdata % readFrom(Input)
  call CEdata % print()


end program test
