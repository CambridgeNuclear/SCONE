program test

  use numPrecision
  use ByIsoNoMT_Data_class, only : byIsoNoMT_Data
  use aceNoMT_class

  use releseLawENDF_class
  use constantRelese_class, only : constantRelese
  use polynomialRelese_class, only : polynomialRelese

  implicit none
  integer(kind=shortInt) :: i
  integer(kind=shortInt),dimension(99),target :: A
  integer(kind=shortInt),dimension(:),pointer :: B
  INTEGER(SHORTiNT),DIMENSION(:), allocatable :: C
  type(ByIsoNoMT_Data)  :: CEdata
  character(len=pathLen)      :: matInput="./testInput"
  character(len=pathLen)      :: isoInput="C:\cygwin64\home\MikolajAdamKowalski\myACE\JEF311.aceXS"
  character(len=99)      :: format
  character(len=99),dimension(2) :: Ach
  character(len=pathLen)    :: acePath = "C:\cygwin64\home\MikolajAdamKowalski\myACE\acedata\92238JF311.ace"
  integer(shortInt)         :: firstLine = 1170286
  type(aceNoMT)             :: isotope
  real(defReal) :: kl

  class(releseLawENDF),pointer  :: relese

  real, pointer :: p1,p2,p3


  !C=[1,2,3,4,5,6,7,8,9,10]
  !print *, C(3:5)

  call CEdata % readFrom(matInput,isoInput)
  call CEdata % print()

  call isotope % init(acePath,firstLine)

  relese => constantRelese(3.9_defReal)
  print *, relese % releseAt(2.0_defReal)
  deallocate(relese)
  relese => polynomialRelese([0.0_8,0.0_8,0.0_8,1.0_8])

  print *, relese % releseAt(2.0_defReal),relese % releseAt(3.0_defReal),relese % releseAt(4.0_defReal)



end program test
