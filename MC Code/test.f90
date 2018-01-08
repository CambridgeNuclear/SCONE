program test

  use numPrecision
  use ByIsoNoMT_Data_class, only : byIsoNoMT_Data
  use aceNoMT_class

  use releaseLawENDF_class
  use constantRelease_class, only : constantRelease
  use polynomialRelease_class, only : polynomialRelease

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
  real(defReal) :: kl

  class(releaseLawENDF),pointer  :: release

  real, pointer :: p1,p2,p3


  !C=[1,2,3,4,5,6,7,8,9,10]
  !print *, C(3:5)

  call CEdata % readFrom(matInput,isoInput)
  call CEdata % print()

  call isotope % init(acePath,firstLine)

  release => constantRelease(3.9_defReal)
  print *, release % releaseAt(2.0_defReal)
  deallocate(release)
  release => polynomialRelease([-4.0_8,0.0_8,0.0_8,1.0_8])

  print *, release % releaseAt(2.0_defReal),release % releaseAt(3.0_defReal),release % releaseAt(4.0_defReal)



end program test
