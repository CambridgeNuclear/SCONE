program test

  use numPrecision
  use genericProcedures
  use ByIsoNoMT_Data_class, only : byIsoNoMT_Data
  use aceNoMT_class

  use releaseLawENDF_class
  use constantRelease_class, only : constantRelease
  use polynomialRelease_class, only : polynomialRelease
  use tabularRelease_class, only : tabularRelease

  implicit none
  integer(kind=shortInt) :: i
  integer(longInt)       :: longI, longI2, rate
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
  real(defReal),dimension(:),allocatable :: R

  class(releaseLawENDF),pointer  :: release
  real(defReal), dimension(10) :: energy, second
  integer(shortInt), dimension(3) :: bounds, ENDF

  real, pointer :: p1,p2,p3


  !C=[1,2,3,4,5,6,7,8,9,10]
  !print *, C(3:5)

  call CEdata % readFrom(matInput,isoInput)
  call CEdata % print()

  call isotope % init(acePath,firstLine)

  R= [1.1_8, 2.3_8, 3.6_8, 9.6_8, 11.9_8]

  !print *, binaryFloorIdxC(R,1.1_8)

  C = [(2*i,i=1,10)]

 ! print *, linearCeilIdxO(C,17)

  !energy = [(real(i),i=1,20)]
  !second = [(i*5.0,i=1,20)]
  !bounds = [ 5, 10, 15, 20]
  !ENDF = [2,2,2,2]

  second = [1.0_8,1.5_8,2.0_8,1.5_8,1.7_8,2.5_8,2.3_8,2.0_8,2.1_8,1.5_8]
  energy = [(real(i),i=1,10)]
  bounds = [3,7,10]
  ENDF = [1,2,1]
  !energy(10)= energy(9)
  energy(4) = energy(3)
  energy(6) = energy(5)

  release => tabularRelease(energy,second,bounds,ENDF)


  do i=1,1000
    kl = (9.999-1.0)/1000.0 * i + 1.0
    print *, kl, release % releaseAt(kl)
  end do

  !print *, binarySearch(energy,1.0_8)

  !print *, linearFloorIdx(C,21)
!  call system_clock(count_rate=rate)
!  call system_clock(count=longI)
!  do i=1,1
!    kl= interpolate(0.00001_8, 1.0_8, 0.01_8, 10.0_8, 0.5_8)
!  end do
!  call system_clock(count=longI2)
!  print *, 'end', (longI2-longI)/real(rate)
end program test

