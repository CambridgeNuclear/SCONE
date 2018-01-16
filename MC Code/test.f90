program test

  use numPrecision
  use RNG_class
  use genericProcedures
  use ByIsoNoMT_Data_class, only : byIsoNoMT_Data
  use aceNoMT_class

  use endfTable_class, only : endfTable

  use releaseLawENDF_class
  use constantRelease_class, only : constantRelease
  use polynomialRelease_class, only : polynomialRelease
  use tabularRelease_class, only : tabularRelease

  use miuEndfPdf_class,   only : miuEndfPdf, miuEndfPdf_ptr
  use isotropicMiu_class, only : isotropicMiu
  use equiBin32Miu_class, only : equiBin32Miu
  use tabularMiu_class,   only : tabularMiu
  use tabularPdf_class, only : tabularPdf


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
  real(defReal) :: kl, eps
  real(defReal),dimension(:),allocatable :: R, x, pdf
  type(RNG) :: random


  class(releaseLawENDF),pointer  :: release
  class(endfTable),pointer       :: eTable

  real(defReal), dimension(20) :: energy, second
  integer(shortInt), dimension(4) :: bounds, ENDF

  real, pointer :: p1,p2,p3

  type(miuEndfPdf_ptr) :: myPtr, myPtr2
  type(tabularPdf) :: table

  !C=[1,2,3,4,5,6,7,8,9,10]
  !B => C(1:8)
  !print *, B(3:5)

  call CEdata % readFrom(matInput,isoInput)
  call CEdata % print()

  call isotope % init(acePath,firstLine)

  !R= [1.1_8, 2.3_8, 3.6_8, 9.6_8, 11.9_8]

  !print *, binaryFloorIdxC(R,1.1_8)

  !C = [(2*i,i=1,10)]

  R = [(-1.0+2.0/32*i,i=0,32)]

  R(29)=0.8

  myPtr = equiBin32Miu(R)
  myPtr2 = myPtr

  energy = [(real(i),i=1,20)]
  second = [(i*5.0,i=1,20)]
  bounds = [ 5, 10, 15, 20]
  ENDF = [5,1,3,4]
   !bounds = [20]
   !ENDF = [1]


!  second = [1.0_8,1.5_8,2.0_8,1.5_8,1.7_8,2.5_8,2.3_8,2.0_8,2.1_8,1.5_8]
!  energy = [(real(i),i=1,10)]
!  bounds = [3,7,10]
!  ENDF = [1,2,1]
!  !energy(10)= energy(9)
!  energy(4) = energy(3)
!  energy(6) = energy(5)
!

!  release => tabularRelease(energy,second,bounds,ENDF)!
!
!  eTable  => endfTable(energy,second,bounds,ENDF)
!
!  do i=1,100
!    kl = (19.99-1.0)/100.0 * i + 1.0
!    print *, kl, release % releaseAt(kl), eTable % at(kl), release % releaseAt(kl)-eTable % at(kl)
!  end do


  !x = [(-1.0+2.0/10*i,i=0,10)]
  !pdf = abs(x)
!
  x = [ -1.0_8, 0.0_8, 1.0_8]
!  pdf = [ 1.0_8/3 ,2.0_8/3, 76876.0_8]

  !call table % init(x,pdf,0)

 ! myPtr = tabularMiu(x,pdf,0)


 ! do i=0,4000
    !kl = 2.0/1000 * i - 1.0
     !eps = random % get()
     !print *, myPtr % sample(random)
    !print *, kl, myPtr % probabilityOf(kl)
  !end do



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

