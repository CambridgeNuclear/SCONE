program testACE

  use numPrecision
  use aceCard_class, only : aceCard

  implicit none

  character(100)    :: filePath ='/home/mak60/myACE/acedata/92238JF311.ace'
  !character(100)    :: filePath ='/home/mak60/myACE/acedata/35081JEF311.ace'
  !character(100)    :: filePath ='/home/mak60/myACE/acedata/8016JEF311.ace'
  !character(100)    :: filePath ='/home/mak60/myACE/acedata/1001JEF311.ace'
  integer(shortInt) :: line = 1170286
  !integer(shortInt) :: line = 14296
  !integer(shortInt) :: line = 503141
  !integer(shortInt) :: line = 7681


  type(aceCard) :: ACE
  real(defReal), dimension(:), allocatable :: eGrid
  real(defReal),dimension(:),allocatable   :: xs
  integer(shortInt) :: i

  call ACE % readFromFile(filePath,line)

  call ACE % print()
  !do i=1,size(ACE % MTdata)
  !  call ACE % MTdata(i) % print()
  !end do

  !call ACE % ESZblock(eGrid,'energyGrid')
  !call ACE % ESZblock(xs,'totalXS')

  !do i=1,size(xs)
  !  print *, eGrid(i), xs(i)
  !end do


end program testACE
