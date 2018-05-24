program testACE

  use numPrecision
  use aceCard_class, only : aceCard

  implicit none

  character(100)    :: filePath ='/home/mak60/myACE/acedata/92238JF311.ace'
  integer(shortInt) :: line = 1170286
  type(aceCard) :: ACE
  real(defReal), dimension(:), allocatable :: eGrid
  real(defReal),dimension(:),allocatable   :: xs
  integer(shortInt) :: i

  call ACE % readFromFile(filePath,line)

  do i=1,size(ACE % MTdata)
    call ACE % MTdata(i) % print()
  end do

  !call ACE % ESZblock(eGrid,'energyGrid')
  !call ACE % ESZblock(xs,'totalXS')

  !do i=1,size(xs)
  !  print *, eGrid(i), xs(i)
  !end do


end program testACE
