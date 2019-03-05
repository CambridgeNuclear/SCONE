program produceRNG
  use numPrecision
  use RNG_class, only : RNG

  implicit none
  type(RNG)         :: rand1
  integer(shortInt) :: i


  call rand1 % init(19_longInt)
  open(99, file='randSeq', form='unformatted', access='direct', recl=8, status='replace')

  do i=1,1

    write(99, rec=i) rand1 % getInt()

  end do

  close(99)
end program produceRNG
