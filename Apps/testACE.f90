program testACE

  use numPrecision
  use aceCard_class, only : aceCard
  use correlatedLawENDF_inter, only : correlatedLawENDF
  use correlatedLawENDFfactory_func, only : new_correlatedLawENDF
  use RNG_class, only : RNG

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
  integer(shortInt) :: i, N, N_bins, idx, NS
  real(defReal)     :: E_in, E_out, E_max
  real(defReal)     :: mu
  class(correlatedLawENDF),allocatable :: correl
  type(RNG) :: random
  integer(defReal),dimension(:),allocatable :: bins

  call random % init(65858758587_8)

  call ACE % readFromFile(filePath,line)

  call ACE % print()
  stop
  allocate(correl, source = new_correlatedLawENDF(ACE,16))
  print *,'A=['

  E_in = 7.5_8
  mu = 0.0_8
  E_max = 0.120_8
  E_out = 0.02240507_8
  N = 100
  NS = 2000000
  N_bins = 32
  do i=1,N
    !E_out = E_max * i/N
    mu = -1.0 + 2.0 * i/N
    print *, correl % probabilityOf(mu,E_out,E_in), mu

  end do
  print *,'];'

  allocate(bins(N_bins))
  bins = 0
  i = 0
  !print *, 'K=['
  do while (i<NS)
    ! sample outgoing stat
    call correl % sample(mu,E_out,E_in,random)

    if( abs(E_out-0.02240507_8) < 0.001 ) then

      idx = ceiling((mu+ONE)/2 * N_bins)
      !idx = ceiling(E_out/E_max * N_bins)
        !print *, E_out, mu
        bins(idx) = bins(idx) + 1
        i = i+1

    end if
  end do
 ! print *, '];'


  print *, 'K=['
  do i=1,N_bins
    print *, bins(i), -ONE+TWO*i/N_bins
  end do
  print *, '];'


  !do i=1,size(ACE % MTdata)
  !  call ACE % MTdata(i) % print()
  !end do

  !call ACE % ESZblock(eGrid,'energyGrid')
  !call ACE % ESZblock(xs,'totalXS')

  !do i=1,size(xs)
  !  print *, eGrid(i), xs(i)
  !end do


end program testACE
