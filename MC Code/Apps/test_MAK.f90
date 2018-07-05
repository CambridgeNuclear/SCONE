program test

  use numPrecision
  use RNG_class
  use genericProcedures
  use hashFunctions_func, only : unSig64Int, FNV_1


  implicit none

  type :: myType
    !real(defReal),allocatable,dimension(:) :: A
    !integer(4) :: B
    !integer(4) :: D
    !integer(8) :: C
    real(defReal),dimension(:),allocatable :: A
    !real(defReal),dimension(9) :: A
  end type

  integer(shortInt)  :: i,j,k
  integer(shortInt)  :: N_uni, N_cell, N_lat
  character(4)      :: Ent
  character(nameLen) :: C1,C2,C3
  integer(16)      :: longlongInt
  integer(longInt) :: hashedKey
  type(myType) :: A

  !allocate(A % A(100))

  print *,storage_size(A)
  stop

  N_uni  = 1000
  N_cell = 1000
  N_lat  = 100000

 ! print '(z20)', FNV_1('123')

  Ent ='aaav'
  ! Loop across universes
  do i=1,N_uni

    ! Loop across cells
    do j=1,N_cell
      !Loop across ijkIdx
      do k=1,N_lat

        ! Calculate and print hash
        !write(ent,'(I4.4, I10.10, I10.10)') i!, j, k
        !print '(A30, z20)',ent, FNV_1(Ent)
        hashedKey = FNV_1(Ent)
        !print *, hashedKey
      end do
    end do
  end do



end program test


