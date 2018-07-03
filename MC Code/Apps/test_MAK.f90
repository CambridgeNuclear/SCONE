program test

  use numPrecision
  use RNG_class
  use genericProcedures
  use hashFunctions_func, only : unSig64Int, FNV_1

  implicit none
  integer(shortInt)  :: i,j,k
  integer(shortInt)  :: N_uni, N_cell, N_lat
  character(12)      :: Ent
  integer(16)      :: longlongInt
  integer(longInt) :: hashedKey

  N_uni  = 10
  N_cell = 1000
  N_lat  = 5000

  ! Loop across universes
  do i=1,N_uni

    ! Loop across cells
    do j=1,N_cell
      !Loop across ijkIdx
      do k=1,N_lat

        ! Calculate and print hash
        Ent = transfer([i,j,k],Ent)
        hashedKey = FNV_1(Ent)
       ! print '(z35, z20)', transfer(Ent,longlongInt), FNV_1(Ent)

      end do
    end do
  end do



end program test


