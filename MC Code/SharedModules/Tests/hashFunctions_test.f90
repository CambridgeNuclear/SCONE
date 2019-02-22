module hashFunctions_test

  use numPrecision
  use hashFunctions_func, only : FNV_1, knuthHash
  use pfUnit_mod

  implicit none

contains

  !!
  !! Test computation of FNV_1 hashes on short and long integer
  !!
@Test
  subroutine testFNV1()
    integer(shortInt) :: hashVal1, hashResult1

    ! Test for short int
    call FNV_1('All hail gfortran!', hashVal1)
    hashResult1 = transfer(z'03deced0', shortInt)
    @assertEqual(hashResult1, hashVal1)

    call FNV_1('MY HASH IS BROKEN', hashVal1)
    hashResult1 = transfer(z'867ff812', shortInt)
    @assertEqual(hashResult1, hashVal1)

  end subroutine testFNV1

  !!
  !! Test knuth hash
  !!
@Test
  subroutine testKnuthHash()
    integer(shortInt), dimension(100) :: keys
    integer(shortInt), dimension(100) :: hashes
    integer(shortInt)                 :: i, m

    ! Generate some random keys of shortInts
    keys(1) = 765875
    do i=2,100
      keys(i) = keys(i-1) * 22695477
    end do

    ! Case 1 range 2**m m = 4
    m = 4
    hashes = knuthHash(keys, m)
    @assertLessThanOrEqual(0, hashes)
    @assertGreaterThan(2**m, hashes)

    ! Case 2 range m = 30
    m = 30
    hashes = knuthHash(keys, m)
    @assertLessThanOrEqual(0, hashes)
    @assertGreaterThan(2**m, hashes)

    ! Case 3 range m = 37 (will be clipped to m = 31)
    m = 37
    hashes = knuthHash(keys, m)
    @assertLessThanOrEqual(0, hashes)
    @assertGreaterThan(huge(m), hashes)

    ! Case 4 -ve m (will be clipped to m=1
    m = -71
    hashes = knuthHash(keys, m)
    @assertLessThanOrEqual(0, hashes)
    @assertGreaterThan(2, hashes)

  end subroutine testKnuthHash

    
end module hashFunctions_test
