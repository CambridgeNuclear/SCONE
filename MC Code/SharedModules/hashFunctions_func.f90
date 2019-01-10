module hashFunctions_func

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private


  !!
  !! Public interface for FNV_1 hash function
  !!
  public :: FNV_1
  interface FNV_1
    module procedure FNV_1_shortInt
    module procedure FNV_1_longInt
  end interface

  !!
  !! Public interface for Knuth Hash
  !!
  public :: knuthHash

contains

  !!
  !! Multiplicative Hash outlined in D. Knuth "The Art of Computer Programming"
  !! Ratio was changed from the popular "golden ratio" to a value based on sphere
  !! maximum packing fraction in Euclidian geometry = pi/3/sqrt(2) [3180339487 for 32 bits]
  !! NOTE: Ratio was changed for fun -> no particular reason
  !!
  pure function knuthHash(key,m) result(hash)
    integer(shortInt), intent(in) :: key
    integer(shortInt), intent(in) :: m
    integer(shortInt)             :: hash
    integer(shortInt),parameter   :: prime = transfer(z'bd90211f',shortInt)
    integer(shortInt)             :: N
    character(100), parameter     :: Here='knuthHash (hashFunctions_func.f90)'

    !! Get bit size of the integer
    N = bit_size(key)

    !! Hash functions should be pure -> Change error checking
!    if(N /= 32) call fatalError(Here,'Implementation assumes 32 bit integer')
!    if(m < 1 )  call fatalError(Here,'At least one bit of hash needs to be kept')
!    if(m > N )  call fatalError(Here,'No more than 32 bits can be kept')

    ! Calculate prime * key modulo word
    hash = prime * key

    ! Keep m uppermost bits
    hash = shiftr(hash,N-m)

  end function knuthHash


  !!
  !!
  !!
  subroutine FNV_1_shortInt(key,hash)
    character(*),intent(in)         :: key
    integer(shortInt), intent(out)  :: hash
    integer(shortInt),parameter     :: FNV_prime  = 16777619_shortInt
    integer(shortInt),parameter     :: FNV_offset = transfer(z'811c9dc5',shortInt)
    integer(shortInt)               :: bajt, i
    character(100),parameter    :: Here ='FNV_1_shortInt (hashFunctions_func.f90)'

    if(storage_size(hash) /= 32) call fatalError(Here,'hash int is not 32bit')

    bajt = ichar(key(1:1),shortInt)
    hash = FNV_offset
    hash = hash * FNV_prime
    hash = ieor(hash,bajt)

    do i=2,len(key)
      bajt = ichar(key(i:i),shortInt)
      hash = hash * FNV_prime
      hash = ieor(hash,bajt)
    end do

  end subroutine FNV_1_shortInt


  !!
  !!
  !!
  subroutine FNV_1_longInt(key,hash)
    character(*),intent(in)      :: key
    integer(longInt),intent(out) :: hash
    integer(longInt),parameter   :: FNV_prime  = 1099511628211 _longInt
    integer(longInt),parameter   :: FNV_offset = transfer(z'cbf29ce484222325',longInt)
    integer(longInt)             :: bajt
    integer(shortInt)            :: i
    character(100),parameter   :: Here ='FNV_1_longInt (hashFunctions_func.f90)'

    if(storage_size(hash) /= 64) call fatalError(Here,'hash int is not 64bit')

    bajt = ichar(key(1:1),shortInt)
    hash = FNV_offset
    hash = hash * FNV_prime
    hash = ieor(hash,bajt)

    do i=2,len(key)
      bajt = ichar(key(i:i),shortInt)
      hash = hash * FNV_prime
      hash = ieor(hash,bajt)
    end do

  end subroutine FNV_1_longInt

    
end module hashFunctions_func
