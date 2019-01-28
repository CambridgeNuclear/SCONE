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
  !!
  !! Implementation uses a 32-bit hash. To avoid complier warning associated with integer
  !! overflow on multiplication a temporary integer of higher precision is used and
  !! then converted to shortInt by ignoring higher bits.
  !!
  !! Ratio was changed from the popular "golden ratio" to a value based on sphere
  !! maximum packing fraction in Euclidian geometry = pi/3/sqrt(2) [3180339487 for 32 bits]
  !! NOTE: Ratio was changed for fun -> no particular reason
  !!
  !!
  pure function knuthHash(key,m) result(hash)
    integer(shortInt), intent(in) :: key
    integer(shortInt), intent(in) :: m
    integer(shortInt)             :: hash
    integer(shortInt)             :: m_loc
    integer(shortInt),parameter   :: kin64 = selected_int_kind(18) ! Choose at least 60-bit integer
    integer(kin64)                :: hash_loc
    integer(kin64),parameter      :: prime = 3180339487_kin64

    ! Constrain m to range 1-32
    m_loc = max(1,m)
    m_loc = min(32,m_loc)

    ! Calculate prime * key modulo word
    hash_loc = modulo(prime * key, 4294967296_kin64)

    ! Keep m uppermost bits
    hash = transfer(shiftr(hash_loc,32-m),shortInt)

  end function knuthHash


  !!
  !!
  !!
  pure subroutine FNV_1_shortInt(key,hash)
    character(*),intent(in)         :: key
    integer(shortInt), intent(out)  :: hash
    integer(shortInt),parameter     :: FNV_prime  = 16777619_shortInt
    integer(shortInt),parameter     :: FNV_offset = transfer(z'811c9dc5',shortInt)
    integer(shortInt)               :: bajt, i
    character(100),parameter    :: Here ='FNV_1_shortInt (hashFunctions_func.f90)'

    !! Hash functions should be pure -> Change error checking
   ! if(storage_size(hash) /= 32) call fatalError(Here,'hash int is not 32bit')

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

    !! Hash functions should be pure -> Change error checking
    !if(storage_size(hash) /= 64) call fatalError(Here,'hash int is not 64bit')

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
