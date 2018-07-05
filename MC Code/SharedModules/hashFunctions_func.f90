module hashFunctions_func

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  !!
  !! Local constants
  !!
  ! Kind of at least 60- bit integer. Assume implementation is 64bit becouse.
  integer(shortInt), parameter  :: int64 = selected_int_kind(18)

  integer(int64),parameter      :: type_int64    = 0  ! An example 64bit int
  integer(shortInt),parameter   :: type_shortInt = 0  ! An example shortInt
  integer(longInt),parameter    :: type_longInt  = 0  ! An example longInt

  integer(shortInt),parameter   :: int64_size    = storage_size(type_int64)
  integer(shortInt),parameter   :: shortInt_size = storage_size(type_shortInt)
  integer(shortInt),parameter   :: longInt_size  = storage_size(type_longInt)

  ! Calculate how many long & short integers fits into 64 bit
  integer(shortInt),parameter   :: longIn64  = int64_size/longInt_size
  integer(shortInt),parameter   :: shortIn64 = int64_size/shortInt_size


  public :: FNV_1

  !!
  !! Unsigned 64-bit integer
  !!
  type, public :: unSig64Int
    integer(int64) :: value
  contains
    generic :: assignment(=)  => assign_shortInt, assign_longInt
    generic :: operator(+)    => add_unSig64Int
    generic :: operator(*)    => mul_unSig64Int

    procedure,private :: assign_shortInt
    procedure,private :: assign_longInt
    procedure,private :: add_unSig64Int
    procedure,private :: mul_unSig64Int
  end type unSig64Int

contains

!  function Knuth_Hash(intKey) result(hash)
!    integer(shortInt)
!
!  end function Knuth_hash


  function FNV_1(key) result(hash)
    character(*),intent(in) :: key
    character(len(key))     :: lKey
    integer(int64)          :: hash
    integer(int64),parameter  :: FNV_prime  = transfer(z'100000001b3',type_int64)
    integer(int64),parameter  :: FNV_offset = transfer(z'cbf29ce484222325',type_int64)
    integer(int64)           :: byt
    !type(unSig64Int)        :: hash_temp
    !type(unSig64Int)        :: FNV_prime
    !type(unSig64Int)        :: byt
    integer(shortInt)       :: N, i

    byt = ichar(key(1:1),int64)
    hash = FNV_offset
    hash = hash * FNV_prime
    hash = ieor(hash,byt)

    do i=2,len(key)
      byt  = ichar(key(i:i),int64)
      hash = hash * FNV_prime
      hash = ieor(hash,byt)
    end do

   ! do i=1,N


!    ! Made local kopy of key and store its length
!    lKey = key
!    N = len(lKey)
!
!    ! Load FNV prime number
!    FNV_prime  = transfer(z'100000001b3',type_int64)
!
!    ! Load offset basis
!    hash_temp  = transfer(z'cbf29ce484222325',type_int64)
!
!    do i=1,N
!      byt = ichar(lKey(i:i),shortInt)
!      hash_temp = hash_temp * FNV_prime
!      hash_temp = ieor(hash_temp % value, byt % value)
!    end do
!
!    hash = hash_temp % value

  end function FNV_1





  !!
  !! Assign unSig64 with value from shortInt
  !!
  subroutine assign_shortInt(LHS,RHS)
    class(unSig64Int), intent(out)         :: LHS
    integer(shortInt), intent(in)          :: RHS
    integer(shortInt),dimension(shortIn64) :: RHS_64
    character(100),parameter :: Here ='assign_shortInt (hashFunctions_func.f90)'

    if( incorrect_bit_assumptions()) then
      call fatalError(Here,'Invalid bit assumptions verify them on your system.')

    else
      ! Read provided value into rightmost elementy of the array
      ! Thus obtain representation of input value in 64 bit spread across multiple elements of array
      RHS_64 = 0
      RHS_64(1) = RHS

      ! Transfer bit pattern
      LHS % value = transfer(RHS_64, LHS % value)

    end if
  end subroutine assign_shortInt

  !!
  !! Assign unSig64 with value from longInt
  !!
  subroutine assign_longInt(LHS,RHS)
    class(unSig64Int), intent(out)        :: LHS
    integer(longInt), intent(in)          :: RHS
    integer(longInt),dimension(longIn64)  :: RHS_64
    character(100),parameter :: Here ='assign_longInt (hashFunctions_func.f90)'

    if( incorrect_bit_assumptions()) then
      call fatalError(Here,'Invalid bit assumptions verify them on your system.')

    else
      ! Read provided value into rightmost elementy of the array
      ! Thus obtain representation of input value in 64 bit spread across multiple elements of array
      RHS_64 = 0
      RHS_64(1) = RHS

      ! Transfer bit pattern
      LHS % value = transfer(RHS_64, LHS % value)

    end if
  end subroutine assign_longInt


  !!
  !! Function that checks assumptions about the integer representations
  !! This is becouse Fortran Standard does not specify bit representation and
  !! avalible kinds. We need to make shure our assumptions are correct
  !!
  function incorrect_bit_assumptions() result(incorrect)
    logical(defBool) :: incorrect

    ! Verify that int64 is 64-bit
    incorrect = (int64_size /= 64)

    ! Verify that natural number of short & long ints fits into int64
    incorrect = incorrect .or. (modulo(int64_size,longInt_size)  /= 0)
    incorrect = incorrect .or. (modulo(int64_size,shortInt_size) /= 0)

  end function incorrect_bit_assumptions

  !!
  !! Unsigned addition that happly overwrites sign bit
  !!
  function add_unSig64Int(LHS,RHS) result(res)
    class(unSig64Int), intent(in) :: LHS
    type(unSig64Int), intent(in)  :: RHS
    type(unSig64Int)              :: res
    type(unSig64Int)              :: a,c

    a   = RHS
    res = LHS
    do
      ! Construct carry vector
      c   = iand(res % value, a % value)

      ! Add numbers
      res = ieor(res % value, a % value)

      ! Shift carry to higher digits
      c % value = shiftl(c % value,1)

      ! Set to add carries in next iteration
      a = c

      if (a % value == 0) return
    end do

  end function add_unSig64Int

  !!
  !! Unsigned multiplication that happly overwrites sign bit
  !!
  function mul_unSig64Int(LHS,RHS) result(res)
    class(unSig64Int), intent(in) :: LHS
    type(unSig64Int), intent(in)  :: RHS
    type(unSig64Int)              :: res
    type(unSig64Int)              :: a,b

    res = LHS % value * RHS % value
    return
    ! Assign a to be the larger of RHS and LHS
    if ( RHS % value > LHS % value) then
      a   = RHS
      b   = LHS
    else
      a  = LHS
      b  = RHS
    end if

    res = 0

    ! Perform multiplication
    do
      if ( iand(b % value, int(1,int64)) /= 0 ) then
        res = res + a
      end if
      a = shiftl(a % value,1)
      b = shiftr(b % value,1)

      if (b % value == 0) return
    end do

!    if ((res % value) /= (LHS % value * RHS % value)) then
!      print *, res,
!    end if

  end function mul_unSig64Int

    
end module hashFunctions_func
