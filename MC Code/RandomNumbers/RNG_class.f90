module RNG_class

  use numPrecision

 implicit none

    type, public :: rng
        private
        integer(longInt) :: g
        integer(longInt) :: c
        integer(longInt) :: m
        integer(longInt) :: bitMask
        real(defReal) :: norm       ! Normalisation factor to [0,1)
        integer(longInt) :: rngSeed
        ! Private data for a single history
        integer(longInt) :: rngCount  !rngNPS
        !common      /rngThread/ rngSeed, rngCount, rngNPS
        !!$OMP THREADPRIVATE( /rnThread/ )

    contains
        procedure :: init
        procedure :: get
        procedure :: skip
    end type rng

contains

  subroutine init(self, seed)
    implicit none
    class(rng), intent(inout) :: self
    integer(longInt), intent(in) :: seed
    self%rngSeed = seed
    self%g = 2806196910506780709_8
    self%m = ibset(0_8,63)
    self%norm = 2._8**(-63)
    self%bitMask = self%m - 1_8
    self%c = 1.0
  end subroutine init

  !
  ! Calculate the new seed for the RNG when skipping forward k steps
  ! Uses formula: Sk = G(k)S0 + C(k) mod 2**m
  !
  subroutine skip(self, k)
    implicit none
    class(rng) :: self
    integer(longInt) :: k          ! places to skip forward
    integer(longInt) :: Gk = 1     ! G**k mod M
    integer(longInt) :: Ck = 0     ! c*(g**k-1)/(g-1) mod M
    integer(longInt) :: h          ! sub for g
    integer(longInt) :: f          ! sub for c
    integer(longInt) :: mask

    h = self%g
    f = self%c
    mask = self%bitMask
    ! Can jump backwards or forwards due to periodicity of RNG
    if(k<0) then
      k = k + self%m
    end if

    k = iand(k + self%m,mask)
    do while( k > 0)
      if(iand(k,2_longInt)==1) then
        Gk = iand(Gk * h,mask)
        Ck = iand(Ck * h + f,mask)
      end if
      f = iand(f*(h+1),mask)
      h = iand(h*h,mask)
      k = ishft(k, -1)
    end do

    self%rngSeed = iand(Gk * self%rngSeed + Ck, mask)

  end subroutine skip

  function get(self) result(rand)
    implicit none
    class(rng) :: self
    integer(longInt) :: seed
    real(longInt) :: rand
    integer(longInt) :: mask

    mask = self%bitMask
    seed = self%rngSeed
    seed = iand(self%g*seed, mask)
    seed = iand(seed + self%c, mask)
    rand = seed*self%norm
    self%rngSeed = seed
    self%rngCount = self%rngCount + 1
    return
  end function get

end module rng_class
