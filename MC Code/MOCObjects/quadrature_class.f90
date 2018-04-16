!!
!! Defines the quadrature for a 2D MOC solver
!! Contains azimuthal and polar angles and weight, and azimuthal widths
!! Polar quadrature may be Gauss-Legendre or Tabuchi-Yamamoto
!!
module quadrature_class

  use numPrecision
  use universalVariables
  use MOCVariables
  use genericProcedures

  implicit none
  private

  type, public                                 :: quadrature
    integer(shortInt)                          :: nAz        !! Number of azimuthal angles
    integer(shortInt)                          :: n_2        !! Half the number of azimuthal angles
    integer(shortInt)                          :: nPol       !! Number of polar angles
    real(defReal), dimension(:), allocatable   :: wAz        !! Azimuthal weights
    real(defReal), dimension(:), allocatable   :: phi        !! Azimuthal angles
    real(defReal), dimension(:), allocatable   :: wPol       !! Polar weights
    real(defReal), dimension(:), allocatable   :: sinTheta   !! Sines of polar angles
    real(defReal), dimension(:), allocatable   :: trackSep   !! Azimuthal track separation
    real(defReal), dimension(:,:), allocatable :: weights    !! Array of weights indexed by azimuthal and polar angle
  contains
    procedure          :: init
    procedure          :: getWeight
    procedure, private :: polarTY
    procedure, private :: polarGL
    procedure, private :: azimuthalCyclic
    procedure, private :: calculateWeights
  end type quadrature

contains

  !!
  !! Constructor for quadrature
  !!
  subroutine init(self, nAz, nPol, trackSep, xWidth, yWidth, polarType)
    class(quadrature), intent(inout) :: self
    integer(shortInt), intent(inout) :: nAz
    integer(shortInt), intent(in) :: nPol
    real(defReal), intent(in) :: trackSep
    real(defReal), intent(in) :: xWidth
    real(defReal), intent(in) :: yWidth
    integer(shortInt), intent(in), optional :: polarType
    integer(shortInt) :: rem

    ! Check nAz is a sensible value
    if (nAz < 1) then
      call fatalError('init, quadrature','Number of azimuthal angles must be greater than zero')
    end if

    ! Check whether the azimuthal angles are divisible by 4
    rem = modulo(nAz,4)
    if (rem .NE. 0) then
      call warning('init, quadrature','Number of azimuthal angles not divisible by 4: number has been adjusted')
      nAz = nAz + 4 - rem
    end if
    self % nAz = nAz
    self % n_2 = nAz/2

    ! Calculate effective azimuthal angles and widths
    call self % azimuthalCyclic(trackSep, xWidth, yWidth)

    ! Ensure nPol is a sensible value
    if (nPol < 1) then
      call fatalError('init, quadrature','Number of polar angles must be greater than zero')
    end if

    ! If polarType isn't supplied, use TY by default
    if((.not.present(polarType)).OR.(polarType == TY)) then
      call self % polarTY(nPol)
    else if (polarType == GL) then
      call self % polarGL(nPol)
    else
      call fatalError('init, quadrature','Invalid polar quadrature type provided')
    end if

    ! Calculate combined azimuthal and polar weights
    call self % calculateWeights()

  end subroutine init

  !!
  !! Return quadrature weight given polar and azimuthal index
  !! Assumes all azimuthal indices are assigned with half symmetry
  !!
  function getWeight(self, idxAz, idxPol)result(w)
    class(quadrature), intent(in) :: self
    integer(shortInt), intent(in) :: idxAz
    integer(shortInt), intent(in) :: idxPol
    real(defReal), intent(out) :: w
    w = self % weights(idxAz, idxPol)
  end function getWeight

  !!
  !! Initialise azimuthal quadrature for cyclic tracking
  !! Calculates only for angles between 0 and Pi due to symmetry
  !!
  subroutine azimuthalCyclic(self, trackSep, xWidth, yWidth)
    class(quadrature), intent(inout) :: self
    real(defReal), intent(in) :: trackSep
    real(defReal), intent(in) :: xWidth
    real(defReal), intent(in) :: yWidth
    integer(shortInt) :: n_2
    integer(shortInt) :: i
    integer(shortInt) :: nx, ny
    real(defReal) :: phi

    n_2 = self % n_2
    allocate(self % phi(n_2))
    allocate(self % wAz(n_2))
    allocate(self % trackSep(n_2))

    ! Calculate effective angles and track separations
    do i=1,n_2/2
      phi = 2. * PI * (i - 0.5)/nAz
      nx = floor(yWidth*abs(cos(phi))/trackSep) + 1
      ny = floor(xWidth*abs(sin(phi))/trackSep) + 1
      self % phi(i) = atan(yWidth*nx/xWidth/ny)
      self % trackSep(i) = xWidth * sin(self % phi(i)) / nx
    end do

    ! Calculate quadrature weights
    do i=1,n_2
      if (i==1) then
        self % wAz(i) = (self % phi(i+1) - self % phi(i))/2. + self % phi(i)
      else if (i==nAz) then
        self % wAz(i) = PI - self % phi(i) + (self % phi(i) - self % phi(i-1))/2.
      else
        self % wAz(i) = (self % phi(i+1) - self % phi(i-1))
      end if
    end do

  end subroutine azimuthalCyclic

  !!
  !! Initialise Tabuchi-Yamamoto polar quadrature
  !!
  subroutine polarTY(self,nPol)
    class(quadrature), intent(inout) :: self
    integer(shortInt), intent(in) :: nPol

    if (nPol < 4) then
      allocate(self % wPol(nPol))
      allocate(self % sinTheta(nPol))
      if (nPol == 1) then
        self % wPol(1) = 1._defReal
        self % sinTheta(1) = 0.798184_defReal
      else if (nPol == 2) then
        self % wPol(1) = 0.212854_defReal
        self % wPol(2) = 0.787146_defReal
        self % sinTheta(1) = 0.363900_defReal
        self % sinTheta(2) = 0.899900_defReal
      else if (nPol == 3) then
        self % wPol(1) = 0.046233_defReal
        self % wPol(2) = 0.283619_defReal
        self % wPol(3) = 0.670148_defReal
        self % sinTheta(1) = 0.166648_defReal
        self % sinTheta(2) = 0.537707_defReal
        self % sinTheta(3) = 0.932954_defReal
      end if
    else
      call fatalError('polarTY, quadrature','TY quadrature can only use up to 3 polar angles')
    end if
  end subroutine polarTY

  !!
  !! Initialise Gauss-Legendre polar quadrature
  !!
  subroutine polarGL(self,nPol)
    class(quadrature), intent(inout) :: self
    integer(shortInt), intent(in) :: nPol
    real(defReal), dimension(:), allocatable :: cosTheta

    if (nPol < 7) then
      allocate(cosTheta(nPol))
      allocate(self % wPol(nPol))
      allocate(self % sinTheta(nPol))
      if (nPol == 1) then
        self % wPol(1) = 1._defReal
        cosTheta(1) = 0.5573502691
      else if (nPol == 2) then
        self % wPol(1) = 0.3478548451_defReal
        self % wPol(2) = 0.6521451549_defReal
        cosTheta(1) = 0.8611363115_defReal
        cosTheta(2) = 0.3399810435_defReal
      else if (nPol == 3) then
        self % wPol(1) = 0.1713244924_defReal
        self % wPol(2) = 0.3607615730_defReal
        self % wPol(3) = 0.4679139346_defReal
        cosTheta(1) = 0.9324695142_defReal
        cosTheta(2) = 0.6612093864_defReal
        cosTheta(3) = 0.2386191860_defReal
      else if (nPol == 4) then
        self % wPol(1) = 0.1012285363_defReal
        self % wPol(2) = 0.2223810344_defReal
        self % wPol(3) = 0.3137066459_defReal
        self % wPol(4) = 0.3626837834_defReal
        cosTheta(1) = 0.9602898564_defReal
        cosTheta(2) = 0.7966664774_defReal
        cosTheta(3) = 0.5255324099_defReal
        cosTheta(4) = 0.1834346424_defReal
      else if (nPol == 5) then
        self % wPol(1) = 0.0666713443_defReal
        self % wPol(2) = 0.1494513492_defReal
        self % wPol(3) = 0.2190863625_defReal
        self % wPol(4) = 0.2692667193_defReal
        self % wPol(5) = 0.2955242247_defReal
        cosTheta(1) = 0.9739065285_defReal
        cosTheta(2) = 0.8650633666_defReal
        cosTheta(3) = 0.6794095682_defReal
        cosTheta(4) = 0.4333953941_defReal
        cosTheta(5) = 0.1488743387_defReal
      else if (nPol == 6) then
        self % wPol(1) = 0.0471753364_defReal
        self % wPol(2) = 0.1069393260_defReal
        self % wPol(3) = 0.1600783286_defReal
        self % wPol(4) = 0.2031674267_defReal
        self % wPol(5) = 0.2334925365_defReal
        self % wPol(6) = 0.2491470458_defReal
        cosTheta(1) = 0.9815606342_defReal
        cosTheta(2) = 0.9041172563_defReal
        cosTheta(3) = 0.7699026741_defReal
        cosTheta(4) = 0.5873179542_defReal
        cosTheta(5) = 0.3678314989_defReal
        cosTheta(6) = 0.1252334085_defReal
      end if
      do i=1,nPol
        self % sinTheta(i) = sqrt(1-cosTheta(i)*cosTheta(i))
      end do
    else
      call fatalError('polarGL, quadrature','GL quadrature can only use up to 6 polar angles')
    end if
  end subroutine polarGL

  !!
  !! Calculates weights used in 2D MOC transport sweep
  !! Combines azimuthal weight, polar weight, track spacing and polar sine
  !!
  subroutine calculateWeights(self)
    class(quadrature), intent(inout) :: self
    integer(shortInt) :: i,j
    integer(shortInt) :: n_2
    integer(shortInt) :: nPol

    n_2 = self % n_2
    nPol = self % nPol
    allocate(self % weights(n_2,nPol))

    do i=1,n_2
      do j=1,nPol
        self % weights(i,j) = self % wAz(i) % self % trackSep(i)
        self % weights(i,j) = 2.0 * self % wPol(j) * self % sinTheta(j) * self % weights(i,j)
      end do
    end do
  end subroutine calculateWeights

end module quadrature_class
