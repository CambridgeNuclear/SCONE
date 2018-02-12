module ByIsoNoMT_class_module

  use numPrecision
  use ByIsoNoMT_Data_class , only : ByIsoNoMT_Data
  use aceNoMT_class,         only : aceNoMT

  implicit none
  private

  integer(shortInt), parameter :: maxMT=4  ! Maximum Number of reaction channels


  type, private :: isotopeMemory
    private
    ! Stored
    integer(shortInt)                          :: nucIdx  = -1
    integer(shortInt)                          :: idx     = -1
    real(defReal)                              :: E       = -1.0
    integer(shortInt)                          :: nReact  = -1
    real(defReal),dimension(:),pointer         :: xs
    integer(shortInt),dimension(:),pointer     :: MT
  !  real(defReal),dimension(maxMT),target      :: xsSpace = -1.0
  !  integer(shortInt),dimension(maxMT),target  :: MTSpce  = -1
    type(aceNoMT),pointer                      :: dataPtr => null()
  contains
  !  procedure :: attach  => attach_isoMemory
  !  procedure :: xsTot   => xsTot_isoMemory
  !  procedure :: xss     => xss_isoMemory
  end type isotopeMemory



!  type, private :: materialMemory
!    ! Stored
!    real(defReal)                              :: E
!    real(defReal)                              :: kT   ! Temperature in MeV
!    real(defReal),dimension(maxIso),target     :: isoXS
!    integer(shortInt),dimension(maxIso),target :: isoMap
!    real(defReal),dimension(maxIso), target    :: isoDen
!  contains
!
!  end type



  type, public :: ByIsoNoMT_class
    private
    ! Material Data Handle State Variabels
    real(defReal)            ::  E
    integer(shortInt)        ::  material
    integer(shortInt)        ::  isotope
    integer(shortInt)        ::  MT

    ! Material Data Handle storage
    type(isotopeMemory), dimension(:), allocatable    :: isoShelf
   ! type(materialMemory),dimension(:), allocatable    :: matShelf
    real(defReal)                                     :: majorant
    type(ByIsoNoMT_Data),pointer                      :: dataBlock => null()

  contains
  !  procedure :: initDebug
  !  procedure :: init
    ! Get majorant XS
    ! Get material XS
    ! Get pointer to isotope total XS and pointer to isotope index map
    ! Get pointer to isotope's micro xs and MT map
    ! Get neutron emission given nuclide index and MT
    ! Sample deflection mu and secondary energy with CF flag given random number generator
  end type ByIsoNoMT_class

contains

  !!
  !! Subroutine to attach
  !!
!  subroutine attach_isoMemory(self,isoIdx,acePtr)
!    class(isotopeMemory), intent(inout) :: self
!    integer(shortInt), intent(in)       :: isoIdx
!    type(aceNoMT),pointer,intent(in)    :: acePtr
!
!    ! Clean data if alrady attached
!    if(associated(self % xs))      deallocate (self % xs)
!    if(associated(self % MT))      self % MT => null()
!    if(associated(self % dataPtr)) self % dataPtr => null()
!
!    ! Read basic data
!    self % nucIdx  =  isoIdx
!    self % dataPtr => acePtr
!    self % nReact  =  acePtr % nReact
!
!    ! Pointe internal pointers to approperiate part of space
!  !  self % xs = self % xsSpace(1:self % nReact)
!!    self % MT = self % MTSPace(1:self % nReact)
!
!    ! Load MT numbers
!    self % MT = acePtr % xsMT
!
!
!
!  end subroutine attach_isoMemory


!  subroutine initDebug(self,dataBlock)
!    class(ByIsoNoMT), intent(inout)           ::  self
!    type(ByIsoNoMT_Data), pointer, intent(in) :: dataBlock
!
! !   self % dataBlock => dataBlock
!    call self % init()
!  end subroutine initDebug


!  subroutine init(self)
!    class(ByIsoNoMT), intent(inout) :: self
!    ! Create
!    ! Here will be
!
!
!  end subroutine

end module ByIsoNoMT_class_module
