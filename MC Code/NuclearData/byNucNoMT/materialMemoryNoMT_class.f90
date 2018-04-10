module materialMemoryNoMT_class

  use numPrecision
  use genericProcedures,       only : fatalError
  use xsMacroSet_class,        only : xsMacroSet
  use materialDataNoMT_class,  only : materialDataNoMT
  use nuclideMemoryNoMT_class, only : nuclideMemoryNoMT
  use matNucCDF_class,         only : matNucCDF


  implicit none
  private

  !!
  !! Object to store material macroscopic XSs.
  !! Remembers energy so no recalculation of XS is necessary if the energy has not changed.
  !! Interface consists of:
  !!
  !! Function   :: getTotal(E) -> returns total XS for energy E
  !! Subroutine :: setTo(E)    -> calculates ALL macro XS for energy E and stores them
  !! Component  :: XS          -> "xsMacroSet" type which stores all Macro XS
  !! Component  :: nucCDF      -> CDF object to select collision nuclide
  !!
  !! WARNING:
  !! Currently it is not rebust. It is possible to access XS without interpolating them first
  !! with setTo(E) subroutine!
  !!
  type, public :: materialMemoryNoMT
    !private
    real(defReal)                  :: E       = -1.0    !! Current Energy of XS
    logical(defBool)               :: isCalc  = .false. !! TRUE if XSs other then TOTAL are calculated
    type(xsMacroSet),public        :: XS                !! Package of macroscopic XS
    type(matNucCDF),public         :: nucCDF            !! CDF to select collision nuclide
    type(materialDataNoMT),pointer :: data    => null() !! Pointer to data for this material

    type(nuclideMemoryNoMT),dimension(:),pointer :: nucShelf => null()  !! Pointer to microscopic xss for all nuclides

  contains
    procedure :: init
    procedure :: getTotal
    procedure :: setTotalToEnergy
    procedure :: setEnergy
    procedure :: checkZZIds

    procedure, private :: calculateTotal
    procedure, private :: calculateTail
    procedure, private :: calculateAll

  end type materialMemoryNoMT

contains

  !!
  !! Initialise material memory by providing it with a pointer to nuclide shelf and
  !! pointer to corresponding material in the data block.
  !!
  subroutine init(self,dataPtr,nucPtr)
    class(materialMemoryNoMT), intent(inout)                :: self
    type(materialDataNoMT), pointer, intent(in)             :: dataPtr
    type(nuclideMemoryNoMT),dimension(:),pointer,intent(in) :: nucPtr
    character(100), parameter                 :: Here ='init (materialMemoryNoMT_class.f90)'

    if(.not.associated(dataPtr)) call fatalError(Here,'dataPtr is null')
    if(.not.associated(nucPtr) ) call fatalError(Here,'nucPtr is null')

    self % data         => dataPtr
    self % nucShelf     => nucPtr

    ! Call subroutine to validate consistancy between nuclide ZZids and indexes
    call self % checkZZIds()

    call self % nucCDF % init(self % data % nucIdx)

  end subroutine


  !!
  !! Return material total XS for energy E.
  !! If the energy has changed calculate new TOTAL XS and set flag for other XSs to .false.
  !!
  function getTotal(self,E) result (totalXS)
    class(materialMemoryNoMT), intent(inout) :: self
    real(defReal), intent(in)                :: E
    real(defReal)                            :: totalXS

!    if ( self % E /= E ) then
!      call self % calculateTotal(E)
!      self % E = E
!      self % isCalc = .false.
!
!    end if
    call self % setTotalToEnergy(E)

    totalXS = self % XS % totalXS

  end function getTotal


  !!
  !! Calculate total macroscopic cross sections for energy E
  !! Do nothing if E is the same as privious value
  !! Change state flags.
  !!
  subroutine setTotalToEnergy(self,E)
    class(materialMemoryNoMT), intent(inout) :: self
    real(defReal), intent(in)                :: E

    if ( self % E /= E ) then
      call self % calculateTotal(E)
      self % E = E
      self % isCalc = .false.

    end if
  end subroutine setTotalToEnergy



  !!
  !! Calculate all macroscopic XS for energy E
  !!
  subroutine setEnergy(self,E)
    class(materialMemoryNoMT), intent(inout) :: self
    real(defReal), intent(in)                :: E
    logical(defBool)                         :: tailIsNotInter
    logical(defBool)                         :: sameEnergy

    tailIsNotInter = .not.( self % isCalc)
    sameEnergy     = (self % E == E)

    if (sameEnergy .and. tailIsNotInter) then
      call self % calculateTail(E)

    elseif(.not. sameEnergy) then
      call self % calculateAll(E)

    end if
    ! If energy has not changed and evrything is calculated
    ! DO NOTHING

  end subroutine setEnergy

  !!
  !! Calculate and store total XS for energy E
  !!
  subroutine calculateTotal(self,E)
    class(materialMemoryNoMT), intent(inout):: self
    real(defReal), intent(in)               :: E
    real(defReal)                           :: tempXS
    real(defReal)                           :: nucDen
    real(defReal)                           :: xsMicro
    real(defReal)                           :: xsMacro
    integer(shortInt)                       :: nucIdx
    integer(shortInt)                       :: i

    tempXS = 0.0

    do i=1,self % data % numNuc
      ! Get nuclide index, density and total micro xs
      nucIdx  = self % data % nucIdx(i)
      nucDen  = self % data % nucDens(i)
      xsMicro = self % nucShelf(nucIdx) % getTotal(E)
      xsMacro = xsMicro * nucDen

      ! Add i-th nuclide contibution to total XS
      tempXS = tempXS + xsMacro

      ! Store nuclide contribution in CDF
      self % nucCDF % nucTotalXS(i) = xsMacro
    end do

    self % XS % totalXS        = tempXS
    self % nucCDF % matTotalXS = tempXS

  end subroutine calculateTotal

  !!
  !! Calculate Tail ( all macroscopic XS EXCEPT TOTAL ) for energy E
  !! tempMacroXS is orgenised as follows:
  !! tempMacroXS(1) -> scattering
  !! tempMacroXS(2) -> capture
  !! tempMacroXS(3) -> fission
  !! xsMicro is orgenised in the same way.
  !!
  subroutine calculateTail(self,E)
    class(materialMemoryNoMT), intent(inout) :: self
    real(defReal)                            :: E
    real(defReal),dimension(3)               :: tempMacroXS
    real(defReal),dimension(3)               :: xsMicro
    real(defReal)                            :: nucDen
    integer(shortInt)                        :: i, nucIdx


    tempMacroXS = 0.0

    do i=1,self % data % numNuc
      ! Find nuclide index and energy
      nucIdx = self % data % nucIdx(i)
      nucDen = self % data % nucDens(i)

      ! Load micoscopic xss into local vector
      call self % nucShelf(nucIdx) % setEnergy(E) ! <<<<< IMPORTANT NOT TO FORGET IT (SET ENERGY)
      xsMicro(1) = self % nucShelf(nucIdx) % xs % scatter
      xsMicro(2) = self % nucShelf(nucIdx) % xs % capture
      xsMicro(3) = self % nucShelf(nucIdx) % xs % fission

      ! Increase Material macroscopic XSs by the nuclide macroscopic XSs
      tempMacroXS = tempMacroXS + xsMicro * nucDen

    end do

    ! Load XS into Macro XS storage
    self % XS % scatterXS = tempMacroXS(1)
    self % XS % captureXS = tempMacroXS(2)
    self % XS % fissionXS = tempMacroXS(3)

  end subroutine calculateTail

  !!
  !! Calculate all macroscopic XSs for energy E
  !! tempMacroXS is orgenised as follows:
  !! tempMacroXS(1) -> total
  !! tempMacroXS(2) -> scattering
  !! tempMacroXS(3) -> capture
  !! tempMacroXS(4) -> fission
  !! xsMicro is orgenised in the same way.
  !!
  subroutine calculateAll(self,E)
  class(materialMemoryNoMT), intent(inout) :: self
    real(defReal)                            :: E
    real(defReal),dimension(4)               :: tempMacroXS
    real(defReal),dimension(4)               :: xsMicro
    real(defReal)                            :: nucDen
    integer(shortInt)                        :: i, nucIdx

    tempMacroXS = 0.0

    do i=1,self % data % numNuc
      ! Find nuclide index and energy
      nucIdx = self % data % nucIdx(i)
      nucDen = self % data % nucDens(i)

      ! Load micoscopic xss into local vector
      call self % nucShelf(nucIdx) % setEnergy(E) ! <<<<< IMPORTANT NOT TO FORGET IT (SET ENERGY)
      xsMicro(1) = self % nucShelf(nucIdx) % xs % total
      xsMicro(2) = self % nucShelf(nucIdx) % xs % scatter
      xsMicro(3) = self % nucShelf(nucIdx) % xs % capture
      xsMicro(4) = self % nucShelf(nucIdx) % xs % fission

      ! Increase Material macroscopic XSs by the nuclide macroscopic XSs
      tempMacroXS = tempMacroXS + xsMicro * nucDen

    end do

    ! Load XS into Macro XS storage
    self % XS % totalXS   = tempMacroXS(1)
    self % XS % scatterXS = tempMacroXS(2)
    self % XS % captureXS = tempMacroXS(3)
    self % XS % fissionXS = tempMacroXS(4)

  end subroutine calculateAll


  !!
  !! Subroutine to check conistancy between ZZids in material definition and coresponding
  !! ZZids in nuclide data
  !!
  subroutine checkZZids(self)
    class(materialMemoryNoMT), intent(in) :: self
    logical(defBool)                      :: isNotInitialised
    character(zzIdLen)                    :: zzFromHere, zzFromNuc
    integer(shortInt)                     :: nucIdx, i
    character(100), parameter             :: Here='checkZZids (materialMemoryNoMT_class.f90)'

    isNotInitialised = .not.(associated(self % data) .and. associated(self % nucShelf))

    if (isNotInitialised) then
      call fatalError(Here,'materialMemoryNoMT is not initialised')

    end if

    ! Loop over all nuclides in material
    do i=1,self % data % numNuc
      ! Get nuclide index
      nucIdx = self % data % nucIdx(i)

      ! Obtain ZZids from here and from nuclide Shelf
      ! adjustl is necessary to remove leading blanks
      zzFromHere = adjustl(self % data % nucNames(i))
      zzFromNuc  = adjustl(self % nucShelf(nucIdx) % getZZId())

      if (zzFromHere /= zzFromNuc) then
        call fatalError(Here,'Inconsistent ZZids Material: '//zzFromHere//' Nuclide: '//zzFromNuc)

      end if
    end do

  end subroutine checkZZIds


end module materialMemoryNoMT_class
