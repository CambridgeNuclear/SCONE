module nuclearData_class

  use numPrecision
  use particle_class
  use RNG_class

  implicit none
  private

  type,abstract, public :: nuclearData
    private
  contains
    ! Change of state procedures
    procedure(setReal),deferred   :: setEnergy
    procedure(setInt),deferred    :: setMaterial
    procedure(setInt),deferred    :: setIsotope
    procedure(setInt),deferred    :: setReaction
    ! State access functions
    procedure(access),deferred    :: Energy
    procedure(access),deferred    :: Material
    procedure(access),deferred    :: Isotope
    procedure(access),deferred    :: Reaction
    ! Nuclear Data Requests
    procedure(access),deferred           :: majorantXS
    procedure(access),deferred           :: totalXS
    procedure(RealIntPointers),deferred  :: macroXS
    procedure(RealIntPointers),deferred  :: microXS
    procedure(access),deferred           :: emits    ! number of neutron emissions for current MT
    procedure(emit),deferred             :: emit     ! takes RNG and resturns energy and angle of emission

  end type nuclearData

  abstract interface

    subroutine setReal(self, R)
      import :: nuclearData, &
                defReal

      class(nuclearData), intent(inout) :: self
      real(kind=defReal), intent(in)    :: R
    end subroutine

    subroutine setInt(self, I)
      import :: nuclearData, &
                shortInt
      class(nuclearData), intent(inout)  :: self
      integer(kind=shortInt), intent(in) :: I
    end subroutine

    function access(self)
      import :: nuclearData
      class(nuclearData), intent(out)  :: self
    end function

    subroutine RealIntPointers(self, R, I)
      import :: nuclearData, &
                defReal, &
                shortInt
      class(nuclearData), intent(inout)           :: self
      real(kind=defReal),dimension(:),pointer     :: R
      integer(kind=shortInt),dimension(:),pointer :: I
    end subroutine

    subroutine emit(self, E, miu, R)
      import :: nuclearData, &
                defReal, &
                RNG
      class(nuclearData), intent(inout)         :: self
      class(RNG),intent(inout)                  :: R
      real(kind=defReal), intent(out)           :: E, miu
    end subroutine


  end interface
    
end module nuclearData_class
