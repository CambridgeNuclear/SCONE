module perNuclideImplicitCaptureCE_class

  use numPrecision
  use endfConstants
  use genericProcedures,             only : fatalError, rotateVector
  use particle_class,                only : particle
  use particleDungeon_class,         only : particleDungeon
  use perNuclideNuclearDataCE_inter, only : perNuclideNuclearDataCE
  use perNuclideCollisionOpCE_class, only : perNuclideCollisionOpCE, &
                                            cutoffs_base  => cutoffs, &
                                            implicit_base => implicit

  ! Cross-section packages to interface with nuclear data
  use xsNucMacroSet_class,    only : xsNucMacroSet_ptr
  use xsMainSet_class,        only : xsMainSet_ptr

  implicit none
  private


  type, public,extends(perNuclideCollisionOpCE) :: perNuclideImplicitCaptureCE
  contains
    ! Overwride procedures to sample colision, perform implicit reaction treatment and
    ! perform cutoffs
    procedure  :: sampleCollision
    procedure  :: implicit
    procedure  :: cutoffs
  end type perNuclideImplicitCaptureCE

contains

  !!
  !! Implementation of the interface procedure.
  !! Overrides procedure in parent class
  !! Samples collision nuclide and a reaction channel
  !!
  subroutine sampleCollision(self,p,thisCycle,nextCycle)
    class(perNuclideImplicitCaptureCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle
    type(xsNucMacroSet_ptr)                       :: nucXSs
    real(defReal)                                 :: r
    character(100),parameter :: Here =' sampleCollision (perNuclideImplicitCaptureCE_class.f90)'

    ! Check if particle is multigroup
    if( p % isMG ) then
      call fatalError(Here,'MG neutron given to CE operator')
    end if

    ! Load collision energy
    self % E = p % E

    ! Select collision nuclide
    call self % xsData % getNucMacroXs(nucXSs, self % E, self % matIdx)

    r = self % pRNG % get()
    self % nucIdx = nucXSs % invert(r)

    ! Select Main reaction channel - always scattering
    self % MT = anyScatter

  end subroutine sampleCollision

  !!
  !! Implementation of the interface procedure.
  !! Is executed just after sampleCollision
  !! Overrides procedure in parent class
  !! Performs implicit fission and capture
  !!
  subroutine implicit(self,p,thisCycle,nextCycle)
    class(perNuclideImplicitCaptureCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle
    type(xsMainSet_ptr)                           :: nuclideXss
    real(defReal)                                 :: sig_tot, sig_abs

    ! Call parent class subroutine
    call implicit_base(self,p,thisCycle,nextCycle)

    ! Perform implicit capture
    call self % xsData % getMainNucXs(nuclideXss, self % E, self % nucIdx)

    sig_abs = nuclideXss % capture()+ nuclideXss % fission()
    sig_tot = nuclideXss % total()

    p % w = p % w * (ONE - (sig_abs/sig_tot) )

  end subroutine implicit

  !!
  !! Implementation of the interface procedure
  !! Overrides procedure in parent class
  !! Is called at the end
  !! Performs russian rulette
  !!
  subroutine cutoffs(self,p,thisCycle,nextCycle)
    class(perNuclideImplicitCaptureCE), intent(inout) :: self
    class(particle), intent(inout)                :: p
    class(particleDungeon),intent(inout)          :: thisCycle
    class(particleDungeon),intent(inout)          :: nextCycle
    real(defReal)                                 :: r1
    real(defReal), parameter              :: W_cut = 0.25_defReal
    real(defReal), parameter              :: W_new = 1.0_defReal

    ! Call base class procedure
    call cutoffs_base(self,p,thisCycle,nextCycle)

    ! Apply russian rullete
    if (p % w < W_Cut ) then
      r1 = self % pRNG % get ()
      if ( r1 > ONE - p % w / W_new ) then
        p % w = W_new

      else
        p % isDead = .true.

      end if
    end if

  end subroutine cutoffs

end module perNuclideImplicitCaptureCE_class
