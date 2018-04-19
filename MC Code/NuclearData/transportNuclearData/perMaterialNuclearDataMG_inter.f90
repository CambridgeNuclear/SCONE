module perMaterialNuclearDataMG_inter

  use numPrecision
  use genericProcedures,          only : fatalError
  use RNG_class,                  only : RNG
  use transportNuclearData_inter, only : transportNuclearData
  use particle_class,             only : particle

  use xsMacroSet_class,           only : xsMacroSet_ptr

  implicit none
  private

  type, public,extends(transportNuclearData),abstract :: perMaterialNuclearDataMG
    private

  contains
    ! Extend generic procedures to accept calls with energy value
    generic :: getTransXS    => getTransXS_G
    generic :: getMajorantXS => getMajorantXS_G
    generic :: getTotalMatXS => getTotalMatXS_G

    ! Adapters to transport procedures
    procedure :: getTransXS_p
    procedure :: getMajorantXS_p
    procedure :: getTotalMatXS_p

    ! Access to transport xs data by Energy Group
    procedure(getTransXS_G),deferred    :: getTransXS_G
    procedure(getMajorantXS_G),deferred :: getMajorantXS_G
    procedure(getTotalMatXS_G),deferred :: getTotalMatXS_G

    procedure(getMatMacroXS), deferred  :: getMatMacroXS
    procedure(releaseAt), deferred      :: releaseAt
    procedure(sampleMuGout), deferred   :: sampleMuGout

  end type perMaterialNuclearDataMG


  abstract interface
    !!
    !! Return transport XS of a given material.
    !! Transport XS is used to determine how far neutron flies. For transport-correced
    !! MG calculations it can be diffrent from material total XS. Usually is equal to material
    !! total XS.
    !!
    function getTransXS_G(self,G,matIdx) result (xs)
      import :: shortInt ,&
                defReal  ,&
                perMaterialNuclearDataMG
      class(perMaterialNuclearDataMG), intent(in) :: self
      integer(shortInt), intent(in)               :: G
      integer(shortInt), intent(in)               :: matIdx
      real(defReal)                               :: xs
    end function getTransXS_G

    !!
    !! Return majorant XS for all active materials
    !!
    function getMajorantXS_G(self,G) result(xs)
      import :: defReal  ,&
                shortInt ,&
                perMaterialNuclearDataMG
      class(perMaterialNuclearDataMG), intent(in) :: self
      integer(shortInt), intent(in)               :: G
      real(defReal)                               :: xs
    end function getMajorantXS_G

    !!
    !! Get total XS of a given material
    !!
    function getTotalMatXS_G(self,G,matIdx) result (xs)
      import :: shortInt ,&
                defReal  ,&
                perMaterialNuclearDataMG
      class(perMaterialNuclearDataMG), intent(in) :: self
      integer(shortInt), intent(in)               :: G
      integer(shortInt), intent(in)               :: matIdx
      real(defReal)                               :: xs
    end function getTotalMatXS_G


    !!
    !! Get set of material macroscopic xross-sections
    !!
    subroutine getMatMacroXS(self,macroXS,G,matIdx)
      import :: xsMacroSet_ptr, &
                defReal, &
                shortInt, &
                perMaterialNuclearDataMG
      class(perMaterialNuclearDataMG), intent(in) :: self
      type(xsMacroSet_ptr),intent(inout)          :: macroXS
      integer(shortInt),intent(in)                :: G
      integer(shortInt),intent(in)                :: matIdx
    end subroutine getMatMacroXS

    !!
    !! Obtain average emission for reaction MT
    !!
    function releaseAt(self,G_in,G_out,MT,matIdx) result(nu)
      import :: defReal, &
                shortInt, &
                perMaterialNuclearDataMG
      class(perMaterialNuclearDataMG), intent(in) :: self
      integer(shortInt), intent(in)               :: G_in
      integer(shortInt), intent(in)               :: G_out
      integer(shortInt), intent(in)               :: MT
      integer(shortInt), intent(in)               :: matIdx
      real(defReal)                               :: nu
    end function releaseAt

    !!
    !! Sample cosine of deflection angle mu and final energy group
    !!
    subroutine sampleMuGout(self,mu,G_out,G_in,rand,MT,matIdx)
      import :: RNG,&
                defReal, &
                shortInt, &
                perMaterialNuclearDataMG
      class(perMaterialNuclearDataMG), intent(in) :: self
      real(defReal), intent(out)                  :: mu
      integer(shortInt), intent(out)              :: G_out
      integer(shortInt), intent(in)               :: G_in
      class(RNG), intent(inout)                   :: rand
      integer(shortInt), intent(in)               :: MT
      integer(shortInt), intent(in)               :: matIdx

    end subroutine sampleMuGout

  end interface


contains


  !!
  !! getTransXS adapter to translate call with particle to a call with ENERGY GROUP
  !! Returns error if multigroup neutron is provided
  !!
  function getTransXS_p(self,p,matIdx) result (xs)
    class(perMaterialNuclearDataMG), intent(in) :: self
    class(particle), intent(in)                 :: p
    integer(shortInt), intent(in)               :: matIdx
    real(defReal)                               :: xs
    character(100), parameter            :: Here='getTransXS_p (perMaterialNuclearDataMG_inter.f90)'

    if (.not.p % isMG) then
      call fatalError(Here,'CE neutron given to MG nuclear Data')
    end if

    xs = self % getTransXS_G(p % G, matIdx)

  end function getTransXS_p

  !!
  !! getMajorantXS adapter to translate call with particle to a call with ENERGY GROUP
  !! Returns error if multigroup neutron is provided
  !!
  function getMajorantXS_p(self,p) result (xs)
    class(perMaterialNuclearDataMG), intent(in) :: self
    class(particle), intent(in)                 :: p
    real(defReal)                               :: xs
    character(100), parameter         :: Here='getMajorantXS_p (perMaterialNuclearDataMG_inter.f90)'

    if (.not.p % isMG) then
      call fatalError(Here,'CE neutron given to MG nuclear Data')
    end if

    xs = self % getMajorantXS_G(p % G)

  end function getMajorantXS_p

  !!
  !! getTotalMatXS adapter to translate call with particle to a call with ENERGY GROUP
  !! Returns error if multigroup neutron is provided
  !!
  function getTotalMatXS_p(self,p,matIdx) result (xs)
    class(perMaterialNuclearDataMG), intent(in) :: self
    class(particle), intent(in)                 :: p
    integer(shortInt), intent(in)               :: matIdx
    real(defReal)                               :: xs
    character(100), parameter         :: Here='getTotalMatXS_p (perMaterialNuclearDataMG_inter.f90)'

    if (.not.p % isMG) then
      call fatalError(Here,'CE neutron given to MG nuclear Data')
    end if

    xs = self % getTotalMatXS_G(p % G, matIdx)

  end function getTotalMatXS_p

    
end module perMaterialNuclearDataMG_inter
