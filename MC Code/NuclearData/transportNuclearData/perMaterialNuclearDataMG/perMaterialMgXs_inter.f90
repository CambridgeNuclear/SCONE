module perMaterialMgXs_inter

  use numPrecision
  use RNG_class,         only : RNG
  use dictionary_class,  only : dictionary
  use nuclearData_inter, only : nuclearData

  ! CrossSectionPackages
  use xsMacroSet_class, only : xsMacroSet

  implicit none
  private

  type, public,extends(nuclearData), abstract :: perMaterialMgXs
    private
  contains

    ! Procedures to obtain XSs
    procedure(getMatMacroXS), deferred   :: getMatMacroXS
    procedure(getTransXS), deferred      :: getTransXS
    procedure(getMajorantXS), deferred   :: getMajorantXS
    procedure(getTotalMatXS), deferred   :: getTotalMatXS

    ! Procedures to obtain emission data
    procedure(releaseAt), deferred       :: releaseAt
    procedure(sampleMuGout), deferred    :: sampleMuGout

    ! Procedures to access material information
   ! procedure(getMatIdx), deferred       :: getMatIdx
   ! procedure(getMatName), deferred      :: getMatName
    procedure(isFissileMat), deferred    :: isFissileMat



  end type perMaterialMgXs


  abstract interface

    !!
    !! Interface to obtain main XS of the material
    !!
    subroutine getMatMacroXS(self,macroXS,G,matIdx)
      import :: perMaterialMgXS,&
                xsMacroSet, &
                defReal, &
                shortInt
      class(perMaterialMgXS), intent(inout)        :: self
      type(xsMacroSet),pointer,intent(inout)       :: macroXS
      integer(shortInt),intent(in)                 :: G
      integer(shortInt),intent(in)                 :: matIdx
    end subroutine getMatMacroXS

    !!
    !! Get transport XS for a material
    !!
    function getTransXS(self,G,matIdx) result(xs)
      import :: perMaterialMgXS,&
                defReal, &
                shortInt
      class(perMaterialMgXS), intent(inout)  :: self
      integer(shortInt), intent(in)          :: G
      integer(shortInt), intent(in)          :: matIdx
      real(defReal)                          :: xs
    end function getTransXS

    !!
    !! Get majorant XS
    !!
    function getMajorantXS(self,G) result (xs)
      import :: perMaterialMgXS,&
                defReal, &
                shortInt
      class(perMaterialMgXS), intent(inout) :: self
      integer(shortInt), intent(in)         :: G
      real(defReal)                         :: xs
    end function getMajorantXS

    !!
    !! Get total XS of the material
    !!
    function getTotalMatXS(self,G,matIdx) result (xs)
      import :: perMaterialMgXS,&
                defReal, &
                shortInt
      class(perMaterialMgXS), intent(inout) :: self
      integer(shortInt), intent(in)         :: G
      integer(shortInt), intent(in)         :: matIdx
      real(defReal)                         :: xs
    end function getTotalMatXS


    !!
    !! Obtain average emission for reaction MT
    !!
    function releaseAt(self,G_in,G_out,MT,matIdx) result(nu)
      import :: perMaterialMgXS,&
                defReal, &
                shortInt
      class(perMaterialMgXS), intent(inout) :: self
      integer(shortInt), intent(in)         :: G_in
      integer(shortInt), intent(in)         :: G_out
      integer(shortInt), intent(in)         :: MT
      integer(shortInt), intent(in)         :: matIdx
      real(defReal)                         :: nu
    end function releaseAt

    !!
    !! Sample cosine of deflection angle mu and final energy group
    !!
    subroutine sampleMuGout(self,mu,G_out,G_in,rand,MT,matIdx)
      import :: perMaterialMgXS,&
                RNG,&
                defReal, &
                shortInt
      class(perMaterialMgXS), intent(inout) :: self
      real(defReal), intent(out)            :: mu
      integer(shortInt), intent(out)        :: G_out
      integer(shortInt), intent(in)         :: G_in
      class(RNG), intent(inout)             :: rand
      integer(shortInt), intent(in)         :: MT
      integer(shortInt), intent(in)         :: matIdx

    end subroutine sampleMuGout

    !!
    !! Get material index associated with material name
    !!
    function getMatIdx(self,matName) result(matIdx)
      import :: perMaterialMgXS,&
                shortInt
      class(perMaterialMgXS), intent(in)  :: self
      character(*), intent(in)            :: matName
      integer(shortInt)                   :: matIdx
    end function getMatIdx

    !!
    !! Get mataerial name associated with material index
    !!
    function getMatName(self,matIdx) result(matName)
      import :: perMaterialMgXS,&
                nameLen, &
                shortInt
      class(perMaterialMgXS), intent(in)  :: self
      integer(shortInt), intent(in)       :: matIdx
      character(nameLen)                  :: matName
    end function getMatName

    !!
    !! Returns .true. if material is fissile
    !!
    function isFissileMat(self,matIdx) result(isIt)
      import :: perMaterialMgXS,&
                defBool, &
                shortInt
      class(perMaterialMgXS), intent(in)  :: self
      integer(shortInt), intent(in)       :: matIdx
      logical(defBool)                    :: isIt
    end function

  end interface

contains

end module perMaterialMgXs_inter
