module transportNuclearData_inter

  use numPrecision
  use nuclearData_inter, only : nuclearData
  use particle_class,    only : particle

  implicit none
  private


  !!
  !! Interface of nuclear Data Module to be seen by transport routines
  !!
  !! This interface allows to use a single Transport Operator for both MG and CE data.
  !!
  type, public,extends(nuclearData) ,abstract :: transportNuclearData
    private
  contains
    ! Generic
    generic :: getTransXS      => getTransXS_p
    generic :: getMajorantXS   => getMajorantXS_p
    generic :: getTotalMatXS   => getTotalMatXS_p
    ! XS Access
    procedure(getTransXS_p), deferred       :: getTransXS_p
    procedure(getMajorantXS_p), deferred    :: getMajorantXS_p
    procedure(getTotalMatXS_p), deferred    :: getTotalMatXS_p

    procedure(isFissileMat), deferred       :: isFissileMat
    procedure(initFissionSite), deferred    :: initFissionSite
    procedure(setActiveMaterials), deferred :: setActiveMaterials

  end type transportNuclearData

  abstract interface
    !!
    !! Return transport XS of a given material.
    !! Transport XS is used to determine how far neutron flies. For transport-correced
    !! MG calculations it can be diffrent from material total XS. Usually is equal to material
    !! total XS.
    !!
    function getTransXS_p(self,p,matIdx) result (xs)
      import :: shortInt ,&
                defReal  ,&
                particle ,&
                transportNuclearData
      class(transportNuclearData), intent(inout) :: self
      class(particle), intent(in)                :: p
      integer(shortInt), intent(in)              :: matIdx
      real(defReal)                              :: xs
    end function getTransXS_p

    !!
    !! Return majorant XS for all active materials
    !!
    function getMajorantXS_p(self,p) result(xs)
      import :: defReal  ,&
                particle ,&
                transportNuclearData
      class(transportNuclearData), intent(inout) :: self
      class(particle), intent(in)                :: p
      real(defReal)                              :: xs
    end function getMajorantXS_p

    !!
    !! Get total XS of a given material
    !!
    function getTotalMatXS_p(self,p,matIdx) result (xs)
      import :: shortInt ,&
                defReal  ,&
                particle ,&
                transportNuclearData
      class(transportNuclearData), intent(inout) :: self
      class(particle), intent(in)                :: p
      integer(shortInt), intent(in)              :: matIdx
      real(defReal)                              :: xs
    end function getTotalMatXS_p

    !!
    !! Returns .true. if material is fissile.
    !! It is fissile if it contains fissile isotopes or MG fission data
    !!
    function isFissileMat(self,matIdx) result(isIt)
      import :: transportNuclearData,&
                defBool, &
                shortInt
      class(transportNuclearData), intent(in) :: self
      integer(shortInt), intent(in)           :: matIdx
      logical(defBool)                        :: isIt

    end function isFissileMat

    !!
    !! Function to generate a fission site from a fissile material.
    !! Necassary in initialisation of an eigenvalue calculation:
    !!
    !! If input particle isMG flag is .true. returns MG neutron
    !! Uses input particle RNG to sample fission site
    !!
    !! If provided material is fissile:
    !!   nuclide nuclear data:
    !!     For multiple fissile nuclides, selection of nuclide to generate a fission
    !!     neutron is determined by a given implementation. Each impelemntation needs
    !!     to be explicit about how it is done.
    !!
    !!   material nuclear data:
    !!     Only one set of data is avalible. No ambiguity.
    !!
    !! If provided material is not fissile:
    !!  Changes isAlive flag of a partice to .false.
    !!
    subroutine initFissionSite(self,p)
      import :: transportNuclearData, &
                particle
      class(transportNuclearData), intent(in) :: self
      class(particle), intent(inout)          :: p

    end subroutine initFissionSite


    !!
    !! Set all materials that are present in geometry
    !! Active materials will be included in evaluation of majorant XS
    !!
    !! Provide an array of all active material indexes
    !! Order is not significant (assume there is none)
    !!
    subroutine setActiveMaterials(self,matIdxList)
      import :: transportNuclearData, &
                shortInt
      class(transportNuclearData), intent(inout) :: self
      integer(shortInt),dimension(:), intent(in) :: matIdxList

    end subroutine

  end interface
    
end module transportNuclearData_inter
