module uniFissSitesField_class

  use numPrecision
  use genericProcedures,     only : fatalError, numToChar
  use universalVariables,    only : OUTSIDE_MAT, VOID_MAT, P_NEUTRON_CE
  use dictionary_class,      only : dictionary
  use particle_class,        only : particle, particleState
  use field_inter,           only : field
  use vectorField_inter,     only : vectorField
  use geometry_inter,        only : geometry
  use RNG_class,             only : RNG

  ! Tally Maps
  use tallyMap_inter,        only : tallyMap
  use tallyMapFactory_func,  only : new_tallyMap

  ! Nuclear Data
  use neutronMaterial_inter, only : neutronMaterial, neutronMaterial_CptrCast
  use nuclearDataReg_mod,    only : ndReg_getNeutronCE => getNeutronCE, &
                                    ndReg_getNeutronMG => getNeutronMG
  use nuclearDatabase_inter, only : nuclearDatabase

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: uniFissSitesField_TptrCast

  !!
  !! Uniform Fission Sites Field
  !!
  !! Returns a 3D vector with:
  !! - fraction of fissionable material volume occupied by the required cell
  !! - percentage of fission sites in the required cell
  !! - unused entry, filled with a contant. It is there to match the interface
  !!
  !! Sample Dictionary Input:
  !!   uniformFissionSites { type uniFissSitesField;
  !!                         #uniformMap 0;#            optional
  !!                         #popVolumes 1.0e7;#        optional
  !!                         map { <map definition> } }
  !!
  !! NOTE: if uniformMap is 0, a routine to estimate the volumes of the map bins
  !!       is run with population popVolumes. This might be inefficient if, for
  !!       example, the volume to be tallied is small compared to the model geometry
  !!
  !! Public Members:
  !!   map ->  map that lays over the geometry. It should be a spatial map, an
  !!           energy map wouldn't make much sense!
  !!   N   ->  total number of map bins
  !!   pop ->  particle population used for the volume estimation
  !!   uniformMap     -> flag to indicate whether the map has bins with uniform volumes
  !!   volFraction    -> array with the volume fraction of each bin
  !!   sourceFraction -> array with the percentage of fission sites in each bin
  !!   buildSource    -> array used to 'tally' fission sites
  !!
  !! Interface:
  !!   vectorField interface
  !!
  type, public, extends(vectorField) :: uniFissSitesField
    private
    class(tallyMap), allocatable :: map
    integer(shortInt)            :: N = 0
    logical(defBool)             :: uniformMap
    integer(shortInt)            :: pop
    real(defReal), dimension(:), allocatable     :: volFraction
    real(defReal), dimension(:), allocatable     :: sourceFraction
    real(defReal), dimension(:), allocatable     :: buildSource
  contains
    ! Superclass interface
    procedure :: init
    procedure :: kill
    procedure :: estimateVol
    procedure :: at
    procedure :: storeFS
    procedure :: updateMap
  end type uniFissSitesField

contains

  !!
  !! Initialise from dictionary
  !!
  !! See field_inter for details
  !!
  subroutine init(self, dict)
    class(uniFissSitesField), intent(inout) :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt), parameter  :: ALL = 0
    character(100), parameter     :: Here = 'init (uniFissSitesField_class.f90)'

    ! Initialise overlay map
    call new_tallyMap(self % map, dict % getDictPtr('map'))
    self % N = self % map % bins(ALL)

    ! Allocate and initialise arrays
    allocate(self % sourceFraction(self % N), self % buildSource(self % N))
    self % sourceFraction = ONE/self % N
    self % buildSource = ZERO

    ! Settings for volume calculation
    call dict % getOrDefault(self % uniformMap,'uniformMap', .true.)
    if (.not. self % uniformMap) call dict % getOrDefault(self % pop,'popVolumes', 1000000)

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(uniFissSitesField), intent(inout) :: self

    call self % map % kill()
    deallocate(self % map)
    deallocate(self % sourceFraction)
    deallocate(self % buildSource)

    self % N = 0

  end subroutine kill

  !!
  !! Generate random points to estimate the volume of the elements on the map
  !!
  subroutine estimateVol(self, geom, rand, type)
    class(uniFissSitesField), intent(inout) :: self
    class(geometry), pointer, intent(in)    :: geom
    class(RNG), intent(inout)               :: rand
    integer(shortInt), intent(in)           :: type
    real(defReal), dimension(6)             :: bounds
    real(defReal), dimension(3)             :: bottom, top
    real(defReal), dimension(3), save       :: rand3, r
    type(particleState), save               :: state
    integer(shortInt)                       :: i
    integer(shortInt), save                 :: j, binIdx, matIdx, uniqueID
    class(nuclearDatabase), pointer         :: nucData
    class(neutronMaterial), pointer, save   :: mat
    character(100), parameter :: Here = 'estimateVol (uniFissSitesField_class.f90)'
    !$omp threadprivate(rand3, r, state, j, binIdx, matIdx, uniqueID, mat)

    allocate(self % volFraction(self % N))

    ! Check if volume estimation is needed or not
    if (self % uniformMap) then
      self % volFraction = ONE/self % N
      return

    else

      ! Get pointer to appropriate nuclear database
      if (type == P_NEUTRON_CE) then
        nucData => ndReg_getNeutronCE()
      else
        nucData => ndReg_getNeutronMG()
      end if
      if(.not.associated(nucData)) call fatalError(Here, 'Failed to retrieve Nuclear Database')

      ! Set bounding region
      bounds = geom % bounds()
      bottom = bounds(1:3)
      top    = bounds(4:6)

      print *, "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>"
      print *, "VOLUME CALCULATION FOR UFS"

      ! Iterate over number of points desired
      !$omp parallel do
      do i = 1, self % pop

        j = 0
        rejection : do
          ! Protect against infinite loop
          j = j +1
          if ( j > 200) then
            call fatalError(Here, 'Infinite loop in sampling of fission sites. Please check that&
                                  & defined volume contains fissile material.')
          end if

          ! Sample Position
          rand3(1) = rand % get()
          rand3(2) = rand % get()
          rand3(3) = rand % get()
          r = (top - bottom) * rand3 + bottom

          ! Find material under position
          call geom % whatIsAt(matIdx, uniqueID, r)

          ! Reject if there is no material
          if (matIdx == VOID_MAT .or. matIdx == OUTSIDE_MAT) cycle rejection

          mat => neutronMaterial_CptrCast(nucData % getMaterial(matIdx))
          if (.not.associated(mat)) call fatalError(Here, "Nuclear data did not return neutron material.")

          ! Resample position if material is not fissile
          if (.not. mat % isFissile()) cycle rejection

          state % r = r

          ! Read map bin index
          binIdx = self % map % map(state)

          ! Return if invalid bin index
          if (binIdx == 0) cycle rejection

          ! Add point to the volume fraction map
          !$omp atomic
          self % volFraction(binIdx) = self % volFraction(binIdx) + 1

          ! Exit the loop
          exit rejection

        end do rejection

      end do
      !$omp end parallel do

      ! Normalise the volume fraction map
      self % volFraction = self % volFraction/sum(self % volFraction)

      print *, "DONE!"
      print *, "<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>"

    end if

  end subroutine estimateVol

  !!
  !! Get value of the vector field given the phase-space location of a particle
  !!
  !! See vectorField_inter for details
  !!
  function at(self, p) result(val)
    class(uniFissSitesField), intent(in) :: self
    class(particle), intent(inout)       :: p
    real(defReal), dimension(3)          :: val
    type(particleState)                  :: state
    integer(shortInt)                    :: binIdx

    ! Get current particle state
    state = p

    ! Read map bin index
    binIdx = self % map % map(state)

    ! Return if invalid bin index
    if (binIdx == 0) then
      val = ONE
      return
    end if

    val(1) = self % volFraction(binIdx)
    val(2) = self % sourceFraction(binIdx)
    val(3) = ZERO

  end function at

  !!
  !! Store the fission sites generated in a vector
  !!
  !! Args:
  !! state [in] -> particle state of the fission site
  !!
  subroutine storeFS(self, state)
    class(uniFissSitesField), intent(inout) :: self
    type(particleState), intent(in)         :: state
    integer(shortInt)                       :: idx

    idx = self % map % map(state)
    if (idx == 0) return
    ! Add fission sites where appropriate
    self % buildSource(idx) = self % buildSource(idx) + state % wgt

  end subroutine storeFS

  !!
  !! Calculates fission site probability distribution
  !! It accounts for possible ares of the map having zero events
  !!
  subroutine updateMap(self)
    class(uniFissSitesField), intent(inout) :: self
    integer(shortInt) :: i

    ! Eliminate zeros in the distribution
    do i = 1, self % N
      if (self % buildSource(i) == ZERO) self % buildSource(i) = ONE
    end do
    ! Normalise to calculate probability
    self % sourceFraction = self % buildSource / sum(self % buildSource)

    self % buildSource = ZERO

  end subroutine updateMap

  !!
  !! Cast field pointer to uniFissSitesField pointer
  !!
  !! Args:
  !!   source [in] -> source pointer of class field
  !!
  !! Result:
  !!   Null is source is not of uniFissSitesField
  !!   Pointer to source if source is uniFissSitesField type
  !!
  pure function uniFissSitesField_TptrCast(source) result(ptr)
    class(field), pointer, intent(in) :: source
    type(uniFissSitesField), pointer  :: ptr

    select type (source)
      type is (uniFissSitesField)
        ptr => source

      class default
        ptr => null()
    end select

  end function uniFissSitesField_TptrCast


end module uniFissSitesField_class
