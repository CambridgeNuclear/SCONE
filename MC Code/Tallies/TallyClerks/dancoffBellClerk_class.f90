module dancoffBellClerk_class

  use numPrecision
  use tallyCodes
  use genericProcedures,          only : fatalError, hasDuplicates
  use dictionary_class,           only : dictionary
  use outputFile_class,           only : outputFile
  use intMap_class,               only : intMap

  use particle_class,             only : particle, phaseCoord, particleState
  use particleDungeon_class,      only : particleDungeon

  ! Nuclear Data
  use nuclearDataRegistry_mod,    only : getMatIdx, getMatName
  use transportNuclearData_inter, only : transportNuclearData

  ! Tally Interfaces
  use tallyEstimator_class,       only : tallyScore, tallyCounter
  use tallyClerk_inter,           only : tallyClerk
  use materialMap_class,          only : materialMap
  use energyMap_class,            only : energyMap

  implicit none
  private

  character(*),parameter :: CLASS_NAME = 'dancoffBellClerk'

  !! Local parameters. Note that OUSIDE is not arbitraty. Choosen to match invalid idx
  !! given by a tallyMap
  integer(shortInt), parameter :: FUEL      = -2, &
                                  MODERATOR = -3, &
                                  OUTSIDE   = 0

  !!
  !! Constructor
  !!
  interface dancoffBellClerk
    module procedure dancoffBellClerk_fromDict
  end interface


  !!
  !! Special tally clerk design to tally combined dancoff and bell factor
  !! multiplied by escepe probability SIGMA_e = 1/ L_bar, where L_bar is
  !! an average cord of a fuel lump.
  !!
  !! TODO: Add additional explenation and refer to Nuclear Engineering Handbook
  !!       for additional details. Also try to find more accesibale reference.
  !!
  !! D_eff = P_e * SIGMA_t / (1-P_e)
  !!
  !! D_eff = wgt_e * SIGMA_t / wgt_mod
  !!
  type, public,extends(tallyClerk) :: dancoffBellClerk
    private
    ! Material for total XS
    ! Materials list for fuel
    ! Materials list for moderator
    character(nameLen)   :: name
    integer(shortInt)    :: cycleCount  = 0                ! Cycles counter
    real(defReal)        :: targetRelSD = 0.0

    ! Maping variables
    logical(defBool)     :: hasMap =.false.
    type(energyMap)      :: map
    type(intMap)         :: materialSet
    integer(shortInt)    :: xsMatIdx

    ! Result bins
    integer(shortInt)                           :: Nbins
    type(tallyScore),dimension(:),allocatable   :: escSigmaT
    type(tallyScore),dimension(:),allocatable   :: modWgt
    type(tallyCounter),dimension(:),allocatable :: D_eff


    contains
    ! Deferred Interface Procedures
    procedure :: validReports
    procedure :: display
    procedure :: isConverged
    procedure :: init
    procedure :: print

    ! Overwrite report procedures
    procedure :: reportTrans
    procedure :: reportCycleEnd

  end type dancoffBellClerk

contains

  !!
  !! Return codes for reports this clerk accepts
  !!
  function validReports(self) result(validCodes)
    class(dancoffBellClerk), intent(in)              :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [trans_CODE, cycleEnd_CODE]

  end function validReports

  !!
  !! Display progress. Print total estimate and maximum relative STD s
  !!
  subroutine display(self)
    class(dancoffBellClerk), intent(in)         :: self
    real(defReal)                            :: resSum, maxSTD

    print '(A)', 'dancoffBellClerk does not support display'

  end subroutine display

  !!
  !! Perform convergance check in the Clerk
  !!
  function isConverged(self) result(isIt)
    class(dancoffBellClerk), intent(in)         :: self
    logical(defBool)                         :: isIt
    real(defReal)                            :: maxSTD

    ! This functionality is not yet implemented
    isIt = .false.

  end function isConverged

  !!
  !! Initialise dancoffBellClerk from dictionary
  !!
  subroutine init(self,dict,name)
    class(dancoffBellClerk),intent(inout)          :: self
    class(dictionary), intent(in)                  :: dict
    character(nameLen), intent(in)                 :: name
    character(nameLen)                             :: tempChar
    character(nameLen),dimension(:),allocatable    :: tempCharArr
    integer(shortInt),dimension(:),allocatable     :: fuelMatIdx
    integer(shortInt),dimension(:),allocatable     :: modMatIdx
    integer(shortInt)                              :: N, i
    character(100),parameter :: Here ='init (dancoffBellClerk_class.f90)'

    ! Assign name
    self % name = name

    ! Store index of XS material
    call dict % get(tempChar,'XSmat')
    self % xsMatIdx = getMatIdx(tempChar)

    ! Get fuel material elements
    call dict % get(tempCharArr, 'fuelMat')
    allocate(fuelMatIdx (size(tempCharArr)))
    do i = 1,size(fuelMatIdx)
      fuelMatIdx(i) = getMatIdx(tempCharArr(i))

    end do

    ! Get moderator material elements
    call dict % get(tempCharArr, 'modMat')
    allocate(modMatIdx (size(tempCharArr)))
    do i = 1,size(modMatIdx)
      modMatIdx(i) = getMatIdx(tempCharArr(i))

    end do

    ! Check for overlap
    if (hasDuplicates([modMatIdx, fuelMatIdx])) then
      call fatalError(Here, 'Same materials are defined as fuel and moderator')

    end if

    ! Create map of matIdx -> matType
    call self % materialSet % init (size(modMatIdx) + size(fuelMatIdx))

    ! Load Fuel
    do i = 1,size(fuelMatIdx)
      call self % materialSet % add(fuelMatIdx(i), FUEL)

    end do

    ! Load Moderator
    do i = 1,size(modMatIdx)
      call self % materialSet % add(modMatIdx(i), MODERATOR)

    end do

    ! Check if has map and load the map
    if( dict % isPresent('energyMap')) then
      self % hasMap = .true.
      self % map = energyMap( dict % getDictPtr('energyMap'))

    else
      self % hasMap = .false.

    end if
    ! Create Space for the results
    ! Find size of map
    if(self % hasMap) then
      N = self % map % bins()
    else
      N = 1
    end if

    ! Create space for results
    allocate(self % escSigmaT(N))
    allocate(self % modWgt(N) )
    allocate(self % D_eff(N)    )

  end subroutine init

  !!
  !! Write contents of the dancoffBellClerk to output file
  !!
  subroutine print(self,outFile)
    class(dancoffBellClerk), intent(in)   :: self
    class(outputFile), intent(inout)      :: outFile
    real(defReal)                         :: val, std
    integer(shortInt)                     :: i, N
    character(nameLen)                    :: name

    ! Begin block
    call outFile % startBlock(self % name)

    if (self % hasMap) then
      ! Write axis data
      call self % map % print(outFile)

      ! Write axis descriptor
      name = 'axisDescriptor'
      call outFile % startArray(name,[1])
      call outFile % addValue(self % map % getAxisName())
      call outFile % endArray()
    end if

    ! Write results
    name = 'Res'
    N = self % map % bins()
    call outFile % startArray(name,[N])

    do i=1,size(self % D_eff)
      call self % D_eff(i) % getEstimate(val, std, self % cycleCount)
      call outFile % addResult(val,std)
    end do

    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

  !!
  !! Process transition report
  !!
  subroutine reportTrans(self,p)
    class(dancoffBellClerk), intent(inout)  :: self
    class(particle), intent(in)             :: p
    real(defReal)                           :: SigmaTot
    integer(shortInt)                       :: T_end, T_start, B
    real(defReal)                           :: w_end
    character(100),parameter :: Here = 'reportTrans (dancoffBellClerk_class.f90)'

    ! Find start material type; Exit if not fuel
    T_start = self % materialSet % getOrDefault(p % preTransition % matIdx, OUTSIDE)
    if( T_start /= FUEL) return

    ! Find bin
    if(self % hasMap) then
      B = self % map % map(p)
    else
      B = 1
    end if
    if(B == 0) return

    ! Find end material type; Exit if not fuel or moderator
    T_end = self % materialSet % getOrDefault(p % matIdx(), OUTSIDE)
    if(T_end == OUTSIDE) return

    ! Obtain starting and ending weights
    w_end   = p % w

    ! Add to approperiate bins
    select case(T_end)
      case(FUEL)
        ! Obtain XSs
        ! Check if it dynamic type is supported
        ! If it is obtain macroscopic XSs
        ! It it isn't throw error
        associate (xsData => p % xsData)
          select type(xsData)
            class is (transportNuclearData)
              SigmaTot = xsData % getTotalMatXS(p, self % xsMatIdx)

            class default
              call fatalError(Here,'Dynamic type of XS data attached to particle is not transportNuclearData')
          end select
        end associate

        call self % escSigmaT(B) % add(w_end * SigmaTot)

      case(MODERATOR)
        call self % modWgt(B) % add(w_end)

      case default
        call fatalError(Here, 'WTF? Impossible state')

    end select

  end subroutine reportTrans


  !!
  !! Process end of cycle report
  !!
  subroutine reportCycleEnd(self,end)
    class(dancoffBellClerk), intent(inout)   :: self
    class(particleDungeon), intent(in)       :: end
    integer(shortInt)                        :: i,N

    self % cycleCount = self % cycleCount + 1

    ! Store results. Behold the glorious Fortran oneliner! Power of the Elemental Procedures!
    ! NOTE: It may tend to create an unnecessary temp array...
    call self % D_eff % addEstimate( self % escSigmaT % get() / self % modWgt % get() )

    call self % escSigmaT % reset()
    call self % modWgt % reset()

  end subroutine reportCycleEnd

  !!
  !! Return an instance of mactoClerk from dictionary and name
  !!
  function dancoffBellClerk_fromDict(dict,name) result(new)
    class(dictionary), intent(in)  :: dict
    character(nameLen), intent(in) :: name
    type(dancoffBellClerk)         :: new

    call new % init(dict,name)

  end function dancoffBellClerk_fromDict
    
end module dancoffBellClerk_class
