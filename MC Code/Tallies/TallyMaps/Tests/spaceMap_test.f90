module spaceMap_test
  use numPrecision
  use universalVariables
  use pFUnit_mod
  use particle_class,       only : particleState
  use dictionary_class,     only : dictionary
  use spaceMap_class,       only : spaceMap

  implicit none


@testParameter
  type, extends(AbstractTestParameter) :: dirPar
    integer(shortInt) :: dir
  contains
    procedure :: toString
  end type dirPar


@testCase(constructor=newTest)
  type, extends(ParameterizedTestCase) :: test_spaceMap
    private
    integer(shortInt) :: dir
    type(spaceMap) :: map_struct
    type(spaceMap) :: map_unstruct

  end type test_spaceMap

contains

  !!
  !! Returns an array of test parameters
  !!
  function getParameters() result(params)
    type(dirPar), dimension(3) :: params

    params(1) % dir = X_AXIS
    params(2) % dir = Y_AXIS
    params(3) % dir = Z_AXIS

  end function getParameters

  !!
  !! Returns only direction x
  !!
  function XdirParameter() result(params)
    type(dirPar), dimension(1) :: params

    params(1) % dir = X_AXIS

  end function XdirParameter

  !!
  !! Write test parameter to string
  !!
  function toString(this) result(string)
    class(dirPar), intent(in) :: this
    character(:), allocatable :: string
    character(1)              :: str

    select case(this % dir)
      case(X_AXIS)
        str ='x'
      case(Y_AXIS)
        str='y'
      case(Z_AXIS)
        str='z'
      case default
        str='?'
    end select
    string = str

  end function toString


  !!
  !! Construct test case
  !!
  function newTest(testParam) result(tst)
    type(dirPar), intent(in) :: testParam
    type(test_spaceMap)      :: tst
    type(dictionary)         :: tempDict
    real(defReal),dimension(*),parameter :: BIN_DIV = [-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, &
                                                        2.0,   4.0,  6.0,  8.0,  10.0]

    ! Load direction
    tst % dir = testParam % dir

    ! Create structured grid
    call tempDict % init(5)
    call tempDict % store('axis',testParam % toString())
    call tempDict % store('grid','lin')
    call tempDict % store('min', -10.0_defReal)
    call tempDict % store('max', 10.0_defReal)
    call tempDict % store('N', 20)

    call tst % map_struct % init(tempDict)
    call tempDict % kill()

    ! Create unstructured grid
    call tempDict % init(3)
    call tempDict % store('axis',testParam % toString())
    call tempDict % store('grid','unstruct')
    call tempDict % store('bins', BIN_DIV)

    call tst % map_unstruct % init(tempDict)
    call tempDict % kill()

  end function newTest

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test structured grid
  !!
@Test(testParameters={getParameters()})
  subroutine testStructuredGrid(this)
    class(test_spaceMap), intent(inout)      :: this
    real(defReal),dimension(2),parameter     :: POS = [0.5_defReal, -10.1_defReal]
    integer(shortInt),dimension(2),parameter :: RES = [11, 0]
    integer(shortInt),dimension(2)           :: idx
    type(particleState),dimension(2)         :: states

    states % r(this % dir) = POS
    idx = this % map_struct % map(states)

    @assertEqual(RES, idx)

  end subroutine testStructuredGrid

  !!
  !! Test unstructured grid
  !!
@Test(testParameters={getParameters()})
  subroutine testUnstructuredGrid(this)
    class(test_spaceMap), intent(inout)      :: this
    real(defReal),dimension(2),parameter     :: POS = [0.5_defReal, -10.1_defReal]
    integer(shortInt),dimension(2),parameter :: RES = [6, 0]
    integer(shortInt),dimension(2)           :: idx
    type(particleState),dimension(2)         :: states

    states % r(this % dir) = POS
    idx = this % map_unstruct % map(states)

    @assertEqual(RES, idx)

  end subroutine testUnstructuredGrid

 !!
 !! Test bin output
 !!
@Test(testParameters ={XdirParameter()})
  subroutine testBins(this)
    class(test_spaceMap), intent(inout) :: this

    ! Structured grid
    @assertEqual(20, this % map_struct % bins(1),'Normal use')
    @assertEqual(20, this % map_struct % bins(0),'All binbs')
    @assertEqual(0, this % map_struct % bins(-2),'Invalid dimension')

    ! Unstructured grid
    @assertEqual(10, this % map_unstruct % bins(1),'Normal use')
    @assertEqual(10, this % map_unstruct % bins(0),'All bins')
    @assertEqual(0, this % map_unstruct % bins(-2),'Invalid dimension')

  end subroutine testBins

end module spaceMap_test
