module collisionClerk_test

  use numPrecision
  use genericProcedures,              only : numToChar
  use collisionClerk_class,           only : collisionClerk
  use particle_class,                 only : particle
  use dictionary_class,               only : dictionary
  use scoreMemory_class,              only : scoreMemory
  use testNeutronDatabase_class,      only : testNeutronDatabase
  use outputFile_class,               only : outputFile
  use ceNeutronCache_mod,             only: cache_init => init, trackingCache
  use funit

  implicit none

@testParameter(constructor = new_testNumber)
  type, extends(AbstractTestParameter) :: testNumber
    integer(shortInt) :: i
  contains
    procedure :: toString
  end type testNumber

@testCase(constructor=newTest)
  type, extends(ParameterizedTestCase) :: test_collisionClerk
    private
    logical(defBool)                            :: hasFilter
    logical(defBool)                            :: hasMap
    logical(defBool)                            :: has2Res
    integer(longInt), dimension(:), allocatable :: bins
    real(defReal), dimension(:), allocatable    :: results
  end type test_collisionClerk

contains

  !!
  !! Build new test parameter form integer
  !!
  function new_testNumber(i) result (tstNum)
    integer(shortInt) :: i
    type(testNumber)  :: tstNum

    tstNum % i = i

  end function new_testNumber

  !!
  !! Write test parameter to string
  !!
  function toString(this) result(string)
    class(testNumber), intent(in) :: this
    character(:), allocatable :: string
    character(nameLen)        :: str

    write (str,*) this % i
    string = str

  end function toString

  !!
  !! Construct test case
  !!
  function newTest(testParam) result(tst)
    type(testNumber), intent(in) :: testParam
    type(test_collisionClerk)    :: tst
    integer(shortInt)            :: testNum, Nbins, i
    real(defReal)                :: score1, score2

    testNum = testParam % i - 1

    ! Set test parameters
    tst % hasFilter = mod(testNum, 2)   == 1
    tst % hasMap    = mod(testNum/2, 2) == 1
    tst % has2Res   = mod(testNum/4, 2) == 1

    ! Load expected results

    ! Determine number of bins
    Nbins = 1
    if(tst % has2Res) Nbins = Nbins * 2
    if(tst % hasMap)  Nbins = Nbins * 7

    ! Allocate result arrays
    tst % bins    = [(int(i,longInt), i=1,Nbins)]
    tst % results = [(ZERO, i=1,Nbins)]

    ! Set appropriate results (wgt * 1/totXs)
    score1 = 0.7_defReal / 0.3_defReal
    score2 = 1.3_defReal / 0.3_defReal

    select case(testParam % i)
      case(1) ! Single Bin, both score to the same
        tst % results(1) = score1 + score2

      case(2) ! Single bin with filter
        tst % results(1) = score1

      case(3) ! Map without filter
        tst % results(1) = score1
        tst % results(6) = score2

      case(4) ! Map with filter
        tst % results(1) = score1

      case(5) ! Single bin with 2nd Response
        tst % results(1) = score1 + score2
        tst % results(2) = score1 * 1.3_defReal + score2 * 1.3_defReal

      case(6) ! Single bin with filter and 2nd Response
        tst % results(1) = score1
        tst % results(2) = score1 * 1.3_defReal

      case(7) ! Map with 2nd response
        tst % results(1)  = score1
        tst % results(2)  = score1 * 1.3_defReal
        tst % results(11) = score2
        tst % results(12) = score2 * 1.3_defReal

      case(8) ! Map with 2nd response and filter
        tst % results(1)  = score1
        tst % results(2)  = score1 * 1.3_defReal
    end select
  end function newTest

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Scoring test
  !!
@Test(cases=[1,2,3,4,5,6,7,8])
  subroutine testScoring(this)
    class(test_collisionClerk), intent(inout) :: this
    logical(defBool)                          :: hasFilter, hasMap, has2Res
    character(:),allocatable                  :: case
    type(collisionClerk)                      :: clerk
    type(scoreMemory)                         :: mem
    type(particle)                            :: p
    type(testNeutronDatabase)                 :: nucData
    type(outputFile)                          :: outF
    type(dictionary)                          :: filterDict, mapDict, res1Dict, res2Dict, clerkDict
    character(nameLen)                        :: res1Name, res2Name, clerkName
    integer(shortInt)                         :: i
    real(defReal)                             :: res
    real(defReal), parameter :: TOL = 1.0E-9

    ! Copy test settings
    hasFilter = this % hasFilter
    hasMap    = this % hasMap
    has2Res   = this % has2Res

    ! Build case description
    case = 'Vanilla case with: '
    if(hasFilter) case = case // ' Filter '
    if(hasMap)    case = case // ' Map '
    if(has2Res)   case = case // ' 2nd Response '

    ! Define filter dictionary
    call filterDict % init(3)
    call filterDict % store('type','testFilter')
    call filterDict % store('minIdx',0)
    call filterDict % store('maxIdx',5)

    ! Define Map dictionary
    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',7)

    ! Define 1st response dictionary and name
    res1Name = 'flux'
    call res1Dict % init(1)
    call res1Dict % store('type','fluxResponse')

    ! Define 2nd response dictionary and name
    res2Name ='testResponse'
    call res2Dict % init(2)
    call res2Dict % store('type','testResponse')
    call res2Dict % store('value', 1.3_defReal)

    ! Configure dictionary for the clerk
    call clerkDict % init(7)
    call clerkDict % store('type','collisionClerk')
    call clerkDict % store('handleVirtual', 0)
    call clerkDict % store(res1Name, res1Dict)
    call clerkDict % store(res2Name, res2Dict)

    ! Store filter or map
    if(hasFilter) call clerkDict % store('filter', filterDict)
    if(hasMap)    call clerkDict % store('map', mapDict)

    ! Store responses used
    if(has2Res) then
      call clerkDict % store('response', [res1Name, res2Name])
    else
      call clerkDict % store('response', [res1Name])
    end if

    ! Build Clerk
    clerkName ='myClerk'
    call clerk % init(clerkDict, clerkName)

    ! Create score memory
    call mem % init(int(clerk % getSize(), longInt) , 1)
    call clerk % setMemAddress(1_longInt)

    ! Build nuclear data
    call nucData % build(0.3_defReal)

    ! Perform scoring
    call p % setMatIdx(1)
    p % w = 0.7_defReal
    call clerk % reportInColl(p, nucData, mem, .false.)

    call p % setMatIdx(6)
    p % w = 1.3_defReal
    call clerk % reportInColl(p, nucData, mem, .false.)

    ! Virtual scoring should not contribute to score in this case
    p % w = 1000.3_defReal
    call clerk % reportInColl(p, nucData, mem, .true.)

    call mem % closeCycle(ONE)

    ! Verify results of scoring
    do i=1,size(this % bins)
      call mem % getResult(res, this % bins(i))
      @assertEqual(this % results(i), res, TOL, case // 'BIN : ' //numToChar(i) )
    end do

    ! Verify that size of memory returned is correct
    @assertEqual(size(this % bins), clerk % getSize(), case // 'Memory size test:')

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call clerk % print (outF, mem)

    @assertTrue(outF % isValid(), case)

    ! Clean up
    call nucData % kill()
    call clerkDict % kill()
    call filterDict % kill()
    call mapDict % kill()
    call res1Dict % kill()
    call res2Dict % kill()

  end subroutine testScoring

@Test(cases=[1,2,3,4,5,6,7,8])
  subroutine testScoringVirtual(this)
    class(test_collisionClerk), intent(inout) :: this
    logical(defBool)                          :: hasFilter, hasMap, has2Res
    character(:),allocatable                  :: case
    type(collisionClerk)                      :: clerk
    type(scoreMemory)                         :: mem
    type(particle)                            :: p
    type(testNeutronDatabase)                 :: nucData
    type(outputFile)                          :: outF
    type(dictionary)                          :: filterDict, mapDict, res1Dict, res2Dict, clerkDict
    character(nameLen)                        :: res1Name, res2Name, clerkName
    integer(shortInt)                         :: i
    real(defReal)                             :: res
    real(defReal), parameter :: TOL = 1.0E-9

    ! Copy test settings
    hasFilter = this % hasFilter
    hasMap    = this % hasMap
    has2Res   = this % has2Res

    ! Build case description
    case = 'Vanilla case with: '
    if(hasFilter) case = case // ' Filter '
    if(hasMap)    case = case // ' Map '
    if(has2Res)   case = case // ' 2nd Response '

    ! Define filter dictionary
    call filterDict % init(3)
    call filterDict % store('type','testFilter')
    call filterDict % store('minIdx',0)
    call filterDict % store('maxIdx',5)

    ! Define Map dictionary
    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',7)

    ! Define 1st response dictionary and name
    res1Name = 'flux'
    call res1Dict % init(1)
    call res1Dict % store('type','fluxResponse')

    ! Define 2nd response dictionary and name
    res2Name ='testResponse'
    call res2Dict % init(2)
    call res2Dict % store('type','testResponse')
    call res2Dict % store('value', 1.3_defReal)

    ! Configure dictionary for the clerk
    call clerkDict % init(6)
    call clerkDict % store('type','collisionClerk')
    call clerkDict % store(res1Name, res1Dict)
    call clerkDict % store(res2Name, res2Dict)

    ! Store filter or map
    if(hasFilter) call clerkDict % store('filter', filterDict)
    if(hasMap)    call clerkDict % store('map', mapDict)

    ! Store responses used
    if(has2Res) then
      call clerkDict % store('response', [res1Name, res2Name])
    else
      call clerkDict % store('response', [res1Name])
    end if

    ! Build Clerk
    clerkName ='myClerk'
    call clerk % init(clerkDict, clerkName)

    ! Create score memory
    call mem % init(int(clerk % getSize(), longInt) , 1)
    call clerk % setMemAddress(1_longInt)

    ! Build nuclear data
    call nucData % build(0.3_defReal)
    ! Build cache
    call cache_init(1, 1, 1)
    trackingCache % xs = 0.3_defReal
    trackingCache % E  = 10.0_defReal

    ! Perform scoring, both virtual and physical should contribute
    call p % setMatIdx(1)
    p % w = 0.7_defReal
    p % E = 10.0_defReal
    p % isMG = .false.
    call clerk % reportInColl(p, nucData, mem, .true.)

    call p % setMatIdx(6)
    p % w = 1.3_defReal
    call clerk % reportInColl(p, nucData, mem, .false.)

    call mem % closeCycle(ONE)

    ! Verify results of scoring
    do i=1,size(this % bins)
      call mem % getResult(res, this % bins(i))
      @assertEqual(this % results(i), res, TOL, case // 'BIN : ' //numToChar(i) )
    end do

    ! Verify that size of memory returned is correct
    @assertEqual(size(this % bins), clerk % getSize(), case // 'Memory size test:')

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call clerk % print (outF, mem)

    @assertTrue(outF % isValid(), case)

    ! Clean up
    call nucData % kill()
    call clerkDict % kill()
    call filterDict % kill()
    call mapDict % kill()
    call res1Dict % kill()
    call res2Dict % kill()

  end subroutine testScoringVirtual

end module collisionClerk_test
