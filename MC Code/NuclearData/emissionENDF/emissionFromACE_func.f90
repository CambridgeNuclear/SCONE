module emissionFromACE_func

  use numPrecision
  use genericProcedures, only : linFind, searchError

  ! Import emission objects
  use emissionENDF_class,             only : emissionENDF
  use uncorrelatedEmissionENDF_class, only : uncorrelatedEmissionENDF
  use correlatedEmissionENDF_class,   only : correlatedEmissionENDF

  ! Import Angular Distributions
  use angleLawENDF_class,             only : angleLawENDF
  use noAngle_class,                  only : noAngle
  use tabularAngle_class,             only : tabularAngle

  ! Import Mu Probability Density Functions (PDF)
  use muEndfPdf_class,                only : muEndfPdf, muEndfPdf_ptr
  use equiBin32Mu_class,              only : equiBin32Mu
  use isotropicMu_class,              only : isotropicMu
  use tabularMiu_class,               only : tabularMiu

  ! Import Energy distributions
  use energyLawENDF_class,            only : energyLawENDF
  use noEnergy_class,                 only : noEnergy
  use contTabularEnergy_class,        only : contTabularEnergy

  ! Import Energy Probability Density Functions (PDF)
  use tabularEnergy_class,            only : tabularEnergy

  ! Import Correlated Angle-Energy Distributions


  ! Import Neutron Release Distributions
  use releaseLawENDF_class,           only : releaseLawENDF
  use constantRelease_class,          only : constantRelease
  use polynomialRelease_class,        only : polynomialRelease
  use tabularRelease_class,           only : tabularRelease

  implicit none
  private

  public   :: emissionFromACE  ! Make the function avalible to import

contains

  function emissionFromACE(NXS,JXS,XSS,MT) result(new)
    !! Function that for a given ACE libraries and MT number build and returns pointer to a
    !! emission ENDF object
    integer(shortInt), dimension(16),intent(in)  :: NXS
    integer(shortInt), dimension(32),intent(in)  :: JXS
    real(defReal), dimension(:),intent(in)       :: XSS
    integer(shortInt), intent(in)                :: MT
    class(emissionENDF),pointer                  :: new
    character(100),parameter                     :: Here='emissionFromACE (emissionFromACE_func.f90)'
    integer(shortInt)                            :: nMT   ! Number of MT reactions in library
    integer(shortInt)                            :: mtIdx ! Index of MT reaction in ACE library
    integer(shortInt)                            :: TY    ! TY value for the MT
    integer(shortInt)                            :: LOCB  ! Location of Angular distribution
    integer(shortInt)                            :: LOCC  ! Locator  of Energy distribution

    !nMT = NXS(4)
    call setToElastic()
    mtIdx = findMtIdx()

    TY = XSS(JXS(5)+mtIdx-1)

    if (MT == 2) then
      LOCB = 1         ! Elastic scattering. Define ENDF reaction numbers as parameters in a module!!!!
    else
      LOCB = XSS(JXS(8)+mtIdx)
    end if

    LOCC = XSS(JXS(10)+mtIdx-1)

    print *, MT,mtIdx,TY,LOCB,LOCC

    contains

      subroutine setToCapture()
        class(angleLawENDF),pointer    :: angleLaw
        class(energyLawENDF), pointer  :: energyLaw
        class(releaseLawENDF), pointer :: releaseLaw

        angleLaw   => noAngle()
        energyLaw  => noEnergy()
        releaseLaw => constantRelease(0.0_defReal)

        new => uncorrelatedEmissionENDF(angleLaw,energyLaw,releaseLaw)

      end subroutine setToCapture

      subroutine setToElastic()
        class(angleLawENDF),pointer    :: angleLaw
        class(energyLawENDF), pointer  :: energyLaw
        class(releaseLawENDF), pointer :: releaseLaw

        energyLaw  => noEnergy()
        releaseLaw => constantRelease(1.0_defReal)
        angleLaw => readAngleArray(JXS(9))

      end subroutine setToElastic


      function readAngleArray(addr) result(angleLaw)
        !! Function returns pointer to angular law that begins at a given address (addr) in XSS
        !! table.
        integer(shortInt),intent(in)                 :: addr
        integer(shortInt)                            :: numE     ! Number of energy points
        real(defReal),dimension(:),allocatable       :: eGrid    ! Energy Grid
        integer(shortInt),dimension(:),allocatable   :: muAddr   ! Addresses of Mu PDFs
        type(muEndfPdf_ptr),dimension(:),allocatable :: muPdfs   ! Pointers to Mu PDFs
        class(angleLawENDF),pointer                  :: angleLaw ! Pointer to Angle Law to return
        integer(shortInt)                            :: i
        real(defReal)                                :: temp

        ! Check if the number of energy points, which should be under addr is an integer

        temp = XSS(addr)


        numE = int( XSS(addr) )
        ! Perform


        allocate( eGrid(numE) )
        allocate( muAddr(numE))
        allocate( muPdfs(numE))

        eGrid =       XSS(addr+1      : addr+1+numE   )
        muAddr = int( XSS(addr+1+numE : addr+1+2*numE ) )

        do i=1,numE
        !  muPdfs(i) = readMuPdf(muAddr(i))
        end do

       ! angleLaw => tabularAngle(eGrid,muPdfs)

      end function readAngleArray






      function findMtIdx()
        !! Function that finds index of the given MT number in the isotope ACE library
        integer(shortInt) :: findMTIdx
        real(defReal)     :: mtReal

        mtReal = real(MT,defReal)

        findMtIdx = linFind(XSS(JXS(3) : JXS(3)+nMT-1),mtReal)
        call searchError(findMtIdx,Here)

      end function findMtIdx


  end function emissionFromACE




end module emissionFromACE_func
