module emissionFromACE_func

  use numPrecision
  use genericProcedures, only : linFind, searchError, isInteger, fatalError
  use endfConstants

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
  use tabularMu_class,                only : tabularMu

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
        angleLaw   => readAngleArray(JXS(9))

        new => uncorrelatedEmissionENDF(angleLaw,energyLaw,releaseLaw)

      end subroutine setToElastic


      subroutine setToMT()
        class(angleLawENDF), pointer     :: angleLaw
        class(energyLawENDF), pointer    :: energyLaw
        class(releaseLawENDF), pointer   :: releaseLaw
        !class(correlatedLawENDF),pointer :: correlatedLaw

      end subroutine setToMT

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

        ! Read number of energy points and verify that it is +ve integer
        temp = XSS(addr)

        if (.not. (temp > 0 .and. isInteger(temp))) then
          call fatalError(Here,'Beginning of angular data record is not +ve int. Adress must be wrong.')
        end if

        numE = temp

        allocate( eGrid(numE) )
        allocate( muAddr(numE))
        allocate( muPdfs(numE))

        ! Read energy points and muPdfs adresses
        eGrid =       XSS(addr+1      : addr+1+numE-1 )
        muAddr = int( XSS(addr+1+numE : addr+1+2*numE ) )

        do i=1,numE
          muPdfs(i) = readMuPdf(muAddr(i))
        end do

        angleLaw => tabularAngle(eGrid,muPdfs)

      end function readAngleArray

      function readMuPdf(addr) result(muPdf)
        integer(shortInt), intent(in)          :: addr
        class(muEndfPdf), pointer              :: muPdf
        real(defReal),dimension(:),allocatable :: muGrid, pdfGrid, cdfGrid
        real(defReal)                          :: temp
        integer(shortInt)                      :: interFlag, numE, locAddr

        integer(shortInt),parameter            :: histogram  = tabPdfHistogram, &
                                                  linLin     = tabPdfLinLin
        character(100),parameter               :: Here='readMuPdf (emissionFromACE_func.f90)'


        if (addr < 0) then
        ! Tabular Energy Distribution
          locAddr = abs(addr)

          ! Read interpolation flag and valide it
          interFlag = XSS(JXS(9)+locAddr-1)

          if (interFlag /= histogram .and. interFlag /= linLin) then
            call fatalError(Here,'Invalid interpolation parameter was read.')
          end if

          ! Read number of energy points and validate it
          temp = XSS(JXS(9)+locAddr)

          if (.not.(isInteger(temp))) call fatalError(Here,'Number of energy points is not'// &
                                                           'integer')

          numE = temp

          if (numE <= 0)              call fatalError(Here,'Number of energy points is not'// &
                                                           'strictly +ve')

          ! Read muGrid and pdfGrid
          muGrid  = XSS( JXS(9)+locAddr+1        : JXS(9)+locAddr+1+numE-1   )
          pdfGrid = XSS( JXS(9)+locAddr+1+numE   : JXS(9)+locAddr+1+2*numE-1 )
          cdfGrid = XSS( JXS(9)+locAddr+1+2*numE : JXS(9)+locAddr+1+3*numE-1 )

          muPdf => tabularMu(muGrid,pdfGrid,cdfGrid,interFlag)

        else if(addr > 0) then
        ! 32 equiprobable bins Energy Distribution
          locAddr = abs(addr)

          allocate(muGrid(33))

          muGrid = XSS(JXS(9)+locAddr-1 : JXS(9)+locAddr+33-1)

          muPdf => equiBin32Mu(muGrid)

        else if(addr == 0) then
        ! Isotropic energy distribution
          muPdf => isotropicMu()
        else
        ! Imppossible state
          call fatalError(Here,'addr in readMuPdf is not -ve, not +ve and not 0. WTF?')
        end if

      end function




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
