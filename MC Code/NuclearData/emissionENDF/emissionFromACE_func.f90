module emissionFromACE_func
  !! This module contains a single function that takes as arguments ACE tables (NXS,JXS,XSS) and
  !! number MT to return pointer to emission distribution object ( emissionENDF ).
  !!
  !! In order to improve readibility function emissionFormACE() uses a Fortran feature called
  !! internal procedures. Basically the function has a "contains" keyword, which is followed  by
  !! definition of a number of functions and subroutines before the "end function emissionFromACE".
  !! Special thing about the internal procedures is that they "see" all local variables defined in
  !! their parent procedure. In this case it allws them to refer to NXS, JXS, XSS and MT.
  !!
  !! For the details about the format of ACE library and NXS,JXS and XSS tables please refer to
  !! MCNP manual Appendix F, which clearly defines the ACE structure.

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
  use isotropicAngle_class,           only : isotropicAngle

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

  integer(shortInt), parameter :: correlated = -1 ,&
                                  isotropic  =  0 ,&
                                  fission    = 19 ,&
                                  absorbtion = 0  ,&
                                  polynomial = 1  ,&
                                  tabular    = 2

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
    class(releaseLawENDF),pointer :: test

    ! Select what data to load
    select case(MT)
      case (N_N_elastic)      ! MT = 2
        call setToElastic()
      case (N_disap : N_da)   ! MT = 1**
        call setToCapture()
      case default
        call setToMT()        ! All other MT numbers
    end select

    ! Return from function. Below after "contains" are internal procedures

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
        !*class(correlatedLawENDF),pointer :: correlatedLaw
        integer(shortInt)                :: nMT   ! Number of MT reactions in library
        integer(shortInt)                :: mtIdx ! Index of MT reaction in ACE library
        integer(shortInt)                :: TY    ! Neutron release (TY value) for the MT
        integer(shortInt)                :: LOCB  ! Locator  of Angular distribution
        integer(shortInt)                :: LOCC  ! Locator  of Energy distribution
        character(100),parameter         :: Here = 'setToMT (emissionFromACE_func.f90)'
        logical(defBool)                 :: isCorrelated
        logical(defBool)                 :: isIsotropic
        logical(defBool)                 :: isInLabFrame
        logical(defBool)                 :: isCapture
        real(defReal)                    :: TY_r,LOCB_r,LOCC_r ! Dummy real temp variables

        ! Read locators and release for the given MT
        nMT    = NXS(4)
        mtIdx  = findMtIdx()
        TY_r   = XSS(JXS(5) +mtIdx-1)
        LOCB_r = XSS(JXS(8) +mtIdx  )
        LOCC_r = XSS(JXS(10)+mtIdx-1)

        ! Convert To integers
        TY   = TY_r
        LOCB = LOCB_r
        LOCC = LOCC_r

        ! Set logical flags depending on the data
        isCorrelated = (LOCB == correlated)
        isIsotropic  = (LOCB == isotropic )
        isInLabFrame = (TY > 0)
        isCapture    = (TY == 0)

        ! Perform Checks
        if (.not.isInteger(TY_r))   call fatalError(Here,'Read neutron release(TY) is not integer')

        ! If the MT reaction is capture then it has no LOCB an LOCC values. Values read are garbage
        ! by default. Checks for sensibility of LOCB and LOCC are necessary
        if (.not. isCapture) then
          if (.not.isInteger(LOCB_r)) call fatalError(Here,'Read angle locator LOCB is not integer')
          if (.not.isInteger(LOCC_r)) call fatalError(Here,'Read energy locator LOCC is not integer')
        end if

        ! Correct sign of release value
        TY = abs(TY)

        ! Read neutron release
        if (TY == fission) then
          !* Call subroutine to read Nu data
          releaseLaw => readNuArray()
        else
          releaseLaw => constantRelease( real(TY,defReal) )
        end if

        ! Create emissionENDF object
        if (isCorrelated) then
          ! Not implemented
          call fatalError(Here,'Correleted emission laws are not yet implemented. Sorry...')

        elseif (isCapture) then
          angleLaw  => noAngle()
          energyLaw => noEnergy()
          new => uncorrelatedEmissionENDF(angleLaw,energyLaw,releaseLaw)

        elseif (isIsotropic) then
          angleLaw  => isotropicAngle()
          energyLaw => readEnergyArray(JXS(11)+LOCC-1)
          new => uncorrelatedEmissionENDF(angleLaw,energyLaw,releaseLaw)

        else
          angleLaw  => readAngleArray( JXS(9)+LOCB-1 )
          energyLaw => readEnergyArray(JXS(11)+LOCC-1)
          new => uncorrelatedEmissionENDF(angleLaw,energyLaw,releaseLaw)

        end if

        ! Set to Labratory coordinate frame
        if (isInLabFrame) call new % setLabFrame()

      end subroutine setToMT



      function readAngleArray(addr) result(angleLaw)
        !! Function returns pointer to angular law that begins at a given address (addr) in XSS
        !! table.
        !! addr should point to JXS(9)+LOCB-1
        integer(shortInt),intent(in)                 :: addr
        integer(shortInt)                            :: numE     ! Number of energy points
        real(defReal),dimension(:),allocatable       :: eGrid    ! Energy Grid
        integer(shortInt),dimension(:),allocatable   :: muAddr   ! Addresses of Mu PDFs
        type(muEndfPdf_ptr),dimension(:),allocatable :: muPdfs   ! Pointers to Mu PDFs
        class(angleLawENDF),pointer                  :: angleLaw ! Pointer to Angle Law to return
        integer(shortInt)                            :: i
        real(defReal)                                :: temp
        character(100),parameter                     :: Here = 'readAngleArray (emissionFromACE_func.f90)'

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



      function readMuPdf(LC) result(muPdf)
        !! Function to read a single mu distribution at an energy point
        !! LC should be locator flag wrt. to JXS(9) read in AND block
        integer(shortInt), intent(in)          :: LC
        class(muEndfPdf), pointer              :: muPdf
        real(defReal),dimension(:),allocatable :: muGrid, pdfGrid, cdfGrid
        real(defReal)                          :: temp
        integer(shortInt)                      :: interFlag, numE, locAddr

        integer(shortInt),parameter            :: histogram  = tabPdfHistogram, &
                                                  linLin     = tabPdfLinLin
        character(100),parameter               :: Here='readMuPdf (emissionFromACE_func.f90)'

        if (LC < 0) then
        ! Tabular Energy Distribution
          locAddr = abs(LC)

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

          if (numE <= 0) call fatalError(Here,'Number of energy points is not strictly +ve')

          ! Read muGrid and pdfGrid
          muGrid  = XSS( JXS(9)+locAddr+1        : JXS(9)+locAddr+1+numE-1   )
          pdfGrid = XSS( JXS(9)+locAddr+1+numE   : JXS(9)+locAddr+1+2*numE-1 )
          cdfGrid = XSS( JXS(9)+locAddr+1+2*numE : JXS(9)+locAddr+1+3*numE-1 )

          muPdf => tabularMu(muGrid,pdfGrid,cdfGrid,interFlag)

        else if(LC > 0) then
        ! 32 equiprobable bins Energy Distribution
          locAddr = abs(LC)

          allocate(muGrid(33))

          muGrid = XSS(JXS(9)+locAddr-1 : JXS(9)+locAddr+33-1)
          print *, muGrid
          muPdf => equiBin32Mu(muGrid)

        else if(LC == 0) then
        ! Isotropic energy distribution
          muPdf => isotropicMu()
        else
        ! Imppossible state
          call fatalError(Here,'addr in readMuPdf is not -ve, not +ve and not 0. WTF?')
        end if

      end function



      function readEnergyArray(addr) result(energyLaw)
        !! Function returns pointer to energy law that begins at a given address (addr) in XSS
        !! table.
        !! addr should point to JXS(11)+LOCC-1
        integer(shortInt),intent(in)   :: addr
        class(energyLawENDF),pointer   :: energyLaw
        integer(shortInt)              :: nLaws                       ! Number of energy Laws
        integer(shortInt)              :: lawType                     ! Type of energy Law
        integer(shortInt)              :: addrLaw                     ! Adress of the Law data
        real(defReal)                  :: nLaws_r,lawType_r,addrLaw_r ! Temp real holders
        character(100),parameter       :: Here='readEnergyArray (emissionFromACE_func.f90)'

        nLaws_r   = XSS(addr)
        lawType_r    = XSS(addr+1)
        addrLaw_r = XSS(addr+2)

        ! Perform Sanity Checks
        if (.not.isInteger(nLaws_r))   call fatalError(Here,'Read number of Laws is not integer')
        if ( nLaws_r < 0)              call fatalError(Here,'-ve number of laws')
        if (.not.isInteger(lawType_r)) call fatalError(Here,'Read Law type is not integer')
        if ( lawType_r < 0)            call fatalError(Here,'Read Law type is -ve integer')
        if (.not.isInteger(addrLaw_r)) call fatalError(Here,'Read addres of energy Law is not int')
        if (addrLaw_r < 0)             call fatalError(Here,'Read adress of energy Law is -ve')

        ! Convert to integers
        nLaws   = nLaws_r
        lawType = lawType_r
        addrLaw = addrLaw_r


        ! Read energy law
        if (nLaws == 0) then
          energyLaw => readSingleEnergyLaw(JXS(11)+addrLaw-1,lawType)
        else
          call fatalError(Here,'Multiple energy laws for a single MT are not yet supported')
        end if

      end function readEnergyArray



      function readSingleEnergyLaw(addr,lawType) result (energyLaw)
        !! Function to select approperiate function for reading a single energy Law
        !! addr should point to LDAT(1) = JXS(11)+IDAT-1
        integer(shortInt), intent(in)  :: addr
        integer(shortInt), intent(in)  :: lawType
        class(energyLawENDF),pointer   :: energyLaw
        character(100),parameter       :: Here='readEnergyLaw (emissionFromACE_func.f90)'

        select case (lawType)
          case (continuousTabularDistribution)
            energyLaw => readEnergyLaw_contTabularEnergy(addr)
          case default
            print *, 'Energy Law Type :', lawType
            call fatalError(Here,'Energy Law Type is not recognised')
        end select

      end function readSingleEnergyLaw


      function readEnergyLaw_contTabularEnergy(addr) result(energyLaw)
        !! Reads energy Law in continuous tabular distribution form.
        !! In the future this subroutine should be moved to object definition
        !! Argument addr should point to LDAT(1) -> NR -> Number of interpolation regions
        integer(shortInt), intent(in)                :: addr
        class(energyLawENDF),pointer                 :: energyLaw
        integer(shortInt)                            :: interRegions  , numE
        real(defReal)                                :: interRegions_r, numE_r
        integer(shortInt)                            :: i
        real(defReal),dimension(:),allocatable       :: eGrid
        integer(shortInt),dimension(:),allocatable   :: ePdfLoc         ! Locators of energy PDFs
        type(tabularEnergy),dimension(:),allocatable :: ePdfs
        character(100),parameter                     :: Here ='readEnergyLaw_contTabularEnergy &
                                                             & (emissionFromACE_func.f90)'
        ! Read and check number of interpolation paramethers
        interRegions_r = XSS(addr)
        if(.not.isInteger(interRegions_r)) then
          call fatalError(Here,'Number of interpolation regions is not an integer')
        end if

        interRegions = interRegions_r

        if(interRegions /= 0) then
          call fatalError(Here,'Many inter. regions on energy distr. table are not supported')
        end if

        ! Read and check number of energy points
        numE_r = XSS(addr+1)
        if(.not.isInteger(numE_r)) call fatalError(Here,'Number of energy points is not an integer')
        if(numE_r < 0 )            call fatalError(Here,'-ve number of energy points')
        numE = numE_r

        ! Read Data
        eGrid   = XSS(addr+2      : addr+2+numE-1   )
        ePdfLoc = XSS(addr+2+numE : addr+2+2*numE-1 )

        allocate(ePdfs(numE))

        do i=1,numE
          ePdfs(i) = readSingleEnergyTabularPdf(JXS(11)+ePdfLoc(i)-1)
        end do

        energyLaw => contTabularEnergy(eGrid,ePdfs)

      end function readEnergyLaw_contTabularEnergy



      function readSingleEnergyTabularPdf(addr) result(tabularEnergyPdf)
        integer(shortInt),intent(in)           :: addr
        type(tabularEnergy)                    :: tabularEnergyPdf
        integer(shortInt)                      :: interType, nPoints
        real(defReal)                          :: interType_r, nPoints_r
        real(defReal),dimension(:),allocatable :: eGrid, PDF, CDF
        character(100),parameter               :: Here ='readSingleEnergyTabularPdf &
                                                        & (emissionFromACE_func.f90)'
        ! Read and check number of points and interpolation type
        interType_r = XSS(addr)
        nPoints_r   = XSS(addr+1)

        if(.not.isInteger(interType_r)) call fatalError(Here,'Interpolation type is not an int')
        if(.not.isInteger(nPoints_r))   call fatalError(Here,'Number of PDF points is not an int')
        if(nPoints_r < 0)               call fatalError(Here,'-ve number of PDF points')

        interType = interType_r
        nPoints   = nPoints_r

        ! Read arrays
        eGrid = XSS(addr+2           : addr+2+nPoints-1  )
        PDF   = XSS(addr+2+nPoints   : addr+2+2*nPoints-1)
        CDF   = XSS(addr+2+2*nPoints : addr+2+3*nPoints-1)

        call tabularEnergyPdf % init(eGrid,PDF,CDF,interType)

      end function readSingleEnergyTabularPdf


      function readNuArray() result (nuDat)
        !! Function reads total vu_bar data (sum of prompt and delayed)
        class(releaseLawENDF), pointer             :: nuDat
        integer(shortInt)                          :: addr             ! Adress of NU data
        integer(shortInt)                          :: nuType           ! Polynomial or Tabular
        real(defReal)                              :: addr_r, nuType_r ! Dummy Reals for tests
        character(100),parameter                   :: Here='readNuArray (emissionFromACE_func.f90)'

        if (JXS(2) == 0) call fatalError(Here,'Nu data does not exist: JXS(2) = 0 !')

        ! Read Adress of NU data
        addr_r = XSS(JXS(2))

        if(.not.isInteger(addr_r)) call fatalError(Here,'Value under JXS(2) is not integer!')
        addr = addr_r

        ! Choose addres of NU array
        if (addr > 0) then               ! Only one table is given (promt or total)
          addr = JXS(2)
        elseif (addr < 0) then           ! Both promt and total tabels are given. Use total.
          addr = JXS(2) + abs(addr) + 1
        else
          call fatalError(Here,'Value XSS(JXS(2)) == 0. This is undefined')
        end if

        ! Read type
        nuType_r = XSS(addr)
        if(.not.isInteger(nuType_r)) call fatalError(Here,'Type of NU array is not an integer.')
        nuType = nuType_r

        if (nuType == polynomial) then
          nuDat => readPolynomialNu(addr+1)

        else if(nuType == tabular) then
          nuDat => readTabularNu(addr+1)
        else
          call fatalError(Here,'Unrecognised type of Nu data')
        end if

      end function readNuArray



      function readPolynomialNu(addr) result (nuDat)
        integer(shortInt), intent(in)          :: addr
        class(releaseLawENDF), pointer         :: nuDat
        integer(shortInt)                      :: nCoeff
        real(defReal)                          :: nCoeff_r
        real(defReal),dimension(:),allocatable :: coeffs
        character(100),parameter               :: Here='readPolynomialNu (emissionFromACE_func.f90)'

        ! Read Number of coefficients and check it
        nCoeff_r = XSS(addr)

        if(.not.isInteger(nCoeff_r)) then
          call fatalError(Here,'Number of polynomial coefficients is not an integer')
        end if

        nCoeff = nCoeff_r

        coeffs = XSS(addr+1 : addr+1+nCoeff-1)

        nuDat => polynomialRelease(coeffs)

      end function readPolynomialNu


      function readTabularNu(addr) result (nuDat)
        integer(shortInt), intent(in)              :: addr
        class(releaseLawENDF), pointer             :: nuDat
        integer(shortInt)                          :: nEne  , nInter, i
        real(defReal)                              :: nEne_r, nInter_r
        real(defReal),dimension(:),allocatable     :: eGrid
        real(defReal),dimension(:),allocatable     :: nuValues
        integer(shortInt),dimension(:),allocatable :: bounds      ! Bounds of ENDF inter regions
        integer(shortInt),dimension(:),allocatable :: interEndf   ! Flags for ENDF interpolation
        character(100),parameter                   :: Here='readtabularNu (emissionFromACE_func.f90)'

        ! Read number of interpolation regions and energies
        nInter_r = XSS(addr)
        if(.not.isInteger(nInter_r)) call fatalError(Here,'Number of inter regions is not an int')
        if(nInter_r < 0 )            call fatalError(Here,'-ve number of inter regions')
        nInter = nInter_r

        nEne_r = XSS(addr+1+2*nInter)
        if(.not.isInteger(nEne_r)) call fatalError(Here,'Number of energy points is not an int')
        if(nEne_r < 0 )            call fatalError(Here,'-ve number of energy points')
        nEne = nEne_r

        ! Construct objects
        if (nInter == 0) then
          i = addr+2
          eGrid    = XSS(i      : i+nEne-1  )
          nuValues = XSS(i+nEne : i+2*nEne-1)

          nuDat => tabularRelease(eGrid,nuValues)

        else
          i = addr +2
          bounds    = XSS(i        : i+nInter-1  )
          interEndf = XSS(i+nInter : i+2*nInter-1)

          i = addr+2+2*nInter
          eGrid    = XSS(i      : i+nEne-1)
          nuValues = XSS(i+nEne : i+2*nEne)

          nuDat => tabularRelease(eGrid,nuValues,bounds,interEndf)

        end if

      end function readTabularNu


      function findMtIdx()
        !! Function that finds index of the given MT number in the isotope ACE library
        integer(shortInt) :: findMTIdx
        integer(shortInt) :: nMT
        real(defReal)     :: mtReal

        nMT = NXS(4)
        mtReal = real(MT,defReal)

        findMtIdx = linFind(XSS(JXS(3) : JXS(3)+nMT-1),mtReal)
        call searchError(findMtIdx,Here)

      end function findMtIdx



  end function emissionFromACE




end module emissionFromACE_func
