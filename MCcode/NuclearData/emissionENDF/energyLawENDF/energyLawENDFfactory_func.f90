module energyLawENDFfactory_func

  use numPrecision
  use endfConstants
  use genericProcedures,  only : fatalError, numToChar
  use aceCard_class,      only : aceCard

  ! Energy Laws
  use energyLawENDF_inter,       only : energyLawENDF
  use contTabularEnergy_class,   only : contTabularEnergy
  use maxwellSpectrum_class,     only : maxwellSpectrum
  use evaporationSpectrum_class, only : evaporationSpectrum
  use levelScattering_class,     only : levelScattering
  use noEnergy_class,            only : noEnergy
  use multipleEnergyLaws_class,  only : multipleEnergyLaws

  implicit none
  private

  ! Define public interface
  public :: new_energyLawENDF

contains

  !!
  !! Allocates a new energyLawENDF from aceCard and MT number
  !! aceCard can be in any poistion. Its position changes at output.
  !!
  subroutine new_energyLawENDF(new, ACE, MT)
    class(energyLawENDF),allocatable, intent(inout) :: new
    type(aceCard), intent(inout)                    :: ACE
    integer(shortInt), intent(in)                   :: MT
    integer(shortInt)                               :: LNW, LAW, loc
    integer(shortInt)                               :: N, i, NR, NEne
    class(multipleEnergyLaws),allocatable           :: multiLaw
    character(100),parameter :: Here='new_energyLawENDF (energyLawENDFfactory_func.f90)'

    ! Protect against reactions with no energy data (absorbtions and eleastic scattering)
    if(MT == N_N_elastic) then
      allocate(new, source = noEnergy())
      return

    else if(ACE % isCaptureMT(MT)) then
      allocate(new, source = noEnergy())
      return

    end if

    ! Set aceCard read head to beginning of energy data for MT
    call ACE % setToEnergyMT(MT)

    ! Read location of next energy law. If LNW == 0 only one law is given
    LNW = ACE % readInt()

    ! Build energy law
    if (LNW == 0) then ! Build single energy law
      LAW = ACE % readInt()
      loc = ACE % readInt()
      call buildENDFLaw(new, LAW, loc, ACE)

    else ! Build multiple energy laws

      ! Calculate number of energy laws
      N = 1
      do while (LNW /= 0 )
        N = N + 1
        call ACE % setToEnergyLaw(LNW)
        LNW = ACE % readINT()
        if (N > 100) call fatalError(Here, 'Infinate Loop. Terminating.')
      end do

      ! Allocate space
      allocate(multiLaw)
      call multiLaw % init(N)

      ! Reset to the begining
      call ACE % setToEnergyMT(MT)
      LNW = ACE % readInt()

      ! Sequentialy read all laws
      do i=1,N
        LAW = ACE % readInt()
        loc = ACE % readInt()

        ! Read number of regions and interpolation parameters
        NR = ACE % readInt()
        associate(bounds => ACE % readIntArray(NR), interENDF => ACE % readIntArray(NR))
          ! Read number of energy regions
          NEne = ACE % readInt()
          associate( eGrid => ACE % readRealArray(NEne), pdf => ACE % readRealArray(NEne))
            ! Build energy law
            call buildENDFLaw(new, LAW, loc, ACE)

            if(NR == 0) then
              call multiLaw % addLaw(new, eGrid, pdf)
            else
              call multiLaw % addLaw(new, eGrid, pdf, bounds, interENDF)
            end if
          end associate
        end associate

        ! Move to the next low and read new LNW
        ! Protect against last iteration
        if (LNW /= 0) then
          call ACE % setToEnergyLaw(LNW)
          LNW = ACE % readInt()
        end if
      end do

      ! Verify that all laws were read
      if(LNW /= 0) call fatalError(Here,'LNW is not 0 after reading all energy laws. It is ' // &
                                         numToChar(LNW) //' Somthing failed')

      ! Move finished multiple laws to new
      call move_alloc(multiLaw, new)
    end if

  end subroutine new_energyLawENDF


  !!
  !! Helper function to allocate an energy law
  !!
  !! See MCNP 4 Manual Appendix F Table F-14 for extra info
  !!
  !! Args:
  !!   lawENDF [inout] -> polymorphic energyLawENDF to be allocated
  !!   LAW [in]    -> Integer indentifier of the ENDF law type
  !!   loc [in]    -> Location of the law data relative to JXS(11)
  !!   ACE [inout] -> ACE Card with head set to the begining of data for law.
  !!
  !! Errors:
  !!   Will crash if loc is set to incorrect location
  !!
  subroutine buildENDFLaw(lawENDF, LAW, loc, ACE)
    class(energyLawENDF),allocatable, intent(inout) :: lawENDF
    integer(shortInt), intent(in)                   :: LAW
    integer(shortInt), intent(in)                   :: loc
    type(aceCard), intent(inout)                    :: ACE
    character(100),parameter :: Here='buildENDFLaw (energyLawENDFfactory_func.f90)'

    ! Deallocate lawENDF if allocated
    if(allocated(lawENDF)) deallocate(lawENDF)

    ! Build approperiate energy Law
    call ACE % setToEnergyLaw(loc)

    ! Allocate new object
    select case(LAW)
      case (continuousTabularDistribution)
        allocate(lawENDF, source = contTabularEnergy(ACE))

      case (simpleMaxwellFissionSpectrum)
        allocate(lawENDF, source = maxwellSpectrum(ACE))

      case (levelScatteringLaw)
        allocate(lawENDF, source = levelScattering(ACE))

      case (evaporationEnergySpectrum)
        allocate(lawENDF, source = evaporationSpectrum(ACE))

      case default
        print *, 'Energy Law Type :', LAW
        call fatalError(Here,'Energy Law Type is not recognised, yet supported or is correlated')

    end select

  end subroutine buildENDFLaw

end module energyLawENDFfactory_func
