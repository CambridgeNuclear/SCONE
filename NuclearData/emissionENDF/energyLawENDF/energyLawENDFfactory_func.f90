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
  !! Allocates a new energyLawENDF from aceCard
  !!
  !! Can be used to build energy distribution for an MT reaction
  !! And for a fission from a delayed precursor group
  !!
  !! Args:
  !!   ACE [inout]  -> ACE Card with the data
  !!   MT [in]      -> MT number. If delayed = .true. then it is precursor group
  !!   delayed [in] -> Optional. If .true. build distribution for a precursor group
  !!
  !! Errors:
  !!   Allocated to noEnergy if MT is Elastic Scttering or Capture reaction
  !!
  subroutine new_energyLawENDF(new, ACE, MT, delayed)
    class(energyLawENDF),allocatable, intent(inout) :: new
    type(aceCard), intent(inout)                    :: ACE
    integer(shortInt), intent(in)                   :: MT
    logical(defBool),optional,intent(in)            :: delayed
    logical(defBool)                                :: del_loc
    integer(shortInt)                               :: LNW, LAW, loc
    integer(shortInt)                               :: N, i, NR, NEne
    integer(shortInt)                               :: root, LOCC
    class(multipleEnergyLaws),allocatable           :: multiLaw
    character(100),parameter :: Here='new_energyLawENDF (energyLawENDFfactory_func.f90)'

    ! Set default value
    if(present(delayed)) then
      del_loc = delayed
    else
      del_loc = .false.
    end if

    ! Set approperiate root and initial offset for the energy law
    if(del_loc) then
      root = ACE % getRootAddress('energyLawsPrecursors')
      LOCC = ACE % LOCCforPrecursor(MT)

    else
      ! Protect against reactions with no energy data (absorbtions and eleastic scattering)
      if(MT == N_N_elastic) then
        allocate(new, source = noEnergy())
        return

      else if(ACE % isCaptureMT(MT)) then
        allocate(new, source = noEnergy())
        return

      end if
      ! Set aceCard read head to beginning of energy data for MT
      root = ACE % getRootAddress('energyLawsMT')
      LOCC = ACE % LOCCforMT(MT)

    end if

    call ACE % setRelativeTo(root, LOCC)


    ! Read location of next energy law. If LNW == 0 only one law is given
    LNW = ACE % readInt()

    ! Build energy law
    if (LNW == 0) then ! Build single energy law
      LAW = ACE % readInt()
      loc = ACE % readInt()
      call buildENDFLaw(new, LAW, root, loc, ACE)

    else ! Build multiple energy laws

      ! Calculate number of energy laws
      N = 1
      do while (LNW /= 0 )
        N = N + 1
        call ACE % setRelativeTo(root, LNW)
        LNW = ACE % readINT()
        if (N > 100) call fatalError(Here, 'Infinate Loop. Terminating.')
      end do

      ! Allocate space
      allocate(multiLaw)
      call multiLaw % init(N)

      ! Reset to the begining
      call ACE % setRelativeTo(root, LOCC)
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
            call buildENDFLaw(new, LAW, root, loc, ACE)

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
          call ACE % setRelativeTo(root, LNW)
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
  !!   root [in]   -> Root location for this energy law data block (e.g. JXS(11))
  !!   offset [in] -> Location of the energy law relative to root
  !!   ACE [inout] -> ACE Card with data
  !!
  !! Errors:
  !!   Will crash if root & offset point to an incorrect location
  !!
  subroutine buildENDFLaw(lawENDF, LAW, root, offset, ACE)
    class(energyLawENDF),allocatable, intent(inout) :: lawENDF
    integer(shortInt), intent(in)                   :: LAW
    integer(shortInt), intent(in)                   :: root
    integer(shortInt), intent(in)                   :: offset
    type(aceCard), intent(inout)                    :: ACE
    character(100),parameter :: Here = 'buildENDFLaw (energyLawENDFfactory_func.f90)'

    ! Deallocate lawENDF if allocated
    if(allocated(lawENDF)) deallocate(lawENDF)

    ! Build approperiate energy Law
    call ACE % setRelativeTo(root, offset)

    ! Allocate new object
    select case(LAW)
      case (continuousTabularDistribution)
        allocate(lawENDF, source = contTabularEnergy(ACE, root))

      case (simpleMaxwellFissionSpectrum)
        allocate(lawENDF, source = maxwellSpectrum(ACE))

      case (levelScatteringLaw)
        allocate(lawENDF, source = levelScattering(ACE))

      case (evaporationEnergySpectrum)
        allocate(lawENDF, source = evaporationSpectrum(ACE))

      case default
        call fatalError(Here,'Energy law type is not recognised or yet &
                              &supported: '//numToChar(LAW))

    end select

  end subroutine buildENDFLaw

end module energyLawENDFfactory_func
