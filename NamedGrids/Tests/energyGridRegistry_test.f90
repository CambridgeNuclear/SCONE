module energyGridRegistry_test

  use numPrecision
  use dictionary_class,       only : dictionary
  use energyGrid_class,       only : energyGrid
  use energyGridRegistry_mod, only : define_energyGrid, &
                                     define_multipleEnergyGrids, &
                                     get_energyGrid, &
                                     kill_energyGrids
  use funit

  implicit none


contains

  !!
  !! Before each test
  !!
@before
  subroutine init()
    type(dictionary)   :: dict_lin
    type(dictionary)   :: dict_log
    type(dictionary)   :: dict_unstruct
    type(dictionary)   :: dict_combined
    character(nameLen) :: name

    ! Define linear grid
    call dict_lin % init(4)
    call dict_lin % store('grid','lin')
    call dict_lin % store('min',0.0_defReal)
    call dict_lin % store('max',20.0_defReal)
    call dict_lin % store('size',20)

    ! Define logarithmic grid
    call dict_log % init(4)
    call dict_log % store('grid','log')
    call dict_log % store('min',1.0E-9_defReal)
    call dict_log % store('max',1.0_defReal)
    call dict_log % store('size',9)

    ! Define unstructured grid
    call dict_unstruct % init(2)
    call dict_unstruct % store('grid','unstruct')
    call dict_unstruct % store('bins',[10.0_defReal, 1.0_defReal, 1.0E-6_defReal])

    ! Define combined dictionary
    call dict_combined % init(3)
    call dict_combined % store('linGridC',dict_lin)
    call dict_combined % store('logGridC',dict_log)
    call dict_combined % store('unstructGridC',dict_unstruct)

    ! Build definitions
    name = 'linGrid'
    call define_energyGrid(name, dict_lin)

    name = 'logGrid'
    call define_energyGrid(name, dict_log)

    name = 'unstructGrid'
    call define_energyGrid(name, dict_unstruct)

    call define_multipleEnergyGrids(dict_combined)

  end subroutine init

  !!
  !! After each test
  !!
@after
  subroutine done()

    call kill_energyGrids()

  end subroutine done

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Proper tests begin here
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test getting user-defined grid and check that bins are as expected
  !!
@Test
  subroutine testGettingUserDefinedGrid()
    type(energyGrid)         :: myGrid
    character(nameLen)       :: name
    real(defReal), parameter :: TOL = 1.0E-9

    ! Get linear grid
    name = 'linGrid'
    call get_energyGrid(myGrid, name)
    @assertEqual(18.0_defReal, myGrid % bin(3), 18.0_defReal * TOL )

    ! Get logarithmic grid
    name = 'logGrid'
    call get_energyGrid(myGrid, name)
    @assertEqual(1.0E-2_defReal, myGrid % bin(3), 1.0E-2_defReal * TOL )

    ! Get unstruct grid
    name = 'unstructGrid'
    call get_energyGrid(myGrid, name)
    @assertEqual(1.0_defReal, myGrid % bin(2), 1.0_defReal * TOL )

    ! Get linear grid from combined dict
    name = 'linGridC'
    call get_energyGrid(myGrid, name)
    @assertEqual(18.0_defReal, myGrid % bin(3), 18.0_defReal * TOL )

    ! Get logarithmic grid from combined dict
    name = 'logGridC'
    call get_energyGrid(myGrid, name)
    @assertEqual(1.0E-2_defReal, myGrid % bin(3), 1.0E-2_defReal * TOL )

    ! Get unstruct grid from combined dict
    name = 'unstructGridC'
    call get_energyGrid(myGrid, name)
    @assertEqual(1.0_defReal, myGrid % bin(2), 1.0_defReal * TOL )


  end subroutine testGettingUserDefinedGrid

  !!
  !! Test getting predefined grids
  !!
@Test
  subroutine testGettingPreDefined()
    type(energyGrid)         :: myGrid
    character(nameLen)       :: name
    real(defReal), parameter :: TOL = 1.0E-9

    ! WIMS 69
    name = 'wims69'
    call get_energyGrid(myGrid, name)
    @assertEqual(6.7340000000E-02_defReal, myGrid % bin(11), 6.7340000000E-02_defReal * TOL)

    ! WIMS 172
    name = 'wims172'
    call get_energyGrid(myGrid, name)
    @assertEqual(7.5239800000E-06_defReal, myGrid % bin(88), 7.5239800000E-06_defReal * TOL)

    ! CASMO 40
    name = 'casmo40'
    call get_energyGrid(myGrid, name)
    @assertEqual(1.1500000000E-06_defReal, myGrid % bin(23), 1.1500000000E-06_defReal * TOL)

  end subroutine testGettingPreDefined

  !!
  !! Test getting defined and undefined grid with error flag
  !!
@Test
  subroutine testGettingUndefinedGrid()
    type(energyGrid)         :: myGrid
    logical(defBool)         :: errV
    character(nameLen)       :: name

    ! Non-Existing grid
    name = 'invalidGrid'
    call get_energyGrid(myGrid,name, err = errV)
    @assertTrue(errV)

    ! User-defined grid
    name = 'linGrid'
    call get_energyGrid(myGrid,name, err = errV)
    @assertFalse(errV)

    ! Pre-defined grid
    name = 'wims172'
    call get_energyGrid(myGrid,name, err = errV)
    @assertFalse(errV)

  end subroutine testGettingUndefinedGrid

    
end module energyGridRegistry_test
