module materialMenu_test

  use numPrecision
  use dictionary_class,   only : dictionary
  use dictParser_func,    only : charToDict

  use materialMenu_mod,   only : init_menu => init, kill_menu => kill, nameMap, materialDefs, &
                                              display, nMat, getMatPtr, materialItem
  use funit

  implicit none

  character(*),parameter :: INPUT_STR = "     &
  mat1 { temp 273;                            &
         composition {                        &
         1001.03 12;                          &
         8016.07 0.00654;                     &
         }                                    &
         xsPath ./A_PATH;                     &
       }                                      &
  mat2 { temp 1;                              &
          composition {                       &
          100253.00 7.0E+4;                   &
          }                                   &
        }                                     "

contains

  !!
  !! Test the initialisation of matarialMenu
  !!
@Test
  subroutine testMaterialMenu()
    type(dictionary)           :: matDict
    type(dictionary)           :: emptyDict
    type(materialItem),pointer :: matPtr
    integer(shortInt)          :: i1, i2, i
    character(nameLen)         :: name
    real(defReal), parameter   :: TOL = 1.0E-6_defReal

    ! Build from empty and see if anything crashes
    call emptyDict % init(1)
    call init_menu(emptyDict)
    call kill_menu()

    ! Build some real definitions
    call charToDict(matDict, INPUT_STR)
    call init_menu(matDict)

    ! Check that materials are present in the nameMap
    name = 'mat1'
    i1 = nameMap % getOrDefault(name,0)
    @assertTrue(i1 == 1 .or. i1 == 2)

    name = 'mat2'
    i2 = nameMap % getOrDefault(name,0)
    @assertTrue(i2 == 1 .or. i2 == 2)

    ! Verify definitions

    ! mat 1
    @assertEqual(273.0_defReal, materialDefs(i1) % T, TOL)
    @assertEqual(i1, materialDefs(i1) % matIdx)

    call materialDefs(i1) % extraInfo % get(name,'xsPath')
    @assertEqual('./A_PATH', trim(name))

    ! Check individual compositions
    do i=1,2
      if(materialDefs(i1) % nuclides(i) % Z == 1) then
        @assertEqual(1, materialDefs(i1) % nuclides(i) % A)
        @assertEqual(3, materialDefs(i1) % nuclides(i) % T)
        @assertEqual(12.0_defReal, materialDefs(i1) % dens(i), TOL*12.0_defReal)
      else if(materialDefs(i1) % nuclides(i) % Z == 8) then
        @assertEqual(16, materialDefs(i1) % nuclides(i) % A)
        @assertEqual(7, materialDefs(i1) % nuclides(i) % T)
        @assertEqual(0.00654_defReal, materialDefs(i1) % dens(i), TOL*0.00654_defReal)
      else
        @assertTrue(.false., 'Error when reading Mat 1 compositions')
      end if
    end do

    ! mat 2
    @assertEqual(1.0_defReal, materialDefs(i2) % T, TOL)
    @assertEqual(i2, materialDefs(i2) % matIdx)
    @assertEqual(7.0E+4_defReal, materialDefs(i2) % dens(1), TOL*7.0E+4_defReal)
    @assertEqual(100, materialDefs(i2) % nuclides(1) % Z)
    @assertEqual(253, materialDefs(i2) % nuclides(1) % A)
    @assertEqual(0,   materialDefs(i2) % nuclides(1) % T)

    ! Get number of materials
    @assertEqual(2, nMat())

    ! Get pointer to material1
    matPtr => getMatPtr(1)
    @assertEqual('mat1', trim(matPtr % name))

    ! Test writing nuclide info back to definition string
    @assertEqual('1001.03', matPtr % nuclides(1) % toChar())

    matPtr => getMatPtr(2)
    @assertEqual('100253.00', matPtr % nuclides(1) % toChar())

    ! Clean
    call kill_menu()

  end subroutine testMaterialMenu


end module materialMenu_test
