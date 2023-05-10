module poly_func_test
  use numPrecision
  use poly_func, only : poly_integrate, poly_eval, poly_solve
  use pFUnit_mod

  implicit none

contains

@Test
  subroutine testPoly()
    real(defReal), dimension(6) :: poly1, poly2, poly3, poly4
    real(defReal)               :: x1, x2
    real(defReal)               :: tol = 1.0E-4_defReal

    ! Test array
    poly1(:) = [23.20_defReal, 1.59_defReal, 0.12_defReal, &
                10.02_defReal, 0.06_defReal, 8.03_defReal ]

    ! Analytical integral
    poly2(:) = [2.10526_defReal, 1.5_defReal,  0.013289_defReal, &
                  11.02_defReal, 1.06_defReal, 9.03_defReal      ]

    ! Integrate using poly_integrate
    poly3 = poly1
    call poly_integrate(poly3)

    ! Check
    @assertEqual(poly3,poly2,tol)

    ! Evaluate using poly_eval
    x1 = poly_eval(poly1, 1.074_defReal)

    ! Check
    @assertEqual(x1,49.2504_defReal,tol)

    ! Solve poly2 = 1 using poly_solve, using x0 = 10
    x2 = poly_solve(poly2, poly1, 10.0_defReal, 1.0_defReal)

    ! Check
    @assertEqual(x2,0.66643669_defReal,tol)

  end subroutine testPoly


end module poly_func_test
