module poly_func

  use universalVariables
  use numPrecision
  use genericProcedures

  implicit none

  !! Module to perform operations on simple polynomials
  !! Polynomials are given as a 1D array containing n coefficients followed by n exponents
  !!
  !! Interface:
  !!   poly_integrate -> Update exponents and coefficients to integrate polynomial
  !!                        - Currently gives error for input exponent of -1
  !!   poly_solve     -> Solves polynomial using Newton-Raphson method
  !!

  contains

    !!
    !! Integrates a simple polynomial given coefficients and exponents,
    !!  returning indefinite integral
    !!
    !! Args:
    !!   equation -> 1D array of n coefficients followed by n exponents
    !!
    !! Errors:
    !!   Input array size is not divisible by 2
    !!   Exponent of -1 is integrated
    !!
    subroutine poly_integrate(equation)
      real(defReal), dimension(:), intent(inout) :: equation
      integer(shortInt)                          :: n, i
      character(100), parameter                  :: Here = "poly_integrate (poly_func.f90)"

      ! Check that array is of even size
      if( modulo(size(equation), 2) /= 0 ) then
        call fatalError(Here, "Array size must be divisible by 2")
      end if

      n = size(equation) / 2

      ! Integrate
      do i=1, n
        ! Update exponents
        equation(i+n) = equation(i+n) + 1
        if( equation(i+n) == 0 ) then
          call fatalError(Here, "Integrating exponent of -1, currently not yet supported")
        end if
        ! Update coefficients
        equation(i) = equation(i) / equation(i+n)
      end do

    end subroutine poly_integrate


    !!
    !! Use Newton-Raphspon method to solve polynomial with m terms
    !!
    !! Args:
    !!   equation ->
    !!   derivative ->
    !!   x0 ->
    !!   const -> For f(x) = const, if not given then solves f(x) = 0
    !!
    !! Errors:
    !!   equation and derivative are different sizes
    !!   Input array sizes are not divisible by 2
    !!
    function poly_solve(equation, derivative, x0, const) result(x)
      real(defReal), dimension(:), intent(in) :: equation
      real(defReal), dimension(:), intent(in) :: derivative
      real(defReal), intent(in)               :: x0
      real(defReal), intent(in), optional     :: const
      real(defReal)                           :: x, x_old, f, f_dash, c, tol
      integer(shortInt)                       :: i, j, m
      character(100), parameter               :: Here = "poly_solve (poly_func.f90)"

      ! Check that input arrays are of same size
      if( size(equation) /= size(derivative) ) then
        call fatalError(Here, "Equation and Derivative array inputs are of different sizes")
      end if

      ! Check that array is of even size
      if( modulo(size(equation), 2) /= 0 ) then
        call fatalError(Here, "Array size must be divisible by 2")
      end if

      x = x0
      m = size(equation) / 2

      ! If no constant present then solving f(x) = 0
      if(  present(const) ) then
        c = const
      else
        c = 0
      end if

      ! Iterate
      i = 0
      iterate: do
        ! Store x for convergence check
        x_old = x

        ! Calculate f(x) and f'(x)
        f = -c
        f_dash = 0
        do j=1,m
          f = f + equation(j) * x ** equation(j+m)
          f_dash = f_dash + derivative(j) * x ** derivative(j+m)
        end do

        ! Update x
        x = x - f / f_dash

        ! Check for convergence
        tol = 0.0000000001
        if( x > (1-tol)*x_old .and. x < (1+tol)*x_old ) exit iterate

        ! Call error if not converged
        if( i >= 1000 ) then
          call fatalError(Here, "Solution has not converged after 1000 iterations")
        end if

        ! Increase counter
        i = i+1

      end do iterate

    end function poly_solve

    !!
    !! Gives output value y for y = f(x)
    !!
    !! Args:
    !!   f -> Array defining polynomial
    !!   x -> Point at which to evaluate
    !!
    function poly_eval(f, x) result(y)
      real(defReal), dimension(:), intent(in) :: f
      real(defReal), intent(in)               :: x
      real(defReal)                           :: y
      integer(shortInt)                       :: n, i
      character(100), parameter               :: Here = "poly_eval (poly_func.f90)"

      ! Check that array is of even size
      if( modulo(size(f), 2) /= 0 ) then
        call fatalError(Here, "Array size must be divisible by 2")
      end if

      n = size(f) / 2

      y = 0
      do i=1,n
        y = y + f(i) * x ** f(i+n)
      end do

    end function poly_eval

end module poly_func
