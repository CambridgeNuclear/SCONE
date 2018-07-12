module vector_class

  use numPrecision

  implicit none
  private

  !!
  !! Type that contains a vector in R^3
  !! It exists so vector is a scalar to Fortran and emelental procedures
  !! with vector arguments can be defined.
  !!
  !! All basic arithmetic operations with scalars are defined)
  !!   Vector can be multiplied and divided by real kind 4&8 and int kind 4&8
  !!   Vector can be added with vector or dimension(3) real kind 4&8
  !!   Vector or dimension(3) real kind 4&8 can be substructed from vector to give vector
  !!   Vector can be assigned from dimension(3) array of real kind 4&8
  !!   Vector is initialised to all ZERO by default
  !!
  type, public :: vector
    real(defReal),dimension(3) :: v = [ZERO, ZERO, ZERO]
  contains
    ! Overloading of normal assignment and arthmetic operators
    generic            :: operator(+)       => addVector_vector   ,&
                                               addVector_RHS_real4,&
                                               addVector_LHS_real4,&
                                               addVector_RHS_real8,&
                                               addVector_LHS_real8

    generic            :: operator(-)       => subVector_vector   ,&
                                               subVector_real4    ,&
                                               subVector_real8

    generic            :: assignment(=)     => assignVector_real8 ,&
                                               assignVector_real4

    generic            :: operator(*)       => scalarMulti_RHS_real8, &
                                               scalarMulti_LHS_real8, &
                                               scalarMulti_RHS_real4, &
                                               scalarMulti_LHS_real4, &
                                               scalarMulti_RHS_int4, &
                                               scalarMulti_LHS_int4, &
                                               scalarMulti_RHS_int8, &
                                               scalarMulti_LHS_int8

    generic            :: operator(/)       => scalarDiv_real8,&
                                               scalarDiv_real4,&
                                               scalarDiv_int4,&
                                               scalarDiv_int8
    ! Vector spefific operators
    generic            :: operator(.dot.)   => dotProduct
    generic            :: operator(.cross.) => crossProduct
    procedure          :: L2norm

    ! Assignment and normal arthmetic functions for diffrent real and int kinds
    procedure, private           :: addVector_vector
    procedure,private            :: addVector_RHS_real4
    procedure,pass(RHS),private  :: addVector_LHS_real4
    procedure,private            :: addVector_RHS_real8
    procedure,pass(RHS),private  :: addVector_LHS_real8
    procedure, private           :: subVector_vector
    procedure,private            :: subVector_real4
    procedure,private            :: subVector_real8
    procedure, private           :: scalarMulti_RHS_real8
    procedure,pass(RHS),private  :: scalarMulti_LHS_real8
    procedure, private           :: scalarMulti_RHS_real4
    procedure,pass(RHS),private  :: scalarMulti_LHS_real4
    procedure, private           :: scalarMulti_RHS_int4
    procedure,pass(RHS),private  :: scalarMulti_LHS_int4
    procedure, private           :: scalarMulti_RHS_int8
    procedure,pass(RHS),private  :: scalarMulti_LHS_int8
    procedure,private            :: scalarDiv_real8
    procedure,private            :: scalarDiv_real4
    procedure,private            :: scalarDiv_int4
    procedure,private            :: scalarDiv_int8
    procedure, private           :: assignVector_real8
    procedure, private           :: assignVector_real4

    ! Vector specific functions
    procedure, private           :: dotProduct
    procedure, private           :: crossProduct
  end type vector

contains

  !!
  !! Add two vectors together
  !!
  elemental function addVector_vector(LHS,RHS) result(res)
    class(vector), intent(in) :: LHS
    type(vector), intent(in)  :: RHS
    type(vector)              :: res

    res % v = LHS % v + RHS % v

  end function addVector_vector

  !!
  !! Add dimension(3) array of real8 from RHS to vector
  !!
  pure function addVector_RHS_real8(LHS,RHS) result(res)
    class(vector), intent(in)         :: LHS
    real(8), dimension(3), intent(in) :: RHS
    type(vector)                      :: res

    res % v = RHS + LHS % v

  end function addVector_RHS_real8

  !!
  !! Add dimension(3) array of real8 from LHS to vector
  !!
  pure function addVector_LHS_real8(LHS,RHS) result(res)
    real(8), dimension(3), intent(in) :: LHS
    class(vector), intent(in)         :: RHS
    type(vector)                      :: res

    res % v = RHS % v + LHS

  end function addVector_LHS_real8

  !!
  !! Add dimension(3) array of real4 from RHS to vector
  !!
  pure function addVector_RHS_real4(LHS,RHS) result(res)
    class(vector), intent(in)         :: LHS
    real(4), dimension(3), intent(in) :: RHS
    type(vector)                      :: res

    res % v = RHS + LHS % v

  end function addVector_RHS_real4

  !!
  !! Add dimension(3) array of real4 from LHS to vector
  !!
  pure function addVector_LHS_real4(LHS,RHS) result(res)
    real(4), dimension(3), intent(in) :: LHS
    class(vector), intent(in)         :: RHS
    type(vector)                      :: res

    res % v = RHS % v + LHS

  end function addVector_LHS_real4

  !!
  !! Substract two vectors
  !!
  elemental function subVector_vector(LHS,RHS) result(res)
    class(vector), intent(in) :: LHS
    type(vector), intent(in)  :: RHS
    type(vector)              :: res

    res % v = LHS % v - RHS % v

  end function subVector_vector

  !!
  !! Substract dimension(3) array of real8 from vector
  !!
  pure function subVector_real8(LHS,RHS) result(res)
    class(vector), intent(in)         :: LHS
    real(8), dimension(3), intent(in) :: RHS
    type(vector)                      :: res

    res % v = LHS % v - RHS

  end function subVector_real8

  !!
  !! Substract dimension(3) array of real4 from vector
  !!
  pure function subVector_real4(LHS,RHS) result(res)
    class(vector), intent(in)         :: LHS
    real(4), dimension(3), intent(in) :: RHS
    type(vector)                      :: res

    res % v = LHS % v - RHS

  end function subVector_real4

  !!
  !! Multiplication of vector by scalar real8 from RHS
  !!
  elemental function scalarMulti_RHS_real8(LHS,RHS) result(res)
    class(vector), intent(in) :: LHS
    real(8), intent(in)       :: RHS
    type(vector)              :: res

    res % v = LHS % v * RHS

  end function scalarMulti_RHS_real8

  !!
  !! Multiplication of vector by scalar real8 from LHS
  !!
  elemental function scalarMulti_LHS_real8(LHS,RHS) result(res)
    real(8), intent(in)       :: LHS
    class(vector), intent(in) :: RHS
    type(vector)              :: res

    res % v = RHS % v * LHS

  end function scalarMulti_LHS_real8

  !!
  !! Multiplication of vector by scalar real from RHS
  !!
  elemental function scalarMulti_RHS_real4(LHS,RHS) result(res)
    class(vector), intent(in) :: LHS
    real(4), intent(in)       :: RHS
    type(vector)              :: res

    res % v = LHS % v * RHS

  end function scalarMulti_RHS_real4

  !!
  !! Multiplication of vector by scalar real from LHS
  !!
  elemental function scalarMulti_LHS_real4(LHS,RHS) result(res)
    real(4),intent(in)        :: LHS
    class(vector), intent(in) :: RHS
    type(vector)              :: res

    res % v = RHS % v * LHS

  end function scalarMulti_LHS_real4

  !!
  !! Multiplication of vector by scalar int4 from RHS
  !!
  elemental function scalarMulti_RHS_int4(LHS,RHS) result(res)
    class(vector), intent(in)     :: LHS
    integer(4), intent(in)        :: RHS
    type(vector)                  :: res

    res % v = LHS % v * RHS

  end function scalarMulti_RHS_int4

  !!
  !! Multiplication of vector by scalar int4 from LHS
  !!
  elemental function scalarMulti_LHS_int4(LHS,RHS) result(res)
    integer(4),intent(in)       :: LHS
    class(vector), intent(in)   :: RHS
    type(vector)                :: res

    res % v = RHS % v * LHS

  end function scalarMulti_LHS_int4

  !!
  !! Multiplication of vector by scalar int8 from RHS
  !!
  elemental function scalarMulti_RHS_int8(LHS,RHS) result(res)
    class(vector), intent(in)     :: LHS
    integer(8), intent(in)        :: RHS
    type(vector)                  :: res

    res % v = LHS % v * RHS

  end function scalarMulti_RHS_int8

  !!
  !! Multiplication of vector by scalar int8 from LHS
  !!
  elemental function scalarMulti_LHS_int8(LHS,RHS) result(res)
    integer(8),intent(in) :: LHS
    class(vector), intent(in)    :: RHS
    type(vector)                 :: res

    res % v = RHS % v * LHS

  end function scalarMulti_LHS_int8


  !!
  !! Division of vector by scalar real8
  !!
  elemental function scalarDiv_real8(LHS,RHS) result(res)
    class(vector), intent(in) :: LHS
    real(8), intent(in)       :: RHS
    type(vector)              :: res

    res % v = LHS % v / RHS

   end function scalarDiv_real8

  !!
  !! Division of vector by scalar real4
  !!
  elemental function scalarDiv_real4(LHS,RHS) result(res)
    class(vector), intent(in) :: LHS
    real(4), intent(in)       :: RHS
    type(vector)              :: res

    res % v = LHS % v / RHS

  end function scalarDiv_real4

  !!
  !! Division of vector by scalar int4
  !!
  elemental function scalarDiv_int4(LHS,RHS) result(res)
    class(vector), intent(in) :: LHS
    integer(4), intent(in)    :: RHS
    type(vector)              :: res

    res % v = LHS % v / RHS

  end function scalarDiv_int4

  !!
  !! Division of vector by scalar int8
  !!
  elemental function scalarDiv_int8(LHS,RHS) result(res)
    class(vector), intent(in) :: LHS
    integer(8), intent(in)    :: RHS
    type(vector)              :: res

    res % v = LHS % v / RHS

  end function scalarDiv_int8

  !!
  !! Assign vector from size(3) array of real8s
  !!
  subroutine assignVector_real8(LHS,RHS)
    class(vector), intent(out)            :: LHS
    real(8),dimension(3),intent(in)       :: RHS

    LHS % v = RHS

  end subroutine assignVector_real8

  !!
  !! Assign vector from size(3) array of real4s
  !!
  subroutine assignVector_real4(LHS,RHS)
    class(vector), intent(out)            :: LHS
    real(4),dimension(3),intent(in)       :: RHS

    LHS % v = RHS

  end subroutine assignVector_real4

  !!
  !! Dot product
  !!
  elemental function dotProduct(LHS,RHS) result(res)
    class(vector), intent(in)  :: LHS
    type(vector), intent(in)   :: RHS
    real(defReal)              :: res
    real(defReal),dimension(3) :: temp

    temp = LHS % v * RHS % v
    res = temp(1) + temp(2) + temp(3)

  end function dotProduct

  !!
  !! Cross product
  !!
  elemental function crossProduct(LHS,RHS) result(res)
    class(vector), intent(in) :: LHS
    type(vector), intent(in)  :: RHS
    type(vector)              :: res

    res % v(1) = LHS % v(2) * RHS % v(3) - LHS % v(3) * RHS % v(2)
    res % v(2) = LHS % v(3) * RHS % v(1) - LHS % v(1) * RHS % v(3)
    res % v(3) = LHS % v(1) * RHS % v(2) - LHS % v(2) * RHS % v(1)

  end function crossProduct

  !!
  !! L-2 Norm of the vector (length
  !!
  elemental function L2norm(self) result(res)
    class(vector), intent(in) :: self
    real(defReal)             :: res

    res = sqrt( self % v(1) * self % v(1) + &
                self % v(2) * self % v(2) + &
                self % v(3) * self % v(3)   )

  end function L2norm

end module vector_class
