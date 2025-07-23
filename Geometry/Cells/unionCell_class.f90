module unionCell_class

  use numPrecision
  use universalVariables, only : INF
  use genericProcedures,  only : fatalError, hasDuplicates, numToChar, charToInt
  use dictionary_class,   only : dictionary
  use surfaceShelf_class, only : surfaceShelf
  use surface_inter,      only : surface
  use cell_inter,         only : cell, surfInfo, kill_super => kill

  implicit none
  private

  ! Define CSG operators
  integer(shortInt), parameter :: union     = huge(shortInt) - 1, &
                                  intersect = union - 1, &
                                  openBra   = union - 2, &
                                  closeBra  = union - 3, &
                                  complem   = union - 4


  !!
  !! CSG cell defined by intersections, unions, and complements of halfspaces
  !! This cell is a generalisation of the simpleCell.
  !!
  !! The implementation follows that given in the paper by Romano et al. (2024).
  !! "Point containment algorithms for constructive solid geometry with
  !!  unbounded primitives"
  !!
  !! It is equipped with its own syntax to perform parentheses and to
  !! denote a union and complement.
  !! These are:
  !!   union = :
  !!   open parenthesis = <
  !!   close parenthesis = >
  !!   complement = #
  !! These are used as entries in the surface list
  !! As with the simple cell, intersections are implicit.
  !!
  !! Sample input dictionary:
  !!   cell {
  !!     type unionCell;
  !!     id 17;
  !!     surfaces [ 3 # < -1 : 2 > ];
  !!     <content info as defined in cellShelf_class >
  !!   }
  !! NOTE: surfaces is not a list but a tokenArray, using [ ], not ( )
  !!
  !! Private Members:
  !!   surfaces  -> Array that contains, halfspace info, surfIdx and a pointer to the surface.
  !!     Negative halfspace is indicated by surfIdx < 0. Positive by surfIdx > 0
  !!   evalArray -> Array containing surfaces and operations necessary to evaluate whether
  !!                a point lies inside the cell
  !!   isSimple  -> Logical identifying whether the cell is only composed of intersections.
  !!                If so, simplifies evaluation logic.
  !!
  !! Interface:
  !!   cell interface
  !!
  type, public, extends(cell) :: unionCell
    private
    type(surfInfo), dimension(:), allocatable    :: surfaces
    integer(shortInt), dimension(:), allocatable :: evalArray
    logical(defBool)                             :: isSimple = .false.
  contains
    procedure :: init
    procedure :: deMorgan
    procedure :: checkPrecedence
    procedure :: inside
    procedure :: insideComplex
    procedure :: shortCircuit
    procedure :: distance
    procedure :: kill
  end type unionCell

contains

  !!
  !! Initialise cell
  !!
  !! See cell_inter for details
  !!
  subroutine init(self, dict, surfs)
    class(unionCell), intent(inout)   :: self
    class(dictionary), intent(in)     :: dict
    type(surfaceShelf), intent(inout) :: surfs
    logical(defBool)                  :: notInt, surfacePrev, closePrev, hasUnion
    integer(shortInt)                 :: surfIdx, id, i, evalSize, num, sCount, depth
    character(nameLen), dimension(:), allocatable :: definition
    integer(shortInt), dimension(:), allocatable  :: tempEval
    type(surfInfo), dimension(:), allocatable     :: tempSurf
    character(100), parameter :: Here = 'init (unionCell_class.f90)'

    ! Get surfaces and operations from dictionary
    call dict % get(definition, 'surfaces')

    ! Allocate initial arrays with a worst-case size
    allocate(tempEval(2*size(definition)))
    tempEval = 0
    allocate(tempSurf(size(definition)))

    ! Check to make sure all entries are either a recognised operation
    ! or a surface ID. 
    ! Also checks how many elements are required to evaluate the cell.
    ! Track depth of nesting to ensure consistency in input.
    surfacePrev = .false.
    closePrev = .false.
    evalSize = 0
    sCount = 0
    depth = 0
    do i = 1, size(definition)
      
      ! Check if definition can be converted to a number
      notInt = .false.
      num = charToInt(definition(i), notInt)
      
      ! Check if definition is a recognised symbol
      if (notInt) then

        evalSize = evalSize + 1
        
        select case(definition(i))
          case('<')
            depth = depth + 1
            if (surfacePrev .or. closePrev) then
              tempEval(evalSize) = intersect
              evalSize = evalSize + 1
            end if
            tempEval(evalSize) = openBra
            closePrev = .false.
          case('>')
            depth = depth - 1
            tempEval(evalSize) = closeBra
            closePrev = .true.
          case(':')
            tempEval(evalSize) = union
            closePrev = .false.
          case('#')
            if (surfacePrev .or. closePrev) then
              tempEval(evalSize) = intersect
              evalSize = evalSize + 1
            end if
            tempEval(evalSize) = complem
            closePrev = .false.
          case default
            call fatalError(Here,'Unrecognised entry in surfaces: '//definition(i))
        end select
        surfacePrev = .false.

      ! Otherwise, it's a surface
      else

        ! Insert intersection operation in the tree
        ! as it is implicit in the input
        if (surfacePrev .or. closePrev) then
          evalSize = evalSize + 1
          tempEval(evalSize) = intersect
        end if
        evalSize = evalSize + 1
        sCount = sCount + 1
        
        tempEval(evalSize) = sCount
        surfIdx = surfs % getIdx(abs(num))
        tempSurf(sCount) % ptr => surfs % getPtr(surfIdx)
        tempSurf(sCount) % surfIdx = sign(surfIdx, num)

        surfacePrev = .true.
        closePrev = .false.

      end if

    end do

    ! Sense checks
    if (depth /= 0) call fatalError(Here,'Inconsistent use of parentheses in cell definition')
    if (tempEval(evalSize) == union) call fatalError(Here, 'Cell definition ends with a union')
    if (tempEval(evalSize) == complem) call fatalError(Here, 'Cell definition ends with a complement')

    allocate(self % surfaces(sCount))
    self % surfaces = tempSurf(1:sCount)
    deallocate(tempSurf)

    ! Apply De Morgan's law to replace complements
    call self % deMorgan(tempEval, evalSize)

    allocate(self % evalArray(evalSize))
    self % evalArray = tempEval(1:evalSize)
    deallocate(tempEval)
    
    ! Check brackets. These must be used appropriately to enforce precedence between
    ! different operations. If not, return an error due to ambiguity of operations.
    call self % checkPrecedence()

    ! Check whether the cell is simple
    hasUnion = .false.
    do i = 1, size(self % evalArray)
      if (self % evalArray(i) == union) hasUnion = .true.
    end do
    self % isSimple = .not. hasUnion

    ! Get cell ID
    call dict % get(id, 'id')
    call self % setId(id)

    ! If a simple cell, check surfaces for duplicates
    if (hasDuplicates(abs(self % surfaces % surfIdx)) .and. self % isSimple) then
      call fatalError(Here, 'There are repeated surfaces in definition of cell: '//numToChar(id))
    end if

  end subroutine init

  !!
  !! Remove complements, applying DeMorgan's law to invert surface senses
  !! and replace unions with intersections and vice versa.
  !! Applies complements to entire parentheses.
  !!
  subroutine deMorgan(self, array, sz)
    class(unionCell), intent(inout)                             :: self
    integer(shortInt), dimension(:), allocatable, intent(inout) :: array
    integer(shortInt), intent(inout)                            :: sz
    integer(shortInt)                                           :: i, j, depth, endPos
    
    i = 1
    do while (i < sz)
      if (array(i) == complem) then
        ! Delete the complement
        array(i:sz-1) = array(i+1:sz)
        sz = sz - 1

        ! If a parenthesis is opened, complement its contents.
        ! Otherwise only complement the current position
        endPos = i
        if (array(i) == openBra) then
          depth = 1
          do j = i + 1, sz
            if (array(j) == openBra) then
              depth = depth + 1
            elseif (array(j) == closeBra) then
              depth = depth - 1
            end if
            if (depth == 0) exit
          end do
          endpos = min(j, sz)
        end if

        ! Proceed until the end of the nesting, applying the complement
        do j = i, endPos
          if (array(j) < complem) then
            self % surfaces(array(j)) % surfIdx = -self % surfaces(array(j)) % surfIdx
          elseif (array(j) == union) then
            array(j) = intersect
          elseif (array(j) == intersect) then
            array(j) = union
          end if
        end do

      end if
      i = i + 1
    end do
    
  end subroutine deMorgan

  !!
  !! Check to make sure that brackets are used appropriately to avoid ambiguous operations.
  !! The use of a union on a surface followed by an intersection is ambiguous as these
  !! operations do not commute. Rather than assuming a precedence, the code will return
  !! an error if this is the case. This can be fixed by use of brackets to enforce a
  !! user-specified order of operations.
  !!
  subroutine checkPrecedence(self)
    class(unionCell), intent(in) :: self
    integer(shortInt)            :: i, j, el1, el2
    character(100), parameter    :: Here = 'checkPrecedence (unionCell_class.f90)'

    ! Go through all elements, checking for unions/intersections which
    ! are not separated by brackets.
    i = 1
    outer: do while(i < size(self % evalArray))
      
      el1 = self % evalArray(i)

      if ((el1 == union) .or. (el1 == intersect)) then

        j = i + 1

        inner: do while(j <= size(self % evalArray))

          el2 = self % evalArray(j)
          if ((el2 == intersect) .or. (el2 == union)) then

            if (el2 /= el1) then
              call fatalError(Here,'The order of operations in a cell is ambiguous')
            else
              i = j
              cycle outer
            end if

          else if((el2 == openBra) .or. (el2 == closeBra)) then
            
            ! A parenthesis is opened or closed, making the operations unambiguous.
            ! Skip to the next element after the bracket.
            i = j + 1
            cycle outer

          end if
          j = j + 1

        end do inner

      end if
      i = i + 1

    end do outer


  end subroutine checkPrecedence

  !!
  !! Return .true. if position is inside the cell
  !!
  !! Branches depending on whether the cell is simple or complex.
  !!
  !! See cell_inter for details
  !!
  pure function inside(self, r, u) result(isIt)
    class(unionCell), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: isIt
    integer(shortInt)                       :: i
    logical(defBool)                        :: halfspace, sense

    if (self % isSimple) then
    
      ! Keep compiler happy (in impossible case of cell with no surfaces)
      isIt = .false.

      do i= 1, size(self % surfaces)
        sense = self % surfaces(i) % surfIdx > 0
        halfspace = self % surfaces(i) % ptr % halfspace(r, u)

        isIt = halfspace .eqv. sense

        ! If halfspace is not equivalent to sense it means that point
        ! is outside the cell.
        if(.not.isIt) return
      end do
    else
      isIt = self % insideComplex(r,u)
    end if

  end function inside

  !!
  !! Return .true. if position is inside the cell
  !!
  !! Follows the infix notation algorithm described by Romano et al.
  !! 'Point containment algorithms for constructive solid geometry with 
  !! unbounded primitives'
  !!
  !! See cell_inter for details
  !!
  pure function insideComplex(self, r, u) result(isIt)
    class(unionCell), intent(in)            :: self
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    logical(defBool)                        :: isIt
    integer(shortInt)                       :: i, idx, depth, el
    logical(defBool)                        :: halfspace, sense

    isIt = .true.
    depth = 0
    i = 1

    do while (i <= size(self % evalArray))

      el = self % evalArray(i)
      
      ! The element is a surface
      if (el < complem) then
        idx = self % evalArray(i)
        sense = self % surfaces(idx) % surfIdx > 0
        halfspace = self % surfaces(idx) % ptr % halfspace(r, u)
        isIt = halfspace .eqv. sense

      elseif (((el == union) .and. isIt) .or. &
               ((el == intersect) .and. .not. isIt)) then

        if (depth == 0) then 
          return
        else
          depth = depth - 1
          call self % shortCircuit(i)
        end if

      elseif (el == openBra) then
        depth = depth + 1

      elseif (el == closeBra) then
        depth = depth - 1

      end if
      i = i + 1
    end do

  end function insideComplex

  !!
  !! Subroutine to short circuit inside cell evaluation.
  !! Skips the array index forward to the end of a parenthesis.
  !! Also described by Romano et al.
  !!
  pure subroutine shortCircuit(self, i)
    class(unionCell), intent(in)     :: self
    integer(shortInt), intent(inout) :: i
    integer(shortInt)                :: d, el

    d = 1
    ! Searches for end of the parenthesis
    do while (d > 0)
      i = i + 1
      el = self % evalArray(i)

      if (el == openBra) then
        d = d + 1
      else if (el == closeBra) then
        d = d - 1
      end if

    end do

  end subroutine shortCircuit

  !!
  !! Return distance to cell boundary
  !!
  !! See cell_inter for details
  !!
  pure subroutine distance(self, d, surfIdx, r, u)
    class(unionCell), intent(in)            :: self
    real(defReal), intent(out)              :: d
    integer(shortInt), intent(out)          :: surfIdx
    real(defReal), dimension(3), intent(in) :: r
    real(defReal), dimension(3), intent(in) :: u
    integer(shortInt)                       :: i
    real(defReal)                           :: test_d

    d = INF
    surfIdx = 0

    do i = 1, size(self % surfaces)
      test_d = self % surfaces(i) % ptr % distance(r, u)

      ! Select minimum distance
      if (test_d < d) then
        d = test_d
        surfIdx = i
      end if
    end do

  end subroutine distance

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(unionCell), intent(inout) :: self

    ! Call superclass procedure
    call kill_super(self)

    ! Clean local
    if(allocated(self % evalArray)) deallocate (self % evalArray)
    if(allocated(self % surfaces)) deallocate (self % surfaces)
    self % isSimple = .false.

  end subroutine kill

end module unionCell_class
