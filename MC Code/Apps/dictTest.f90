program dictTest

  use numPrecision
  use dictionary_class,   only : dictionary, charLen
  use IOdictionary_class, only : IOdictionary

  implicit none

  ! Variables declaration
  type(IOdictionary) :: IOdict
  type(dictionary)   :: dict

  integer(shortInt)                            :: i,k
  integer(shortInt), dimension(:), allocatable :: i_alloc
  integer(shortInt), dimension(:), pointer     :: i_ptr => null()

  real(defReal)                            :: r
  real(defReal), dimension(:), allocatable :: r_alloc
  real(defReal), dimension(:), pointer     :: r_ptr => null()

  character(nameLen)                            :: c
  character(nameLen), dimension(:), allocatable :: c_alloc
  character(nameLen), dimension(:), pointer     :: c_p => null()

  character(nameLen), dimension(:), allocatable :: keys1

  ! Main body

  ! Read IO dictionary
  call IOdict % initFrom("./testDictInput")
  dict = IOdict



  print *, "RANK 0 GET TESTS FOR ORIGINAL AND COPY "

  call IOdict % get(i,'keyword7')
  print *, "IOdict integer keyword7 : ", i
  call dict % get(i,'keyword7')
  print *, "dict integer keyword7 : ", i

  call IOdict % get(r,'keyword7')
  print *, "IOdict real keyword7 : ", r
  call dict % get(r,'keyword7')
  print *, "dict real keyword7 : ", r

  call IOdict % get(c,'keyword3')
  print *, "IOdict char keyword3 : ", c
  call dict % get(c,'keyword3')
  print *, "dict char keyword3 : ", c

  print *, "RANK 1 GET TESTS FOR ORIGINAL AND COPY "

  ! RANK 1 INT
  call IOdict % get(i_alloc,'listint')
  call IOdict % get(i_ptr,'listint' )
  print *, "IOdict int array allocatable : ", i_alloc," and ptr ", i_ptr
  call dict % get(i_alloc,'listint')
  call dict % get(i_ptr,'listint' )
  print *, "dict int array allocatable : ", i_alloc," and ptr ", i_ptr

  ! RANK 1 REAL
  call IOdict % get(r_alloc,'list1')
  call IOdict % get(r_ptr,'list1' )
  print *, "IOdict real array allocatable : ", r_alloc," and ptr ", r_ptr
  call dict % get(r_alloc,'list1')
  call dict % get(r_ptr,'list1' )
  print *, "dict real array allocatable : ", r_alloc," and ptr ", r_ptr

  ! RANK 1 CHAR
  call IOdict % get(c_alloc,'list3')
  call IOdict % get(c_p,'list3')
  print*, "IOdict char array alloc : ", c_alloc, " and ptr ", c_p


  ! getOrDefault tests
  print *, "GET-OR-DEFAULT TESTS "
  call Iodict % getOrDefault(c,'dummy','abc')
  call Iodict % getOrDefault(c_alloc,'dummy',['abc'])
  print *, len(c_alloc), charLen
  call IOdict % getOrDefault(c_p,'dummy',['abc'])
  print *, "Character tests scalar : ", c, " Alloc ", c_alloc, " Ptr ", c_p

  call Iodict % getOrDefault(r,'dummy',1.1_defReal)
  call Iodict % getOrDefault(r_alloc,'dummy',[1.1_defReal,1.1_defReal])
  call IOdict % getOrDefault(r_ptr,'dummy',[1.1_defReal,1.1_defReal])
  print *, "Real tests scalar : ", r, " Alloc ", r_alloc, " Ptr ", r_ptr

  call Iodict % getOrDefault(i,'dummy',2)
  call Iodict % getOrDefault(i_alloc,'dummy',[2,2])
  call IOdict % getOrDefault(i_ptr,'dummy',[2,2])
  print *, "Int tests scalar : ", i, " Alloc ", i_alloc, " Ptr ", i_ptr

  ! Sub dictionary
  print *, "NESTED DICTIONERY IS COPIED"
  call IOdict % get(dict,'dictionary1')
  call dict % get(r,'keyword')
  call dict % get(i_alloc,'list2')
  print *, i_alloc, r

  call dict % keys(keys1)
  print *, keys1

  call dict % keysDict(keys1)
  print *, keys1

  call dict % keysReal(keys1)
  print *, keys1


  ! Avoid false report of memory leak from valgrind
  deallocate(c_alloc)
  deallocate(c_p)
  deallocate(i_alloc)
  deallocate(i_ptr)
  deallocate(r_alloc)
  deallocate(r_ptr)

end program dictTest



