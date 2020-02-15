program test

  use numPrecision
  use charMap_class, only : charMap

  implicit none
    character(nameLen) :: temp
    type(charMap)      :: map
    integer(shortInt)  :: i

    ! Load map
    temp = '94239.03'
    call map % add(temp, 1)

    temp = '94240.03'
    call map % add(temp, 2)

    temp = '94240.03'
    call map % add(temp, 2)

    temp = '94241.03'
    call map % add(temp, 3)

    temp = '31000.03'
    call map % add(temp, 4)

    !! Print for debug
    print *, map % map % key
    print *, map % map % status

    i = map % begin()
    do while ( i /= map % end())

      print *, i, map % atKey(i), map % atVal(i)

      i = map % next(i)
    end do

    !print *,

end program test




