import cream.ENDF as endf
from pytest import approx, raises

def test_endf_float_conversions():

    # Test ordinary float conversion
    assert 1.4567 == endf.endf2float('1.4567')
    assert 13.13 == endf.endf2float(' 1.313E+001')

    # ENDF-style float conversion
    assert -10.0 == endf.endf2float('-1.000+1')
    assert 1.23456E+38 == endf.endf2float('  +1.23456+38    ')

    # Test invalid cases
    with raises(ValueError):
        endf.endf2float('This is Not A Number')


def test_endf_int_conversions():
    # Test oridnary int conversions
    assert 2 == endf.endf2int('  +2  ')
    assert -3 == endf.endf2int('-3')

    # Test ENDF-style float conversion
    assert 0 == endf.endf2int(' ')

    # Test invalid cases
    with raises(ValueError):
        endf.endf2int('')
    with raises(ValueError):
        endf.endf2int('This is not int')


def test_endfTape():
    """ TODO: Investigate separating into a number of tests using fixture.
    Be carefull about the state of endfTape that can change
    """
    # Load File
    tape = endf.ENDFTape('./files/sample.endf')

    #######################################
    # Check content checking
    assert tape.has(MAT=1122)
    assert not tape.has(MAT=9238)

    assert tape.has(MAT=1122, MF=2)
    assert not tape.has(MAT=600, MF=2)   # MF is wrong
    assert not tape.has(MAT=6001, MF=3)  # MAT is wrong

    assert tape.has(MAT=600, MF=27, MT=505)
    assert not tape.has(MAT=600, MF=27, MT=1)    # Fail on MT
    assert not tape.has(MAT=1122, MF=27, MT=51)  # Fail on MF
    assert not tape.has(MAT=601, MF=27, MT=505)  # Fail on MAT

    # Check error
    with raises(ValueError):
        tape.has(MAT=600, MT=3)


