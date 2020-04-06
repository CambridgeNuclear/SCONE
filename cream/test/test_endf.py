import cream.endf as endf
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

    #############################################
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

    # Find by MAT
    c16 = tape.get(MAT=600)
    assert (600, 1, 451) in c16
    assert (600, 23, 535) in c16

    # Find by MF
    mf23 = tape.get(MF=28)
    assert (800, 28, 533) in mf23
    assert (600, 27, 506) not in mf23

    # Find by MT
    mt2 = tape.get(MT=2)
    assert (1122, 3, 2) in mt2

    # Find non-existant
    non = tape.get(MAT=-123, MF=3, MT=13)
    assert non == []

    #############################################
    # Test HEAD, CONT, TEXT & DIR records
    # Read General Info Block
    tape.set_to(1122, 1, 451)
    data = tape.readHEAD()
    assert (1.102200E+4, 2.180550E+1, 1, 0, 0, 0,) == data
    data = tape.readCONT()
    assert (0.0, 1.0, 0, 0, 0, 6) == data
    data = tape.readCONT()
    assert (1, 2.0E+7, 1, 0, 10, 7) == data
    tape.readCONT()

    data = tape.readTEXT()
    ref = " 11-Na- 22 NEA        RCOM-JUN83 Scientific Co-ordination Group   "
    assert ref == data
    assert not tape.readSEND()
    assert not tape.readFEND()
    assert not tape.readMEND()

    # Skip to DIR Record
    for _ in range(80):
        tape.readTEXT()
    data = tape.readDIR()
    assert (3, 22, 7, 0) == data

    #############################################
    # Test HEAD and LIST records
    # Read Relaxation
    tape.set_to(800, 28, 533)
    data = tape.readHEAD()
    assert (8.0E+3, 15.8619530, 0, 0, 4, 0) == data
    data = tape.readLIST()
    assert (1.0, 0.0, 0, 0, 54, 8) == data[0:6]
    assert [5.3728E+2, 2.0] == data[6][0:2]  # 1st Line Results
    assert [1.668090E-1, 0.0, 0.0] == data[6][51:54]

    #############################################
    # TEST TAB2
    # Read Na-22 Elastic Scattering mu-distribution header
    tape.set_to(1122, 4, 2)
    data = tape.readHEAD()
    assert (1.102200E+4, 2.180550E+1, 0, 1, 0, 0) == data
    data = tape.readCONT()
    assert (0.0, 2.18055E+1, 0, 2, 0, 0) == data
    data = tape.readTAB2()
    assert (0.0, 0.0, 0, 0, 1, 57) == data[0:6]
    assert [57] == data[6]  # NBT
    assert [2] == data[7]  # INT

    #############################################
    # TEST TAB1
    # Read Photointeraction Data. Pair Production XS
    tape.set_to(600, 23, 516)
    data = tape.readHEAD()
    assert (6000.0, 1.190780E+1, 0, 0, 0, 0) == data
    data = tape.readTAB1()
    assert (0.0, 0.0, 0, 0, 2, 126) == data[0:6]
    assert [2, 126] == data[6]
    assert [2, 5] == data[7]
    assert 1.022E+6 == data[8][0]
    assert 0.0 == data[9][0]
    assert 1.0E+11 == data[8][125]
    assert 3.5893E-1 == data[9][125]

    #############################################
    # Test Setting to invalid location
    with raises(ValueError):
        tape.set_to(601, 23, 516)

    #############################################
    # Test misc
    repr(tape)
