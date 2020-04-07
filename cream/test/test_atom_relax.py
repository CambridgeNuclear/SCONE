from cream.endf import ENDFTape
from cream.atom_relax import AtomRelax
from pytest import approx, raises


def test_reading_data():
    # Load ENDF Tape
    tape = ENDFTape('./files/sample.endf')
    atom = AtomRelax()
    atom.fromENDF(8, tape)

    # Load with ELement Symbol
    atom = AtomRelax()
    atom.fromENDF('O', tape)

    # Try loading non-present data
    atom = AtomRelax()
    with raises(ValueError):
        atom.fromENDF('Re', tape)


def test_repr():
    # Load ENDF Tape
    tape = ENDFTape('./files/sample.endf')
    atom = AtomRelax()
    atom.fromENDF(8, tape)
    repr(atom)


def test_checking_data():
    # Load ENDF Tape
    tape = ENDFTape('./files/sample.endf')
    atom = AtomRelax()
    atom.fromENDF(8, tape)

    # Basic check
    atom.check()

    # Change charge and fail on wrong electron sum
    atom._charge -= 1
    with raises(ValueError):
        atom.check()
    atom._charge += 1

    # Verify sumation of probabilities to 1.0 by setting unreasonable tolerance
    with raises(ValueError):
        atom.check(tol=-1.0e-12)

    # Remove one shell and verify missing shells
    atom._subshells[17] = atom._subshells.pop(3)
    with raises(ValueError):
        atom.check()
    atom._subshells[3] = atom._subshells.pop(17)

    # Verify for missmatch between number of subshells in header and dictionary
    atom._subshellNum -= 1
    with raises(ValueError):
        atom.check()
    atom._subshellNum += 1


def test_print():
    # Load ENDF Tape
    tape = ENDFTape('./files/sample.endf')
    atom = AtomRelax()
    atom.fromENDF(8, tape)
    Res1 = ('Element 8\n'
            'NSS 1\n'
            'SUBI 1\n'
            'SUBI 8\n'
            'EBI   5.372800E+02\n'
            'ELN 2\n'
            '   3   0  5.230900E+02  1.907680E-03\n'
            '   4   0  5.231300E+02  3.800270E-03\n'
            '   2   2  4.788200E+02  1.786440E-01\n'
            '   2   3  4.938600E+02  1.162240E-01\n'
            '   2   4  4.939000E+02  2.304180E-01\n'
            '   3   3  5.089000E+02  1.108220E-02\n'
            '   3   4  5.089400E+02  2.911150E-01\n'
            '   4   4  5.089800E+02  1.668090E-01\n')
    Res2 = ('Element 8\n'
            'NSS 4\n'
            ' 2 537.280000\n'
            ' 2 29.230000\n'
            ' 1 14.190000\n'
            ' 3 14.150000\n')
    assert Res1 == atom.printTrans()
    assert Res2 == atom.printGround()