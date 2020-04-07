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
