from cream.endf import ENDFTape
import cream.atom_relax as ar
from pytest import approx, raises


def test_reading_data():
    # Load ENDF Tape
    tape = ENDFTape('./files/sample.endf')
    atom = ar.AtomRelax()
    atom.fromENDF(8, tape)
