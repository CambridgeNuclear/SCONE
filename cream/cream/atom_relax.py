""" Atomic Relaxation Data

This module contains classes and procedures related to storing and exporting
atomic realaxation data to diffrent output formats.

Attributes:
    ATOMIC_SUBSHELL (dict): Dictionary that maps subshell index to subshell
        type (e.g. K, L3 etc.). Taken from [1]


[1] Perkins, S.T., and Cullen, D.E. ENDL type formats for the LLNL Evaluated
    Atomic Data Library, EADL, for the Evaluated Electron Data Library, EEDL,
    and for the Evaluated Photon Data Library, EPDL. United States: N. p.,
    1994. Web. doi:10.2172/10172308.
"""
from cream.endf import ENDFTape

ATOMIC_SUBSHELL = {1: 'K', 2: 'L', 3: 'L1', 4: 'L23', 5: 'L2', 6: 'L3', 7: 'M',
                   8: 'M1', 9: 'M23', 10: 'M2', 11: 'M3', 12: 'M45', 13: 'M4',
                   14: 'M5', 15: 'N', 16: 'N1', 17: 'N23', 18: 'N2', 19: 'N3',
                   20: 'N45', 21: 'N4', 22: 'N5', 23: 'N67', 24: 'N6',
                   25: 'N7', 26: 'O', 27: 'O1', 28: 'O23', 29: 'O2', 30: 'O3',
                   31: 'O45', 32: 'O4', 33: 'O5', 34: 'O67', 35: 'O6',
                   36: 'O7', 37: 'O89', 38: 'O8', 39: 'O9', 40: 'P', 41: 'P1',
                   42: 'P23', 43: 'P2', 44: 'P3', 45: 'P45', 46: 'P4',
                   47: 'P5', 48: 'P67', 49: 'P6', 50: 'P7', 51: 'P89',
                   52: 'P8', 53: 'P9', 54: 'P1011', 55: 'P10', 56: 'P11',
                   57: 'Q', 58: 'Q1', 59: 'Q23', 60: 'Q2', 61: 'Q3'}

SUBSHELL_IDX = dict([(name, idx) for idx, name in ATOMIC_SUBSHELL.items()])

EL_Z = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
        'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
        'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
        'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
        'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
        'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
        'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
        'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
        'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
        'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
        'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
        'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
        'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,
        'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99,
        'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105,
        'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111,
        'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117,
        'Og': 118}


class AtomRelax:
    """ Stores Atomic Relaxation Data

    Contains atomic relaxation data related to a single atom.
    Atomic relaxation data consists of:
        * Number of subshells
        * Binding energy for a subshell
        * Number of electrons in a subshell in ground state
        * List vacancy transisitions with associated probability and energy
            change for each subshell.
    """

    def __init__(self):
        self._element = None
        self._charge = None
        self._source = None
        self._subshellNum = None
        self._subshells = dict()

    def __repr__(self):
        string = 'Atomic Relaxation Data\n'
        string += '\tFrom: {}\n'.format(self._source)
        string += '{: <16}=\t{}\n'.format('\tElement', self._element)
        string += '{: <16}=\t{}\n'.format('\tZ', self._charge)
        string += '{: <16}=\t{}\n'.format('\t# Subshells', self._charge)
        return string

    def fromENDF(self, elem, tape):
        """ Read Atomic Relaxation Data from ENDF

        Procedure uses variable names that originate from ENDF specification
        in order to improve readibility note:
            EBI: Subshell binding energy [eV]
            ELN: Number of electron in subshell for ground state
            SUB_T: Target subshell for any transition
            SUB_T2: Subshell of electron emission for non-radiative transistion
            ETR: Energy of transition [eV]
            FTR: Fractional probability of transition

        Args:
            elem (int,str): Element Z number (1, 13) or symbol ('H', 'U')
            tape (ENDFTape): ENDF tape to search for data

        Raises:
            ValueError: required data is not present on ENDFtape
        """
        # Find element atomic number
        if isinstance(elem, str):
            Z = EL_Z.get(elem)
        else:
            Z = int(elem)
        MAT = Z*100

        # Check if data is available
        if not tape.has(MAT=MAT, MF=28, MT=533):
            msg = ('Could not find atomic relaxation data (MF=28, MT=533) for '
                   'the element with Z={1} on ENDFTape from '
                   'file {}'.format(Z, tape.name))
            raise ValueError(msg)

        # Read data
        # Header
        tape.set_to(MAT=MAT, MF=28, MT=533)
        self._charge, _, _, _, self._subshellNum, _ = tape.readHEAD()

        # Subshell data
        for _ in range(self._subshellNum):
            ss, _, _, _, length, _, data = tape.readLIST()

            # Change rank of the list so it is list of lists of length 6
            # if the ENDF is formatted correctly length of the list WILL BE
            # a multiple of 6
            it = iter(data)
            data = [[next(it) for _ in range(6)] for _ in range(int(length/6))]

            # Read the special first line with binding energy
            # and electron number
            EBI, ELN, _, _, _, _ = data.pop(0)
            trans = list()

            # Read all transitions for subshell ss
            for line in data:
                SUB_T = int(line[0])
                SUB_T2 = int(line[1])
                ETR = line[2]
                FTR = line[3]
                trans.append((SUB_T, SUB_T2, ETR, FTR))
            self._subshells[ss] = (EBI, ELN, trans)

        # For correctly formatted ENDF should be now at SEND record
        if not tape.readSEND():
            msg = ('ENDF tape: {} may be ill formated. Reading of atomic '
                   'relaxation data for MAT: {} did not end with SEND record. '
                   'Check that number of subshells matches number of '
                   'provided lists'.format(tape.name, MAT))
            raise ValueError(msg)
