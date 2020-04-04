import re

ENDF_FLOAT_REG = re.compile(r"\s*(([-+]?\d+.\d*)([-+]\d\d?))\s*")

ENDF_NLIB = {0: 'ENDF/B',
             1: 'ENDF/A',
             2: 'JEFF',
             3: 'EFF',
             5: 'CENDL',
             6: 'JENDL',
             21: 'SG-23',
             31: 'INDL/V',
             32: 'INDL/A',
             33: 'FENDL',
             34: 'IRDF',
             35: 'BROND-IAEA',
             36: 'INGDB-90',
             37: 'FENDL/A',
             41: 'BROND'}

ENDF_IPART = {0: 'Photon',
              1: 'Neutron',
              11: 'Electron'}

ENDF_SUBLIB = {0: 'Photo-Nuclear Data',
               1: 'Photo-Induced Fission Product Yields',
               3: 'Photo-Atomic Interaction Data',
               4: 'Radioactive Decay Data',
               5: 'Spontaneous Fission Product Yields',
               6: 'Atomic Relaxation Data',
               10: 'Incident Neutron Data',
               11: 'Neutron-Induced Fission Product Yields',
               12: 'Thermal Neutron Scattering Data',
               113: 'Electro-Atomic Interaction Data',
               10010: 'Incident-Proton Data',
               10011: 'Proton-Induced Fission Product Yields',
               10020: 'Incident-Deuteron Data'}


ENDF_MF = {1: 'General information',
           2: 'Resonance parameter data',
           3: 'Reaction cross sections',
           4: 'Angular distributions for emitted particles',
           5: 'Energy distributions for emitted particles',
           6: 'Energy-angle distributions for emited particles',
           7: 'Thermal neutron scattering law data',
           8: 'Radioactivity and fission-product yield data',
           9: 'Multiplicities for radioactive nuclide production',
           10: 'Cross sections for radioactive nuclide production',
           12: 'Multiplicities for photon production',
           13: 'Cross sections for photon production',
           14: 'Angular distributions for photon production',
           15: 'Energy distributions for photon production',
           23: 'Photo- or electro-atomic interaction cross sections',
           26: 'Electro-atomic angle and energy distribution',
           27: ('Atomic form factors or scattering functions for photo-atomi'
                ' interactions'),
           28: 'Atomic relaxation data',
           30: 'Data covariances obtained from covariances and sensitivities',
           31: 'Data covariances for nu(bar)',
           32: 'Data covariances for resonance parameters',
           33: 'Data covariances for reaction cross-sections',
           34: 'Data covariances for angular distributions',
           35: 'Data covariances for energy distributions',
           39: 'Data covariances for radionuclide production yields',
           40: 'Data covariances for radionuclide production cross-sections'}


def endf2float(st):
    """ Convert string from ENDF file into Float

    Args:
        st (str): Character string that contains ENDF float number
            e.g. "-1.23456+6; 1.23456+38"

    Result:
       float: Converted number from the string

    Raises:
        ValueError: If conversion fails
    """
    # Match Regular Expression
    res = ENDF_FLOAT_REG.match(st)
    if res:
        return float(res[2] + 'E' + res[3])
    else:
        return float(st)


def endf2int(st):
    """ Convert sring from ENDF finto Integer

    Args:
        st (str): CHaracter string that contains ENDF integer
            e.g. "     3"; "    "; "1"
    Result:
        int: Converted integer. If string contains only blanks gives 0

    Raises:
        ValueError: If conversion fails
    """
    if st.isspace():
        return 0
    else:
        return int(st)


class ENDFTape:
    """ Class wrapper around an ENDF file on disk

    Contains a single ENDF Tape (ASCII File od disk).
    It can hold multiple evaluations. Each with MAT ID e.g 9238 for U-238
    It is assumed that combinations of MAT ID, MF file number and MT section
    numer are UNIQUE. There is no check if this is the case at the moment.

    For more information about ENDF-6 format. In particular documentation
    of diffrent data records (HEAD, TAB1, ...) refer to:

    A. Trkov, M. Herman, and D. A. Brown, ‘ENDF-6 formats manual’,
    Data Formats and Procedures for the Evaluated Nuclear Data Files ENDF/B-VI
    and ENDF/B-VII, National Nuclear Data Center Brookhaven National
    Laboratory, Upton, NY, pp. 11973–5000, 2012.

    Example:
        tape = ENDFTape(path)
        if (tape.has(MAT=600, MF=2, MT=13)):
            # Read data
            data = tape.readHEAD()
            # Do smth with data
            data = tape.readTAB1()
            ...
        else:
            raise ValueError('MAT 600, MF = 2 and MT=13 was not found')

    """
    def __init__(self, path):
        """ Initialise ENDF file from path in filesystem

        Args:
            path (str): Path to the ENDF file

        """
        # Open file for reading
        self._handle = open(path)
        self._map = dict()
        self._lineMap = dict()

        # Go through file and create line map
        # Skip first line following advice that TPID record is often
        # ill-formated
        self._handle.readline()

        # Iterate over lines saving MAT, MF, MT Locations
        save = True
        while True:
            # Save position and read a line
            pos = self._handle.tell()
            line = self._handle.readline()

            # Get parameters from line
            MAT = endf2int(line[66:70])
            MF = endf2int(line[70:72])
            MT = endf2int(line[72:75])

            # Update state
            if MAT == -1:
                # Terminate if end is reached. Must be checked first
                break

            elif save and 0 not in {MAT, MF, MT}:
                # Need to save and line is not terminator
                self._map[(MAT, MF, MT)] = pos
                save = False

            elif 0 in {MAT, MF, MT}:
                # Read terminato, needs to save next file location
                save = True
        # Rewind file to the top
        self._handle.seek(0)

        # Create dictionary with inventory information
        for key, pos in self._map.items():
            # Restructure self._lineMap into tree of dictionaries with
            # a horrible one-liner
            self._lineMap.setdefault(key[0], dict()) \
                .setdefault(key[1], dict())[key[2]] = pos

    def __del__(self):
        self._handle.close()

    def __repr__(self):
        st = ""
        st += "Opened File: " + self._handle.name
        return st

    def has(self, MAT, MF=None, MT=None):
        """ Enquire about contents of a tape

        MF and MT set to None (default) are ignored. For example:
        tape.has(600)   # Returns True if there is MAT=600 in the tape
        tape.has(600,1) # Returns True if there if MAT=600 that has MF=1

        Args:
            MAT (int): MAT=Z*100+(A mod 100) Identifier e.g
                612 for C-12; 9238 fro U-238; 100 for natural or atomic H;
            MF (int, optional): File Identifier e.g. 2 for Resonance Data
            MT (int, optional): Section ID in specific MF file e.g. 1
                in MF=3 for total XS.

        Result:
            True if the combination is present. False otherwise.

        Raises:
            ValueError: If MT is provided without MF

        """
        if MT is not None and MF is None:
            raise ValueError('MT was given without MF. You must define MF.')

        res = self._lineMap.get(MAT)
        if MF is not None and res is not None:
            res = res.get(MF)
        if MT is not None and res is not None:
            res = res.get(MT)
        return res is not None

    def get(self, MAT=None, MF=None, MT=None):
        """ Return all combinations that match MAT, MF, MT pattern

        Args:
             MAT (int, optional): MAT=Z*100+(A mod 100) Identifier e.g
                612 for C-12; 9238 fro U-238; 100 for natural or atomic H;
            MF (int, optional): File Identifier e.g. 2 for Resonance Data
            MT (int, optional): Section ID in specific MF file e.g. 1
                in MF=3 for total XS.
        Result:
            Set of (MAT, MF, MT) tuples that match the pattern.
            Empty set if no matches were found.
        """
        res = list()
        for trip in self._map:
            save = True
            save = save and (trip[0] == MAT if MAT is not None else True)
            save = save and (trip[1] == MF if MF is not None else True)
            save = save and (trip[2] == MT if MT is not None else True)
            if save:
                res.append(trip)
        return res

    def set_to(self, MAT, MF, MT):
        """ Set Position in the file to beginning of MAT-MF-MT

        Puts the file reader at the beginning of a section of an ENDF tape
        which is specified by MAT material identifier (e.g. 612), MF file
        ID (e.g. 1-general Info), and MT ID of file section (e.g. 451 for
        evaluation information)

        Args:
            MAT (int): MAT=Z*100+(A mod 100) Identifier e.g
                612 for C-12; 9238 fro U-238; 100 for natural or atomic H;
            MF (int): File Identifier e.g. 2 for Resonance Data
            MT (int): Section ID in specific MF file e.g. 1 in MF=3 for
                total XS.

        Raises:
            ValueError: If combination of MAT, MF and MT in not present in the
                tape.
        """
        if not self.has(MAT, MF, MT):
            msg = ('Cannot set to MAT: {0} MF: {1} MT: {2} becouse this'
                   ' combination is not present in file: {3}'
                   .format(MAT, MF, MT, self._handle.name))
            raise ValueError(msg)

        pos = self._lineMap.get(MAT).get(MF).get(MT)
        self._handle.seek(pos, 0)  # Rewind file to the location

    def readSEND(self):
        """ Read SEND record.

        Result:
            True if next line is SEND record.
        """
        line = self._handle.readline()
        MAT = endf2int(line[66:70])
        MF = endf2int(line[70:72])
        MT = endf2int(line[72:75])

        return (MAT != 0) and (MF != 0) and (MT == 0)

    def readFEND(self):
        """ Read FEND record

        Result:
            True if next line is FEND record.
        """
        line = self._handle.readline()
        MAT = endf2int(line[66:70])
        MF = endf2int(line[70:72])
        MT = endf2int(line[72:75])

        return (MAT != 0) and (MF == 0) and (MT == 0)

    def readMEND(self):
        """ Read MEND record

        Result:
            True if next line is MEND record.
        """
        line = self._handle.readline()
        MAT = endf2int(line[66:70])
        MF = endf2int(line[70:72])
        MT = endf2int(line[72:75])

        return (MAT == 0) and (MF == 0) and (MT == 0)

    def readTEXT(self):
        """ Read TEXT record

        Result:
            String with text characters
        """
        line = self._handle.readline()
        return line[0:66]

    def readCONT(self):
        """ Read CONT Record

        Result:
            (tuple): tuple containing:
                C1,C2 (float):
                L1, L2, N1, N2 (float):
        """
        line = self._handle.readline()
        return (endf2float(line[0:11]),
                endf2float(line[11:22]),
                endf2int(line[22:33]),
                endf2int(line[33:44]),
                endf2int(line[44:55]),
                endf2int(line[55:66]))

    def readHEAD(self):
        """ Read HEAD Record

        Has the same form as CONT

        Result:
            (tuple): tuple containig:
                ZA (float): ZZAAA ID of the nuclide e.g. 92238 for U-238
                AWR (float): Atomic Mass
                L1, L2, N1, N2 (int):
        """
        return self.readCONT()

    def readDIR(self):
        """ Read DIR Record

        Result:
            (tuple): typle containing:
                L1, L2, N1, N2 (int)
        """
        line = self._handle.readline()
        return (endf2int(line[22:33]),
                endf2int(line[33:44]),
                endf2int(line[44:55]),
                endf2int(line[55:66]))

    def readLIST(self):
        """ Read LIST Record

        Result:
            (tuple): tuple containing:
                C1,C2 (float):
                L1,L2 (int):
                NPL (int): Number of elements in the list
                N2 (int):
                list (list): List of floats
        """
        # Read first line
        head = self.readCONT()
        NPL = head[4]
        lis = list()
        for i in range((NPL - 1)//6 + 1):
            line = self._handle.readline()
            n = min(6, NPL - 6*i)
            lis += list(map(endf2float, line[0:11*n].split()))
        return (*head, lis)

    def readTAB1(self):
        """ Read TAB1 Record

        Result:
            (tuple): tuple containing:
                C1, C2 (float):
                L1, L2 (int):
                NR (int): Number of diffrent interpolation regions
                NP (int): Number of points in the table
                NBT (list): List of integers. Interpolation region boundaries
                INT (list): List of integers. Interpolation flags
                x (list): List of floats. X-values in the table
                y (list): List of floats. Y-values in the table
        """
        data = self.readCONT()
        NR = data[4]
        NP = data[5]

        # Read interpolation regions
        NBT = list()
        INT = list()
        if NR != 0:
            for i in range((NR - 1)//3 + 1):
                line = self._handle.readline()
                n = min(3, NR - 3*i)
                nums = list(map(endf2int, line[0:22*n].split()))
                NBT += nums[::2]
                INT += nums[1::2]

        # Read the table
        x = list()
        y = list()
        for i in range((NP - 1)//3 + 1):
            line = self._handle.readline()
            n = min(3, NP - 3*i)
            nums = list(map(endf2float, line[0:22*n].split()))
            x += nums[::2]
            y += nums[1::2]
        return (*data, NBT, INT, x, y)

    def readTAB2(self):
        """ Read TAB2 record

        Result:
            (tuple): tuple containing:
                C1,C2 (float):
                L1,L2 (int):
                NR (int): Number of interpolation regions in z-axis
                NZ (int): Number of points on z-axis
                NBT (list): List of integers. Location of interpolation regions
                INT (list): List of integers. Interpolation flags
        """
        data = self.readCONT()
        NR = data[4]

        # Read interpolation regions
        NBT = list()
        INT = list()
        if NR != 0:
            for i in range((NR - 1)//3 + 1):
                line = self._handle.readline()
                n = min(3, NR - 3*i)
                nums = list(map(endf2int, line[0:22*n].split()))
                NBT += nums[::2]
                INT += nums[1::2]
        return (*data, NBT, INT)
