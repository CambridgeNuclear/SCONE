import re

ENDF_FLOAT_REG = re.compile(r"\s*(([-+]?\d+.\d*)([-+]\d\d?))\s*")


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

    Wrapes around file handle and allows reading records from an ENDF-6
    formated file.

    """
    def __init__(self, path):
        """ Initialise ENDF file from path in filesystem

        Args:
            path (str): Path to the ENDF file

        """
        # Open file for reading
        self._handle = open(path)
        self._lineMap = dict()

        # Go through file and create line map
        # Skip first line following advice that TPID record is often
        # ill-formated
        self._handle.readline()

        # Iterate over lines saving MAT, MF, MT Locations
        temp = dict()
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
                temp[(MAT, MF, MT)] = pos
                save = False

            elif 0 in {MAT, MF, MT}:
                # Read terminato, needs to save next file location
                save = True
        # Rewind file to the top
        self._handle.seek(0)

        # Create dictionary with inventory information
        for key, pos in temp.items():
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
        pos = self._lineMap.get(MAT).get(MF).get(MT)

        if pos is not None:
            self._handle.seek(pos, 0)  # Rewind file to the location
        else:
            msg = ('Cannot set to MAT: {0} MF: {1} MT: {2} becouse this'
                   ' combination is not present in file: {3}'
                   .format(MAT, MF, MT, self._handle.name))
            raise ValueError(msg)

    def readline(self):
        """ Helper Proc for development and debug
        """
        return self._handle.readline()
