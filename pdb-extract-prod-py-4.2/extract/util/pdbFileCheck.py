#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2022-08-01
# Updates:
# =============================================================================
"""
check PDB format and essential content but not dictionary compliance
"""
import re
import logging
logger_name = '.'.join(["PDB_EX", __name__])
logger = logging.getLogger(logger_name)

class Check():
    """ Class to check PDB format
    """    
    def __init__(self, filepath):
        """
        Attempt to read PDB file

        Parameters
        ----------
        filepath : TYPE
            DESCRIPTION.
        Returns
        -------
        None.

        """
        self.l_split = ["head", "serial", "spacer1", "atomName", "altLoc",
                        "resName", "spacer2", "chainID", "resSeq", "iCode",
                        "spacer3", "x", "y", "z", "occupancy", "tempFactor",
                        "spacer4", "element", "charge"]
        self.d_split = {}
        self.d_split["head"] = (1, 6)
        self.d_split["serial"] = (7, 11)
        self.d_split["spacer1"] = (12, 12)
        self.d_split["atomName"] = (13, 16)
        self.d_split["altLoc"] = (17, 17)
        self.d_split["resName"] = (18, 20)
        self.d_split["spacer2"] = (21, 21)
        self.d_split["chainID"] = (22, 22)
        self.d_split["resSeq"] = (23, 26)
        self.d_split["iCode"] = (27, 27)
        self.d_split["spacer3"] = (28, 30)
        self.d_split["x"] = (31, 38)
        self.d_split["y"] = (39, 46)
        self.d_split["z"] = (47, 54)
        self.d_split["occupancy"] = (55, 60)
        self.d_split["tempFactor"] = (61, 66)
        self.d_split["spacer4"] = (67, 76)
        self.d_split["element"] = (77, 78)
        self.d_split["charge"] = (79, 80)

        self.l_lines = []
        try:
            with open(filepath) as file:
                self.l_lines = file.read().splitlines()
        except IOError as e:
            logger.exception(e)
        self.i_cryst1 = None
        self.checkCryst1()  # get row index of CRYST1
        self.i_1st_atom = None
        self.checkAtom()  # get row index of the 1st ATOM/HETATM row
        self.l_index_ter = []  # get list of row indices of TER
        self.l_index_model = []  # get list of row indices of MODEL

    def checkCryst1(self):
        """
        Check the presence and format of CRYST1

        Returns
        -------
        bool
            True if all essential categories are present

        """
        re_cryst1 = re.compile(r"^CRYST1\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+[A-Z]\s")
        for i in range(len(self.l_lines)):
            line = self.l_lines[i]
            if re_cryst1.search(line):
                self.i_cryst1 = i
                logger.debug("input PDB has CRYST1 record in line %s of" % i)
                logger.debug(line)
                return True
        self.i_cryst1 = -1
        logger.debug("input PDB file has no CRYST1 record")
        return False

    def checkAtom(self):
        """
        Check the presence and format of ATOM
        based on PDB format guide 3.3
        For simplicity, check only one row of ATOM/HETATM presence
        and assume the rest are fine. 

        Returns
        -------
        bool
        True if all essential categories are present

        """
        for i in range(self.i_cryst1+1, len(self.l_lines)):
            line = self.l_lines[i].strip()
            if len(line) >= 66:
                d_line = {}
                for item in self.d_split:
                    j_start = self.d_split[item][0] - 1
                    j_end = self.d_split[item][1]
                    try:
                        d_line[item] = line[j_start:j_end].strip()
                    except IndexError:
                        d_line[item] = ''
                if self.checkSplit(d_line):  # whether is ATOM row
                    self.i_1st_atom = i
                    logger.debug("input PDB has ATOM/HETATM in line %s of" % i)
                    logger.debug(line)
                    return True
        self.i_1st_atom = -1
        logger.debug("input PDB has no ATOM/HETATM record")
        return False

    def checkSplit(self, d_line):  # doesn't check chain ID presence
        if d_line["head"].strip() not in ("ATOM", "HETATM"):
            return False

        try:
            int(d_line["serial"].strip())
        except ValueError:
            return False

        re_atomName = re.compile(r"\w+")
        if not re_atomName.search(d_line["atomName"].strip()):
            return False

        re_resName = re.compile(r"\w+")
        if not re_resName.search(d_line["resName"].strip()):
            return False

        try:
            int(d_line["resSeq"].strip())
        except ValueError:
            return False

        for item in ["x", "y", "z", "occupancy", "tempFactor"]:
            try:
                float(d_line[item].strip())
            except ValueError:
                return False
        return True

    def checkTER(self):
        """
        Check the presence and format of CRYST1

        Returns
        -------
        bool
            True if all essential categories are present

        """
        re_ter = re.compile(r"^TER\s*")
        for i in range(self.i_cryst1+1, len(self.l_lines)):
            line = self.l_lines[i]
            if re_ter.search(line):
                self.l_index_ter.append(i)
        logger.debug("input PDB has TER in line %s" % self.l_index_ter)

    def checkMODEL(self):
        """
        Check the presence and format of CRYST1

        Returns
        -------
        bool
            True if all essential categories are present

        """
        re_model = re.compile(r"^MODEL\s+\d")
        for i in range(self.i_cryst1+1, len(self.l_lines)):
            line = self.l_lines[i]
            if re_model.search(line):
                self.l_index_model.append(i)
        logger.debug("input PDB has MODEL in line %s" % self.l_index_model)

    def checkChainID(self):
        for i in range(self.i_cryst1+1, len(self.l_lines)):
            line = self.l_lines[i].strip()
            if len(line) >= 66:
                d_line = {}
                for item in self.d_split:
                    j_start = self.d_split[item][0] - 1
                    j_end = self.d_split[item][1]
                    try:
                        d_line[item] = line[j_start:j_end].strip()
                    except IndexError:
                        d_line[item] = ''
                if self.checkSplit(d_line):
                    if d_line["chainID"].strip(): ## OK if find at least one row with chain ID, because solvent/ligand does not need chain ID.
                        logger.debug("input PDB has chain ID in line %s" % i)
                        return True
        logger.debug("input PDB has no chain ID in any line")
        return False

    def checkFormat(self):
        """
        Check whether the file is proper PDB format
        The only one to checked is ATOM/HETATM presence
        Not checking TER card because of NMR entries.

        Returns
        -------
        bool
            True if in PDB format.

        """
        # if method in ("XRAY", "ED"):
        #     if self.i_cryst1 >=0:  # check presence of CRYST1
        #         if self.i_1st_atom >= 0:  # check presence of 1st good ATOM/HETATM row
        #             return True
        #     return False
        # else:
        #     if self.i_1st_atom >= 0:  # check presence of 1st good ATOM/HETATM row
        #         return True
        #     return False
        if self.i_1st_atom >= 0:
            return True
        else:
            return False

def main():
    ## filepath = sys.argv[1]
    filepath = "minimal.pdb"
    filepath = "missing_ATOM.pdb"
    print(filepath)
    checker = Check(filepath)
    print(checker.i_cryst1)
    print(checker.i_1st_atom)
    print(checker.checkFormat())
    print(checker.checkFormat())
    ## print(checker.checkChainID())
    ## print(checker.checkTER())
    ## print(checker.l_index_ter)


if __name__ == "__main__":
    main()
