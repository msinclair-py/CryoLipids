#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File:    convertModel.py
# Author:  Chenghua Shao
# Date:    2021-12-07
# Version: 1.0

# Updates:
#   2022-06-06 CS refactor, use direct args input, add tracking dictionary
"""
run pre-installed maxit command tool
"""
__docformat__ = "restructuredtext en"
__author__ = "Chenghua Shao"
__email__ = "chenghua.shao@rcsb.org"
__license__ = "Apache 2.0"

import os
import subprocess



class Maxit:
    """ Class to run maxit
    """
    def __init__(self):
        """
        Set up maxit path, command, and run env

        Returns
        -------
        None.

        """
        TOP_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        MAXIT_DIR = os.path.join(TOP_DIR, "packages", "maxit-v11.100-prod-src")
        os.environ["RCSBROOT"] = MAXIT_DIR
        self.maxit_command = os.path.join(MAXIT_DIR, "bin", "maxit")

    def pdb2cif(self, filepath_input, filepath_output, filepath_log):
        """
        Run PDB to CIF conversion

        Parameters
        ----------
        filepath_input : str
            PDB model input filepath.
        filepath_output : str
            mmCIF model output filepath.
        filepath_log : str
            maxit log filepath.

        Returns
        -------
        integer
            return status code, only zero indicates successful maxit run.

        """
        maxit_run = subprocess.run([os.path.join(self.maxit_command),
                                    "-input", filepath_input,
                                    "-output", filepath_output,
                                    "-o", "1",
                                    "-log", filepath_log])
        return maxit_run.returncode

    def cif2cif(self, filepath_input, filepath_output, filepath_log):
        """
        Run CIF to CIF conversion

        Parameters
        ----------
        filepath_input : str
            CIF model input filepath.
        filepath_output : str
            mmCIF model output filepath.
        filepath_log : str
            maxit log filepath.

        Returns
        -------
        integer
            return status code, only zero indicates successful maxit run.

        """
        maxit_run = subprocess.run([os.path.join(self.maxit_command),
                                    "-input", filepath_input,
                                    "-output", filepath_output,
                                    "-o", "8",
                                    "-log", filepath_log])
        return maxit_run.returncode


def main():
    maxit = Maxit()
    pdb_extract_folder = os.getenv("PDB_EXTRACT" )
    #folder_test = "tests/test_data/PDB_HKL2000_Phaser_Refmac_2021_9_28_13_26368"
    #folder_test = "tests/test_data/PDB_Buster_2021_10_14_11_20818/"
    #folder_test = "tests/test_data/PDB_Buster_2021_10_14_11_20818/"
    #folder_test = "tests/test_data/PDB_Xscale_Phaser_Buster_2021_10_15_16_9523/"
    #folder_test = "tests/test_data/PDB_CNSsf_Denzo_Scalepack_CNS_2021_10_19_13_26023/"
    #folder_test = "tests/test_data/PDB_XDS_Aimless_Phaser_Phenix_2021_10_8_18_16442/"
    #folder_test = "tests/test_data/PDB_Template_Scalepack_1"
    folder_test = os.path.join(pdb_extract_folder, "tests/test_data/PDB_Template_Scalepack_1")
    filepath_input = os.path.join(folder_test, "deposited_XRAY.pdb")
    filepath_output = os.path.join(folder_test, "maxit.cif")
    filepath_log = os.path.join(folder_test, "maxit.log")
    maxit_returncode = maxit.pdb2cif(filepath_input, filepath_output, filepath_log)
    print(maxit_returncode)


if __name__ == "__main__":
    main()
