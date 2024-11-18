#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2021-12-07
# Updates:
# =============================================================================
"""
convert PDB model to mmCIF format by invoking pre-installed maxit command
"""
import os
from extract.process_model.runMaxit import Maxit
import extract.util.pdbFileCheck as pdbFileCheck

import logging
logger_name = '.'.join(["PDB_EX", __name__])
logger = logging.getLogger(logger_name)


class PdbModel():
    """
    Class to handle PDB model conversion
    """

    def convert(self, filepath_in, filepath_out, filepath_log):
        """
        convert PDB/CIF model to mmCIF format by invoking pre-installed maxit

        Parameters
        ----------
        filepath_in : str
            filepath for model input.
        filepath_out : str
            filepath for model output, requires write privilege.
        filepath_log : str
            filepath for maxit log, requires write privilege.

        Returns
        -------
        bool
            True for conversion success, with mmCIF generated at filepath_out.

        """
        # check input model file before conversion
        logger.debug("""Class PdbModel runs, filepath_in=%s; filepath_out=%s,
                     filepath_log=%s
                     """ % (filepath_in, filepath_out, filepath_log))
        pdb_checker = pdbFileCheck.Check(filepath_in)
        if not pdb_checker.checkFormat():
            logger.error("Failed PDB format check on %s" % filepath_in)
            return False

        # convert model file by invoking maxit command line
        maxit = Maxit()
        logger.info("Run Maxit pdb2cif within PdbModel Class")
        maxit_returncode = maxit.pdb2cif(filepath_in, filepath_out, filepath_log)

        # check maxit return code
        if maxit_returncode != 0:  # zero indicates successful maxit run
            logger.error("Failed maxit run on converting %s" % filepath_in)
            return False

        # check output model file existence
        if not os.path.isfile(filepath_out):
            logger.error("Failure, No maxit output file of %s" % filepath_out)
            return False

        # check outpout log file existence
        if not os.path.isfile(filepath_log):
            logger.error("Failure, No maxit log file %s" % filepath_log)
            return False

        # check maxit log file
        # if not self.fileCheckMaxitLog(filepath_log)
        #     logger.error("maxit log error in log file" % filepath_log)
        #     return False

        # pass all checkings and run above
        return True

    def fileCheckMaxitLog(self, filepath_log):
        """
        Check maxit log file

        Parameters
        ----------
        filepath_log : str
            maxit log file path.

        Returns
        -------
        None.

        """
        logger.debug("check maxit log: %s" % filepath_log)
        with open(filepath_log) as file:
            line_1st = file.readline()
            if line_1st.strip() == "Finished!":
                return True
            else:
                return False


def main():
    pdb_extract_folder = os.getenv("PDB_EXTRACT" )
    run_folder = os.path.join(pdb_extract_folder, "tests/test_data/PDB4test")
    filepath_in = os.path.join(run_folder, "in.pdb")
    filepath_out = os.path.join(run_folder, "maxit_out.cif")
    filepath_log = os.path.join(run_folder, "maxit.log")

    model = PdbModel()
    conversion_status = model.convert(filepath_in, filepath_out, filepath_log)
    print(conversion_status)


if __name__ == "__main__":
    main()
