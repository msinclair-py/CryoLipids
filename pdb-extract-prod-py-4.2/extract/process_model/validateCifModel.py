#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2021-12-07
# Updates:
# =============================================================================
"""
Validate mmCIF
"""
import os
import extract.util.cifFileCheck as cifFileCheck
import extract.util.cifDictCheck as cifDictCheck

def validateCif(filepath):
    """
    validate mmCIF format

    Parameters
    ----------
    filepath_in : str
        filepath for model input.

    Returns
    -------
    bool
        True for validation OK.

    """
    cif_checker = cifFileCheck.Check(filepath)
    if cif_checker.checkFormat():
        if cif_checker.checkMandatoryCat():
            return True
    return False

def main():
    # run_folder = "/Users/chenghua/Projects/pdb_extract/tests/test_data/CIF4test"
    # filepath = os.path.join(run_folder, "D_1000267334_model-upload_P1.cif.V1")
    filepath = "/Users/chenghua/Projects/pdb_extract/tests/test_temp/pdb_extract_00voip3w/maxit_out.cif"
    validate_status = validateCif(filepath)
    print(validate_status)


if __name__ == "__main__":
    main()
