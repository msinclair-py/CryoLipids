#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2022-08-01
# Updates:
# =============================================================================
"""
check PDB format and essential content but not dictionary compliance
"""
import os
import re
import sys
import string
from extract.util.pdbFileCheck import Check

import logging
logger_name = '.'.join(["PDB_EX", __name__])
logger = logging.getLogger(logger_name)


def addChainID(filepath_in, filepath_out):
    checker = Check(filepath_in)
    checker.checkTER()  # get row indices of TER to mark the ends of polymer
    checker.checkMODEL()  # get row indices of MODEL if available for NMR models
    ## MODEL indices l_index_model will be used primarily as the start of chain ID assignment cycle (A/B/C...)
    if not checker.l_index_model:  # if no MODEL record is found (e.g. conventional X-ray structure)
        if checker.i_cryst1 >= 0:  # if CRYST1 record is present, not i_cryst1=-1 if not found
            checker.l_index_model = [checker.i_cryst1]  # assign CRYST1 index;
        elif checker.i_1st_atom >= 0:  # if no CRYST1 record, but 1st ATOM/HETATM row is found
            checker.l_index_model = [checker.i_1st_atom - 1]  # assign the row index before 1st ATOM/HETATM
        else: 
            return False  # file format must be wrong without ATOM/HETATOM

    ## get all 62 chain ID chars
    l_chars = list(string.ascii_uppercase)  # 26 upper, A/B/C...
    l_chars.extend(list(string.ascii_lowercase))  # 26 lower, a/b/c...
    for num in range(10):
        l_chars.append(str(num))  # 10 digits 0/1/2/...

    with open(filepath_in) as file_in:
        with open(filepath_out, 'w') as file_out:
            ## read file into a list of rows; row indices will be used as primary identifiers
            l_file_in = file_in.read().splitlines() 

            ## write header before 1st-MODEL, or before CRYST1 if their is no MODEL record
            for i in range(checker.l_index_model[0]+1):  
                line_in = l_file_in[i].strip()
                line_out = line_in
                file_out.write(line_out)
                file_out.write('\n')

            ## create list of indices for different sections of the coordinates
            l_index = checker.l_index_ter  # create list beased on TER as start/end of sections
            l_index.extend(checker.l_index_model)  # add row indices of MODELs 
            l_index.append(len(l_file_in)-1)  # add file ending index
            l_index.sort()  # order the indices

            k = 0 # k for chain ID index of l_chars, initiate k in case l_index_model = [-1] (No CRYST1/MODEL record present)
            for j in range(len(l_index)-1): # j for section index, start with 1st MODEL, or CRYST1 if there is no MODEL
                if l_index[j] in checker.l_index_model:
                    k = 0  #return chain ID back to 'A' for the start of each MODEL section for NMR
                try:
                    chain_id = l_chars[k]  # A/B/C../a/b/c.../0/1/2...
                except IndexError: # if exceeding 62 single char limit
                    return False  # failed adding chain ID
                
                k += 1  # go to next chain ID char for next section
                
                for i in range(l_index[j]+1, l_index[j+1]+1):  # i for row index
                    line_in = l_file_in[i].strip()
                    if len(line_in) >= 66:  # check if row is long enough to be xyz coordinates row
                        d_line = {}
                        for item in checker.d_split:  
                            l_start = checker.d_split[item][0] - 1  # l_start for start of column index within row
                            l_end = checker.d_split[item][1]  # l_end for end of column index within row
                            d_line[item] = line_in[l_start:l_end]  # keep the white space, no strip()
                        if checker.checkSplit(d_line):  # verify if the row is in xyz coordinates format
                            d_line["chainID"] = chain_id  # add chain ID
                            l_line_out = []
                            for item in checker.l_split:
                                l_line_out.append(d_line[item])
                            line_out = ''.join(l_line_out)
                        else:
                            line_out = line_in  # long non-xyz row stays the same
                    else:  
                        line_out = line_in  # short non-xyz row stays the same
                    file_out.write(line_out)
                    file_out.write('\n')
    return True
                  

def main():
    ## filepath_in = "2or2_two_chains_different_no_chainID.pdb"
    ## filepath_in = "3bn6_no_chainID_no_TER.pdb"
    ## filepath_in = "7euc_NMR_no_chainID.pdb"
    filepath_in = "7euc_NMR_no_chainID_noLigand.pdb"
    ## filepath_in = "minimal_no_chainID.pdb"
    filepath_out = "add.pdb"
    addChainID(filepath_in, filepath_out)

if __name__ == "__main__":
    main()
    
