# =============================================================================
# Author:  Chenghua Shao
# Date:    2023-12-20
# Updates:  
# =============================================================================
"""
Utility to find pattern in file list object. Used for log file parsing.
"""

import re
import logging

logger = logging.getLogger(__name__)

def getStartingIndex(l_file, pattern, n_occur="last"):
    """
    Find the starting index of the l_file for a particular string parttern

    Parameters
    ----------
    l_file  : list
        list of each row of file, or a simple list of anything
    pattern : str
        text pattern of the section start
    n_occur : int, optional
        take the 1st/2nd/3rd/... appearance of the pattern. The default is the last.

    Returns
    -------
    int
        starting line number (index) of the pattern.

    """
    re_pattern = re.compile(pattern)
    index_line = -1
    if n_occur == "last":
        for i in range(len(l_file)):
            line = l_file[i]
            if re_pattern.search(line):  # read through file and find last match
                index_line = i
    else:
        i_occur = 0
        for i in range(len(l_file)):
            line = l_file[i]
            if re_pattern.search(line):
                i_occur += 1
                if i_occur == n_occur:  # stop at the n_th match
                    index_line = i
                    break
    return index_line
