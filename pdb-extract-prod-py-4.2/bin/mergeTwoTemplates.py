#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2022-08-01
# Updates: 2023-02-01 Use mmCIF parser to replace pdbx_v2 parser
# =============================================================================
"""
merge two templates, specially for merging author-input sequence from Server to template
"""
import sys
import logging
# from extract.pdbx_v2.PdbxReader import PdbxReader
# from extract.pdbx_v2.PdbxWriter import PdbxWriter
from mmcif.io.IoAdapterCore import IoAdapterCore

logger = logging.getLogger(__name__)

def mergeTwoTemplates(filepath_to, filepath_from, filepath_merged, l_cat_to_merge, tf_replace_if_exist=True):
    io1 = IoAdapterCore()
    try:
        l_dc_to = io1.readFile(filepath_to)
        dc0_to = l_dc_to[0]
    except Exception as e:
        logger.exception(e)

    io2 = IoAdapterCore()
    try:
        l_dc_from = io2.readFile(filepath_from)
        dc0_from = l_dc_from[0]
    except Exception as e:
        logger.exception(e)

    for cat_name in l_cat_to_merge:
        if cat_name in dc0_from.getObjNameList():
            if cat_name in dc0_to.getObjNameList():
                if tf_replace_if_exist:
                    dc0_to.append(dc0_from.getObj(cat_name))  # append function will essentially do replace
            else:
                dc0_to.append(dc0_from.getObj(cat_name))
    
    try:
        io1.writeFile(filepath_merged, l_dc_to)
        logger.info("wrote merged model data to file %s" % filepath_merged)
    except Exception as e:
        logger.exception(e)


def main():
    # filepath_to = "data_template_1"
    # filepath_from = "entity_poly.cif"
    # filepath_merged = "data_template_2"
    filepath_to = sys.argv[1]
    filepath_from = sys.argv[2]
    filepath_merged = sys.argv[3]
    l_cat_to_merge = ["entity_poly"]
    mergeTwoTemplates(filepath_to, filepath_from, filepath_merged, l_cat_to_merge)


if __name__ == "__main__":
    main()
