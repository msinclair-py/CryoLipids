#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2021-12-06
# Updates:
# =============================================================================
"""
Convert Cat from obj to d_cat, and from d_cat to obj
"""
from mmcif.api.DataCategory import DataCategory
# from extract.pdbx_v2.DataCategory import DataCategory

import logging
logger_name = '.'.join(["PDB_EX", __name__])
logger = logging.getLogger(logger_name)


def convertDictToCatObj(d_cat):
    """
    Convert cat dictionary to data category obj of pdbx parser

    Parameters
    ----------
    d_cat : dict
        cat dictionary with keys of mmcif items for the category, and
        values as list for each key.

    Returns
    -------
    cat_obj : class
        data category class for pdbx parser.

    """
    if not d_cat:
        logger.debug("cannot convert void cat dict to cat obj")
        return None

    try:
        cat_name = list(d_cat.keys())[0].split('.')[0].strip('_')
        cat_obj = DataCategory(cat_name)  # create a new cat obj
        l_items = list(d_cat.keys())
        n_rows = len(list(d_cat.values())[0])
        for item in l_items:
            attr = item.split('.')[1]
            cat_obj.appendAttribute(attr)  # add cat obj attr
        for i in range(n_rows):
            l_row = []
            for item in l_items:
                l_row.append(d_cat[item][i])
            cat_obj.append(l_row)  # add data, one row as a list
        # logger.debug("convert cat dict to cat obj for %s" % cat_name)
        return cat_obj
    except Exception as e:
        logger.debug("cannot convert cat dict to cat obj")
        logger.error(e)
        return None


def convertCatObjToDict(cat_obj):
    """
    Convert data category obj of pdbx parser to cat dictionary

    Parameters
    ----------
    cat_obj : class
        data category class for pdbx parser.

    Returns
    -------
    d_cat : dict
        cat dictionary with keys of mmcif items for the category, and
        values as list for each key.

    """
    if not cat_obj:
        logger.warning("cannot convert None type cat obj to cat dict")
        return {}

    d_cat = {}
    cat_name = cat_obj.getName()
    l_items = cat_obj.getItemNameList()
    for item in l_items:
        d_cat[item] = []
    for l_row in cat_obj.data:
        for i in range(len(l_items)):
            item = l_items[i]
            value = l_row[i]
            d_cat[item].append(value)
    # logger.debug("convert cat obj to cat dict for %s" % cat_name)
    return d_cat
