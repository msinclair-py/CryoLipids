#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2022-08-01
# Updates:
# =============================================================================
"""
check file against wwPDB mmCIF dictionary
Designed for checking on mmCIF-formatted input and output of PDB_extract
Skip the parent-child relationship check
"""
import sys
import os
import re
from mmcif.io.IoAdapterCore import IoAdapterCore
from mmcif.io.IoAdapterPy import IoAdapterPy
#from extract.pdbx_v2.PdbxReader import PdbxReader
#from extract.pdbx_v2.DictionaryApi  import DictionaryApi
from mmcif.api.DictionaryApi  import DictionaryApi
from extract.util.exceptions import *
from extract.util.convertCatDataFormat import convertCatObjToDict

import logging
logger_name = '.'.join(["PDB_EX", __name__])
logger = logging.getLogger(logger_name)


class Dict:
    """ Dictionary check 
    """
    
    def __init__(self, filepath_dict):
        """
        Read dicionary into a dictionary API
        Initialize class variables for model file checking

        Returns
        -------
        None.

        """
        __lfh=sys.stderr
        __verbose=True
        __containerList = []
        
        if not os.path.isfile(filepath_dict):
            return False
        logger.info("use dictionary file at %s" % filepath_dict)
        # with open(filepath_dict) as file:
        #     reader = PdbxReader(file)
        #     reader.read(__containerList)   
        io = IoAdapterPy()
        __containerList = io.readFile(inputFilePath=filepath_dict, enforceAscii=False)
        self.dApi=DictionaryApi(containerList=__containerList,consolidate=True,verbose=__verbose,log=__lfh)
        self.catIndex = self.dApi.getCategoryIndex()
        self.l_dc = []
        self.l_cat_not_in_dict = []
        self.l_item_not_in_dict = []
        self.d_value_failRE = {}
        self.d_value_failEnum = {}
        self.d_value_failBoundary = {}
        self.l_item_with_duplicate_value = []

    def readModelCif(self, filepath):
        """
        Read model cif file into self.l_dc

        Parameters
        ----------
        filepath : str
            filepath of model cif file.

        Returns
        -------
        None.

        """
        logger.info("check %s against dictionary" % filepath)
        io = IoAdapterCore()
        self.l_dc = io.readFile(filepath)
        # with open(filepath) as file:
        #     reader = PdbxReader(file)
        #     reader.read(self.l_dc)    
    

    def valueInBoundary(self, value, l_boundary):
        """
        Check whether a value is within a defined dictionary boundary
        The input value must be either number or string that with code
        of ("int","float","non_negative_int","positive_int") in dictionary


        Parameters
        ----------
        value : str or number
            value to check.
        l_boundary : list of tuples, can be of the following scenarios:

            _refine.ls_R_factor_R_free
            [('1.0', '1.0'), ('0.0', '1.0'), ('0.0', '0.0')]
            where ('1.0', '1.0') means higher boundary can be equal;
            ('0.0', '1.0') means 0.0<x<1.0;
            ('0.0', '0.0') means lower boundary can be equal.

            _refine.ls_number_reflns_obs
            [('0', '.'), ('0', '0')]
            where ('0', '.') means higher boundary is positive infinity;
            ('0.0', '0.0') means lower boundary can be equal.

        Returns
        -------
        bool
            True if in the boundary.

        """
        try:
            value = float(value)
        except ValueError:
            return False

        l_equal = []
        l_range = []
        for (low, high) in l_boundary:
            if low == high:
                l_equal.append(low)
            else:
                l_range = [low, high]

        if not l_range:
            return True

        if l_range[0] == ".":  # lower boundary is negative infinity
            if value < float(l_range[1]):
                return True
            else:
                if l_equal:  # if value can be equal
                    if value == float(l_equal[0]):
                        return True
                return False
        elif l_range[1] == ".":  # higher boundary is positive infinity
            if value > float(l_range[0]):
                return True
            else:
                if l_equal:
                    if value == float(l_equal[0]):
                        return True
                return False
        else:
            if value > float(l_range[0]) and value < float(l_range[1]):
                return True
            else:
                if l_equal:
                    for equal_limit in l_equal:
                        if value == float(equal_limit):
                            return True
                return False

    def check(self):
        """
        Check model file against the dictionary

        Returns
        -------
        None.

        """
        for dc in self.l_dc:
            for cat in dc.getObjNameList():
                if cat in self.catIndex:
                    cat_obj = dc.getObj(cat)
                    self.checkCat(cat_obj)
                else:
                    self.l_cat_not_in_dict.append(cat)
    
    def getLowerCaseEnum(self, l_enum):
        d_enum_lower = {}
        for each in l_enum:
            d_enum_lower[each.lower()] = each
        return d_enum_lower
            
    
    def checkCat(self, cat_obj):
        """
        Check category object against the dictionary

        Parameters
        ----------
        cat_obj : class object
            data category object from cif parser.

        Returns
        -------
        None.

        """
        cat = cat_obj.getName()
        if cat == "diffrn":
            l_diffrn_id = cat_obj.getAttributeValueList("id")
            if len(l_diffrn_id) != len(set(l_diffrn_id)):
                self.l_item_with_duplicate_value.append("_diffrn.id")
        d_cat = convertCatObjToDict(cat_obj)
        for item in cat_obj.getItemNameList():
            attr = item.split('.')[1]            
            if attr in self.catIndex[cat]:
                code = self.dApi.getTypeCode(cat, attr)
                pattern = self.dApi.getTypeRegex(cat, attr)
                re_pattern = re.compile(r'%s' % pattern)
                l_enum = self.dApi.getEnumListAlt(cat, attr)
                if l_enum:
                    d_enum_lower = self.getLowerCaseEnum(l_enum)
                if code in ("int","float","non_negative_int","positive_int"):
                    l_boundary = self.dApi.getBoundaryList(cat, attr)
                else:
                    l_boundary = []
                l_values = d_cat[item]
                for value in l_values:
                    if value in ("?", "."):
                        continue
                    if re_pattern.search(value):
                        if l_enum:
                            if value.lower() not in d_enum_lower:
                                if item in self.d_value_failEnum:
                                    self.d_value_failEnum[item].append(value)
                                else:
                                    self.d_value_failEnum[item] = [value]
                        if l_boundary:
                            if not self.valueInBoundary(value, l_boundary):
                                if item in self.d_value_failBoundary:
                                    self.d_value_failBoundary[item].append(value)
                                else:
                                    self.d_value_failBoundary[item] = [value]
                    else:
                        if item in self.d_value_failRE:
                            self.d_value_failRE[item].append(value)
                        else:
                            self.d_value_failRE[item] = [value]
            else:
                self.l_item_not_in_dict.append(item)

    def reportError(self, filepath):
        """
        Write report to a log file dedicated for cif dictionary check

        Parameters
        ----------
        filepath : str
            filepath for log file.

        Returns
        -------
        None.

        """
        logger.info("report error in %s" % filepath)
        file = open(filepath, 'w')
        file.write("## List of data category names not in the dictionary:")
        file.write("\n")
        if self.l_cat_not_in_dict:
            for cat in self.l_cat_not_in_dict:
                file.write(cat)
                file.write("\n")
        else:
            file.write("None")
            file.write("\n")
        file.write("\n")

        file.write("## List of data item names not in the dictionary:")
        file.write("\n")
        if self.l_item_not_in_dict:
            for item in self.l_item_not_in_dict:
                file.write(item)
                file.write("\n")
        else:
            file.write("None")
            file.write("\n")
        file.write("\n")

        file.write("## List of values with wrong data format:")
        file.write("\n")
        if self.d_value_failRE:
            for item in self.d_value_failRE:
                file.write("Item name: %s" % item)
                file.write("\n")
                file.write("Values that are incorrect:")
                file.write("\n")
                for value in self.d_value_failRE[item]:
                    file.write(value)
                    file.write("\n")
        else:
            file.write("None")
            file.write("\n")
        file.write("\n")

        file.write("## List of values not in enumeration:")
        file.write("\n")
        if self.d_value_failEnum:
            for item in self.d_value_failEnum:
                file.write("Item name: %s" % item)
                file.write("\n")
                file.write("Values that are incorrect:")
                file.write("\n")
                for value in self.d_value_failEnum[item]:
                    file.write(value)
                    file.write("\n")
        else:
            file.write("None")
            file.write("\n")
        file.write("\n")

        file.write("## List of values not within boundary:")
        file.write("\n")
        if self.d_value_failBoundary:
            for item in self.d_value_failBoundary:
                file.write("Item name: %s" % item)
                file.write("\n")
                file.write("Values that are incorrect:")
                file.write("\n")
                for value in self.d_value_failBoundary[item]:
                    file.write(value)
                    file.write("\n")
        else:
            file.write("None")
            file.write("\n")
        file.write("\n")

        file.write("## List of key items with duplicate value:")
        file.write("\n")
        if self.l_item_with_duplicate_value:
            for item in self.l_item_with_duplicate_value:
                file.write("Item name: %s" % item)
                file.write("\n")
        else:
            file.write("None")
            file.write("\n")
        file.write("\n")

        file.close()

# # move all testing to unit test
# def main():
#     pdb_extract_folder = "/Users/chenghua/Projects/pdb-extract-prod-py"
#     folder_dict = os.path.join(pdb_extract_folder, "data/dictionary")
#     filename_dict = "mmcif_pdbx_v5_next.dic"
#     filepath_dict = os.path.join(folder_dict, filename_dict)
#     print(filepath_dict)

#     dict = Dict(filepath_dict)

#     folder = os.path.join(pdb_extract_folder, "tests/test_data/CIF4test")
#     filename = "wrong_value_type.cif"
#     filepath = os.path.join(folder, filename)
#     dict.readModelCif(filepath)
#     dict.check()
#     filepath_report = os.path.join(folder, "cifDictCheck.log")
#     dict.reportError(filepath_report)


# if __name__ == '__main__':
#     main()
