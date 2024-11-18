#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 11:51:18 2021

@author: chenghua
"""
import os
#from extract.pdbx_v2.PdbxReader import PdbxReader
from mmcif.io.IoAdapterCore import IoAdapterCore
import extract.util.cifFileCheck as cifFileCheck

## os.environ["PDB_EXTRACT"] = "/Users/chenghua/Projects/pdb_extract"

class Template():
    def __init__(self):
        """
        Create d_template to record values from template;
        Create tracking dictionary that can be used to locate the fail step

        Returns
        -------
        None.

        """
        self.d_template = {}  # ditionary for parsed template values
        self.d_track = {}  # tracking dictionary
        self.d_track["checkOK_format"] = None
        # self.d_track["openFile_OK"] = None
        self.d_track["parseFile_OK"] = None
        self.d_track["processValues_OK"] = None

    def parse(self, filepath):
        """
        process to parse cif-format template

        Parameters
        ----------
        filepath : str
            cif template file path.

        Returns
        -------
        bool
            parse processs succeeds or fails.

        """
        self.fileCheck(filepath)
        if not self.d_track["checkOK_format"]:
            return False

        if self.parseCifTemplate(filepath):
            return True
        else:
            return False

    def fileCheck(self, filepath):
        """
        Check model file's format both input and output

        Parameters
        ----------
        filepath : str
            template file path.

        Returns
        -------
        None.

        """
        cif_checker = cifFileCheck.Check(filepath)
        self.d_track["checkOK_format"] = cif_checker.checkFormat()

    def parseCifTemplate(self, filepath):
        """
        parse template, parse only the 1st data block by assumption

        Parameters
        ----------
        filepath : TYPE
            DESCRIPTION.

        Returns
        -------
        bool
            parse succceeds or fails

        """
        io = IoAdapterCore()
        try:
            l_dc = io.readFile(filepath)
            dc0 = l_dc[0]
            self.d_track["parseFile_OK"] = True
        except Exception:
            self.d_track["parseFile_OK"] = False
            return False

        try:
            for cat_name in dc0.getObjNameList():
                self.d_template[cat_name] = {}
                cat = dc0.getObj(cat_name)
                for item_name in cat.getItemNameList():
                    self.d_template[cat_name][item_name] = []
                for i in range(cat.getRowCount()):
                    d_row = cat.getRowItemDict(i)
                    for item_name in d_row:
                        self.d_template[cat_name][item_name].append(d_row[item_name])
            self.d_track["processValues_OK"] = True
        except Exception:
            self.d_track["processValues_OK"] = False
            return False

        return True


def main():
    filepath = "/Users/chenghua/Projects/pdb-extract-prod-py/tests/test_data/Templates/template_new_XRAY.cif"
    template = Template()
    template.parse(filepath)
    print(template.d_template)
    print(template.d_track)


if __name__ == "__main__":
    main()
