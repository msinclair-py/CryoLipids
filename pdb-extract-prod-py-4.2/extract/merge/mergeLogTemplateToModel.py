#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2021-12-06
# Updates:
#   2022-06-07 CS refactor
# =============================================================================
"""
Merge log data and template data into model
Priority: model > template > log
"""
import os
import json

TOP_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mmcif.io.IoAdapterCore import IoAdapterCore
from extract.util.convertCatDataFormat import convertDictToCatObj, convertCatObjToDict
from extract.util.mergeCats import mergeCatItemsUpdateValuesByRowIndex, mergeCatItemsUpdateValuesByKeys, mergeCatItemsUnionValuesByKeys
# import dict with software class
from data.em_nmr_software.em_nmr_software_class import d_em_software_lower, d_nmr_software_lower

import logging
logger_name = '.'.join(["PDB_EX", __name__])
logger = logging.getLogger(logger_name)


class Merger():
    """ Class to merge log data and template data into converted model file
    """

    def __init__(self):
        self.l_dc = None  # container list for model file
        self.dc0 = None  # container for the 1st data section of model file

    def readModel(self, filepath_model):
        """
        Read model file data into class variables

        Parameters
        ----------
        filepath_model : str
            filepath of converted model file after running maxit.

        Returns
        -------
        None.

        """
        try:
            # with open(filepath_model) as file:
            # reader = PdbxReader(file)
            self.io = IoAdapterCore()
            logger.debug("Merge starts by opening converted file at %s" %
                         filepath_model)
            try:
                # self.l_dc = []
                # reader.read(self.l_dc)
                self.l_dc = self.io.readFile(filepath_model)
                self.dc0 = self.l_dc[0]
                logger.info("Read converted model file into merger")
            except Exception as e:
                logger.exception(e)
        except IOError as e:
            logger.exception(e)

    def addCat(self, d_cat):
        """
        Add data category to model container, based on cat dictionary
        If the category object exists in model, replace it with the d_cat

        Parameters
        ----------
        d_cat : dict
            cat dictionary with keys of mmcif items for the category, and
            values as list for each key.

        Returns
        -------
        bool
            True for successful addition.

        """
        if not d_cat:
            return False
        cat_obj = convertDictToCatObj(d_cat)
        if cat_obj:
            self.dc0.append(cat_obj)  # add cat obj to data container
            cat_name = cat_obj.getName()
            logger.debug("added cat %s to model" % cat_name)
            return True
        else:
            logger.debug("cannot add cat %s to model")
            return False

    def mergeCat(self, d_cat_from, tf_replace=False):
        """
        Merge cat dictionary data into the current data category in model

        Parameters
        ----------
        d_cat_from : dict
            cat dictionary with keys of mmcif items for the category, and
            values as list for each key.

        Returns
        -------
        bool
            True for successful merge.

        """
        if not d_cat_from:
            return False

        cat_name = list(d_cat_from.keys())[0].split('.')[0].strip('_')
        logger.info("merge log/template data into category %s" % cat_name)
        if cat_name not in self.dc0.getObjNameList():
            self.addCat(d_cat_from)
            return True
        cat_obj_in_model = self.dc0.getObj(cat_name)
        d_cat_in_model = convertCatObjToDict(cat_obj_in_model)
        if d_cat_in_model:
            if cat_name == "software":
                d_cat_updated = mergeCatItemsUnionValuesByKeys(
                    d_cat_in_model, d_cat_from, ["_software.name", "_software.classification"])
            elif cat_name == "entity_poly":
                d_cat_updated = mergeCatItemsUpdateValuesByRowIndex(
                    d_cat_in_model, d_cat_from, True)
            else:
                d_cat_updated = mergeCatItemsUpdateValuesByRowIndex(
                    d_cat_in_model, d_cat_from, tf_replace)
        else:
            d_cat_updated = d_cat_from
        cat_obj_updated = convertDictToCatObj(d_cat_updated)
        if not cat_obj_updated:
            logger.debug("cannot converted merged cat %s" % cat_name)
        self.dc0.replace(cat_obj_updated)
        logger.debug("merged cat %s to model OK")
        return True

    def mergeLog(self, d_log):
        """
        Merge parsed log meta data of all process into model container.
        Choose model data over log data if in conflict

        Parameters
        ----------
        d_log : dict
            log data in 3-level dictionary.
            d_log[process]: 1st-level key is the process.
            d_log[process][category]: 2nd-level key is the category name without underscore
            d_log[process][category] is cat dictionary with keys of mmcif items
            for the category, and values as list for each key.


        Returns
        -------
        bool
            True for successful process of log merge process for at least 1 cat

        """
        if not d_log:
            logger.info("No log data parsed, no merge")
            return None

        tf_log_data_exist = False
        n_merged_cat = 0
        for process in d_log:
            d_process = d_log[process]
            if d_process:
                tf_log_data_exist = True
            else:
                logger.info("No log data parsed for process %s" % process)
                continue

            for cat_name in d_process:
                logger.info("processing log data for %s" % cat_name)
                d_cat = d_process[cat_name]
                if cat_name in self.dc0.getObjNameList():
                    if self.mergeCat(d_cat):
                        logger.info("merged log data for %s OK" % cat_name)
                        n_merged_cat += 1
                    else:
                        logger.info(
                            "failed merging log data for %s" % cat_name)
                else:
                    if self.addCat(d_cat):
                        logger.info("added log data for %s OK" % cat_name)
                        n_merged_cat += 1
                    else:
                        logger.info("failed adding log data for %s" % cat_name)
        if tf_log_data_exist:
            if not n_merged_cat:
                logger.info(
                    "no data category was merged/added from log to model")
                return False
            else:
                logger.info(
                    "merged/added %s categories from log to model" % n_merged_cat)
                return True
        else:
            return None

    def mergeTemplate(self, d_template):
        """
        Merge parsed template meta data into model container.
        Choose model data over template data if in conflict.

        Parameters
        ----------
        d_template : dict
            template data in 2-level dictionary.
            d_template[category]: key is the category name without underscore
            d_template[category] is cat dictionary with keys of mmcif items
            for the category, and values as list for each key.

        Returns
        -------
        None.

        """

        if not d_template:
            logger.info("No template data parsed, no merge")
            return None
        n_merged_cat = 0
        for cat_name in d_template:
            d_cat = d_template[cat_name]
            if cat_name in self.dc0.getObjNameList():
                if self.mergeCat(d_cat):
                    logger.info("merged template data for %s OK" % cat_name)
                    n_merged_cat += 1
                else:
                    logger.info(
                        "failed merging template data for %s" % cat_name)
            else:
                if self.addCat(d_cat):
                    logger.info("added template data for %s OK" % cat_name)
                    n_merged_cat += 1
                else:
                    logger.info(
                        "failed adding template data for %s" % cat_name)
        if not n_merged_cat:
            logger.info(
                "no data category was merged/added from template to model")
            return False
        else:
            logger.info(
                "merged/added %s categories from template to model" % n_merged_cat)
            return True

    def truncateMaxitCat(self):
        l_cat_to_keep = ["atom_site",
                         "atom_site_anisotrop",
                         "atom_sites",
                         "cell",
                         "symmetry",
                         "refine",
                         "refine_ls_restr",
                         "refine_ls_shell",
                         "pdbx_refine_tls",
                         "pdbx_refine_tls_group",
                         "refine_analyze",
                         "refine_hist",
                         "refine_ls_restr_ncs",
                         "software",
                         "em_software",
                         "pdbx_nmr_software"
                         ]
        # l_cat_to_remove = ["pdbx_database_status",
        #                    "entity_poly_seq",
        #                    "entity",
        #                    "pdbx_poly_seq_scheme",
        #                    "pdbx_nonpoly_scheme",
        #                    "chem_comp",
        #                    "struct_asym",
        #                    "pdbx_validate_close_contact",
        #                    "pdbx_validate_torsion",
        #                    "pdbx_unobs_or_zero_occ_atoms",
        #                    "atom_type",
        #                    "pdbx_entity_nonpoly"
        #                    ]
        l_cat = []  # must create a separate list to go though during deletion
        for each in self.dc0.getObjNameList():
            l_cat.append(each)

        for cat_name in l_cat:  # cannot use "for cat_name in self.dc0.getObjNameList() because of deletion will mess up iteration"
            if cat_name not in l_cat_to_keep:
                self.dc0.remove(cat_name)
        # print(self.dc0.getObjNameList())

    def processSoftwareXRAY(self, d_software_author):
        """Process X-ray software
        Merge author's software
        Merge is in the form or union, i.e. taking all softwares from model and authors

        Args:
            d_software_author (dict): software from author's commandline or webpage input

        Returns:
            bool: merge succeeds or fails
        """
        if not d_software_author:
            return False
        if "software" in self.dc0.getObjNameList():
            d_software_model = convertCatObjToDict(self.dc0.getObj("software"))
            d_software = mergeCatItemsUnionValuesByKeys(d_software_model, d_software_author, [
                "_software.name", "_software.classification"])
            self.dc0.remove("software")
            self.addCat(d_software)
        else:
            self.addCat(d_software_author)
        return True

    def processSoftwareEM(self, d_em_software_author):
        """Process EM software
        Merge author's software
        Before merge, convert/merge software category to em_software, 
        because maxit-converted mmCIF file from PDB format file has 
        software category for all methods

        Args:
            d_em_software_author (dict): software from author's commandline or webpage input

        Returns:
            bool: merge succeeds or fails
        """
        if "software" in self.dc0.getObjNameList():
            c_software = self.dc0.getObj("software")
            l_attr = c_software.getAttributeList()
            if "name" not in l_attr:  # check if name attr is present
                d_software_model_converted = {} 
            else:
                # add classification and version if not present
                if "classification" not in l_attr:
                    c_software.appendAttributeExtendRows("classification", "?")
                if "version" not in l_attr:
                    c_software.appendAttributeExtendRows("version", "?")

                d_software_model = convertCatObjToDict(c_software)
                d_software_model_converted = {"_em_software.name": [
                ], "_em_software.category": [], "_em_software.version": []}
                for i in range(len(d_software_model["_software.name"])):
                    if d_software_model["_software.classification"][i].lower() == "refinement":
                        d_software_model_converted["_em_software.name"].append(
                            d_software_model["_software.name"][i])
                        d_software_model_converted["_em_software.category"].append(
                            "MODEL REFINEMENT")
                        d_software_model_converted["_em_software.version"].append(
                            d_software_model["_software.version"][i])
                    else:
                        name = d_software_model["_software.name"][i]
                        if name.lower() in d_em_software_lower:
                            category = d_em_software_lower[name.lower()]
                        else:
                            category = '?'
                        d_software_model_converted["_em_software.name"].append(
                            name)
                        d_software_model_converted["_em_software.category"].append(
                            category)
                        d_software_model_converted["_em_software.version"].append(
                            d_software_model["_software.version"][i])

            self.dc0.remove("software")
        else:
            d_software_model_converted = {}

        if "em_software" in self.dc0.getObjNameList():
            d_em_software_model = convertCatObjToDict(
                self.dc0.getObj("em_software"))
        else:
            d_em_software_model = {}

        # merge _software and _em_software if both are in model use d_software_model_converted as primary source because:
        # if the input is PDB, then d_software_model_converted is from PDB header, while d_em_software_model is from template
        # if the input if CIF, then d_software_model_converted is supposed to be empty
        d_em_software_model_temp = mergeCatItemsUnionValuesByKeys(d_software_model_converted, d_em_software_model, [
            "_em_software.name", "_em_software.category"])
        d_em_software = mergeCatItemsUnionValuesByKeys(d_em_software_model_temp, d_em_software_author, [
                                                       "_em_software.name", "_em_software.category"])

        self.addCat(d_em_software)
        return True

    def processSoftwareNMR(self, d_nmr_software_author):
        """Process NMR software
        Merge author's software
        Before merge, convert/merge software category to pdbx_nmr_software, 
        because maxit-converted mmCIF file from PDB format file has 
        software category for all methods

        Args:
            d_nmr_software_author (dict): software from author's commandline or webpage input

        Returns:
            bool: merge succeeds or fails
        """
        if "software" in self.dc0.getObjNameList():
            c_software = self.dc0.getObj("software")
            l_attr = c_software.getAttributeList()
            if "name" not in l_attr:  # check if name attr is present
                d_software_model_converted = {} 
            else:
                # add classification and version if not present
                if "classification" not in l_attr:
                    c_software.appendAttributeExtendRows("classification", "?")
                if "version" not in l_attr:
                    c_software.appendAttributeExtendRows("version", "?")

                d_software_model = convertCatObjToDict(self.dc0.getObj("software"))
                d_software_model_converted = {"_pdbx_nmr_software.name": [
                ], "_pdbx_nmr_software.classification": [], "_pdbx_nmr_software.version": []}
                for i in range(len(d_software_model["_software.name"])):
                    if d_software_model["_software.classification"][i].lower() == "refinement":
                        d_software_model_converted["_pdbx_nmr_software.name"].append(
                            d_software_model["_software.name"][i])
                        d_software_model_converted["_pdbx_nmr_software.classification"].append(
                            "refinement")
                        d_software_model_converted["_pdbx_nmr_software.version"].append(
                            d_software_model["_software.version"][i])
                    else:
                        name = d_software_model["_software.name"][i]
                        if name.lower() in d_nmr_software_lower:
                            classification = d_nmr_software_lower[name.lower()]
                        else:
                            classification = '?'
                        d_software_model_converted["_pdbx_nmr_software.name"].append(
                            name)
                        d_software_model_converted["_pdbx_nmr_software.classification"].append(
                            classification)
                        d_software_model_converted["_pdbx_nmr_software.version"].append(
                            d_software_model["_software.version"][i])

            self.dc0.remove("software")
        else:
            d_software_model = {}
            d_software_model_converted = {}

        if "pdbx_nmr_software" in self.dc0.getObjNameList():
            d_nmr_software_model = convertCatObjToDict(
                self.dc0.getObj("pdbx_nmr_software"))
        else:
            d_nmr_software_model = {}

        # merge _software and _pdbx_nmr_software if both are in model use d_software_model_converted as primary source because:
        # if the input is PDB, then d_software_model_converted is from PDB header, while d_nmr_software_model is from template
        # if the input if CIF, then d_software_model_converted is supposed to be empty
        d_nmr_software_model = mergeCatItemsUnionValuesByKeys(d_software_model_converted, d_nmr_software_model, [
                                                              "_pdbx_nmr_software.name", "_pdbx_nmr_software.classification"])
        d_nmr_software = mergeCatItemsUnionValuesByKeys(d_nmr_software_model, d_nmr_software_author, [
                                                        "_pdbx_nmr_software.name", "_pdbx_nmr_software.classification"])
        if d_nmr_software:
            try:
                filepath_nmr_software_author = os.path.join(
                    TOP_DIR, "data/em_nmr_software/nmr_software_author.json")
                with open(filepath_nmr_software_author) as file_nmr:
                    d_software_author = json.load(file_nmr)

                d_software_author_lower = {}
                for each in d_software_author:
                    d_software_author_lower[each.lower(
                    )] = d_software_author[each]

                n_rows = len(list(d_nmr_software.values())[0])
                if not "_pdbx_nmr_software.authors" in d_nmr_software:
                    d_nmr_software["_pdbx_nmr_software.authors"] = [
                        '?'] * n_rows

                for i in range(n_rows):
                    software = d_nmr_software["_pdbx_nmr_software.name"][i]
                    if software.lower() in d_software_author_lower:
                        d_nmr_software["_pdbx_nmr_software.authors"][i] = d_software_author_lower[software.lower(
                        )]

                self.addCat(d_nmr_software)
                return True
            except Exception as e:
                logger.exception(e)
                return False
        else:
            return False

    def addRefine(self, phasing_method=''):
        """Add refine category if not exist
        Merge phasing_method from author's command line or webpage input

        Args:
            phasing_method (str, optional): phasing method

        Returns:
            bool: process succeeds or fails
        """
        d_refine = {}
        d_refine["_refine.pdbx_refine_id"] = ["X-RAY DIFFRACTION"]
        d_refine["_refine.pdbx_diffrn_id"] = ['1']
        if phasing_method:
            d_refine["_refine.pdbx_method_to_determine_struct"] = [
                phasing_method]

        if "refine" in self.dc0.getObjNameList():
            if self.mergeCat(d_refine):
                logger.info("merge refine category succeeds")
                return True
            else:
                logger.warning("merge refine category fails")
                return False
        else:
            if self.addCat(d_refine):
                logger.info("add refine category for X-ray succeeds")
                return True
            else:
                logger.warning("add refine category for X-ray fails")
                return False

    def write(self, filepath_out):
        """
        Write the current data container list into a file

        Parameters
        ----------
        filepath_out : str
            output filepath.

        Returns
        -------
        None.

        """
        # try:
        #     file = open(filepath_out, 'w')
        # except IOError as e:
        #     logger.exception(e)
        # else:
        #     try:
        #         writer = PdbxWriter(file)
        #         writer.write(self.l_dc)
        #     except Exception as e:
        #         logger.exception(e)
        # finally:
        #     file.close()
        #     logger.info("wrote merged model data to file %s" % filepath_out)
        try:
            self.io.writeFile(filepath_out, self.l_dc)
        except IOError as e:
            logger.exception(e)


# # moved test to unit test
# def main():
#     filepath_model = '/Users/chenghua/Projects/pdb-extract-prod-py/tests/test_data/PDB_XDS_Aimless_Phaser_Phenix_2021_10_8_18_16442/maxit.cif'
#     filepath_out = '/Users/chenghua/Projects/pdb-extract-prod-py/tests/test_data/PDB_XDS_Aimless_Phaser_Phenix_2021_10_8_18_16442/merged.cif'

#     d_log = {'scaling':
#              {'reflns':
#               {'_reflns.entry_id': ['UNNAMED'],
#                '_reflns.pdbx_diffrn_id': ['1'],
#                '_reflns.pdbx_ordinal': ['1'],
#                '_reflns.d_resolution_low': ['50.00'],
#                '_reflns.d_resolution_high': ['1.67'],
#                '_reflns.test': ["9999"]},
#                   }
#                  }

#     d_template = {'pdbx_database_status':
#                   {'_pdbx_database_status.entry_id': ['ToBeAssigned'],
#                    '_pdbx_database_status.dep_release_code_coordinates': ['HOLD FOR PUBLICATION'],
#                    '_pdbx_database_status.dep_release_code_sequence': ['HOLD FOR RELEASE']
#                    },
#                   '_pdbx_audit_support':
#                       {'_pdbx_audit_support.funding_organization': ['National Institutes of Health/National Institute of General Medical Sciences'],
#                        '_pdbx_audit_support.country': ['United States'],
#                        '_pdbx_audit_support.grant_number': ['XX000000']
#                        }
#                       }
#     d_software = {"_software.name": ["DENZO", "SCALEPACK"],
#                   "_software.classification":["data reduction", "data scaling"]
#                   }
#     merger = Merger()
#     merger.readModel(filepath_model)
#     merger.truncateMaxitCat()
#     merger.addCat(d_software)

#     merger.mergeLog(d_log)
#     merger.mergeTemplate(d_template)

#     merger.addRefine()
#     merger.write(filepath_out)


# if __name__ == "__main__":
#     main()
