#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2022-09-11
# Updates: 2023-02-01 Use mmCIF parser to replace pdbx_v2 parser
# =============================================================================
"""
pdb_extract_cgi is an abbreviated pdb_extract program for web server
implementation. It simply converts structure file and generates cgi_value that
is used subsequently to create html pages
"""
import os
import sys
import warnings
import subprocess
import argparse

from extract.process_model.convertPdbModel import PdbModel
from extract.process_model.validateCifModel import validateCif
from extract.util.exceptions import *
from mmcif.io.IoAdapterCore import IoAdapterCore
from extract.util.convertCatDataFormat import convertCatObjToDict
import extract.util.pdbFileCheck as pdbFileCheck
from extract.util.addChainID import addChainID
from extract.merge.mergeLogTemplateToModel import Merger

import logging
log_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(module)s - %(funcName)s:%(lineno)d - %(message)s')
f_handler = logging.FileHandler('pdb_extract_cgi.log')
f_handler.setLevel(logging.DEBUG)
f_handler.setFormatter(log_format)

c_handler = logging.StreamHandler()
c_handler.setLevel(logging.WARNING)
c_handler.setFormatter(log_format)

logger = logging.getLogger("PDB_EX")
logger.setLevel(logging.DEBUG)
logger.addHandler(f_handler)
logger.addHandler(c_handler)


# def generateErrorLog(error, filename_error_log="pdb_extract_cgi_error_log"):
#     file = open(filename_error_log, 'w')
#     file.write(error)
#     file.write("\n")
#     file.close()
#     return True


def generateErrorLog(error, filename_error_log="pdb_extract_cgi_error_log"):
    contact_email = "deposit@deposit.rcsb.org"
    session_folder = os.getcwd().split('/')[-1]
    file = open(filename_error_log, 'w')
    file.write(error)
    file.write("\n")
    file.write("""Please re-run, and if the problem remains, please contact us 
               at %s and indicate the identifier of %s
               """ % (contact_email, session_folder))
    file.write("\n")
    file.close()
    return True


def fileCopy(filepath_source, filepath_destination):
    fcopy = subprocess.run(["cp", filepath_source, filepath_destination])
    if fcopy.returncode == 0:
        if os.path.isfile(filepath_destination):
            return True
        else:
            return False
    else:
        return False


def create_parser():
    parser = argparse.ArgumentParser(description="pdb_extract_cgi is an abbreviated \
    pdb_extract program for web server implementation. It simply converts structure \
    file and generates cgi_value that is used subsequently to create html pages")

    # add mutually exclusive group for -iPDB and -iCIF
    group_input = parser.add_mutually_exclusive_group(required=True)
    group_input.add_argument('-iPDB', dest="pdb_input", metavar='PDB-input',
                             help="input filename, if in PDB format")
    group_input.add_argument('-iCIF', dest="cif_input", metavar='mmCIF-input',
                             help="input filename, if in mmCIF format")
    return parser


def parseArgs(args):
    if args.pdb_input:
        file_format = "PDB"
        filename_author_in = args.pdb_input
    elif args.cif_input:
        file_format = "CIF"
        filename_author_in = args.cif_input
    return (file_format, filename_author_in)


def processPdbModel(filepath_in, filepath_maxit_out, filepath_out):
    filepath_log = "maxit.log"
    logger.info("Check input PDB format")
    try:
        pdb_checker = pdbFileCheck.Check(filepath_in)
    except UnicodeDecodeError:
        logger.error("PDB input file is not Unicode file, likely in wrong format or with special character. STOP!!!")
        generateErrorLog("Error: PDB input file is not Unicode file, likely in wrong format or with special character.")
        return False
    if not pdb_checker.checkFormat():
        logger.error("Failed PDB format check")
        generateErrorLog("Error: uploaded file is not standard PDB format file. Please make sure the file complies with PDB format.")
        return False
    logger.info("Passed prelimnary PDB format check, chain ID not checked yet")
    if not pdb_checker.checkChainID():
        logger.warning("Missing Chain ID. Try adding chain ID")
        filepath_checked = filepath_in + "_add_chainID"
        if addChainID(filepath_in, filepath_checked):
            logger.info("Added chain ID")
            pdb_checker_2 = pdbFileCheck.Check(filepath_checked)
            if not pdb_checker_2.checkFormat():
                logger.info("Failed format check after adding chain ID")
                generateErrorLog("Error: uploaded file is not standard PDB format file, even after chain ID(s) addition has been attempted. Please add chain ID(s) yourself and standardize the file format")
                return False
        else:
            logger.info("Failed to add chain ID")
            generateErrorLog("Error: uploaded file is not standard PDB format file, and chain ID(s) cannot be added. Please add chain ID(s) yourself and standardize the file format")
            return False
    else:
        logger.info("Passed format check after adding chain ID")
        filepath_checked = filepath_in

    logger.info("Converting PDB starts")
    pdb_model = PdbModel()  # create model conversion processs instance
    if pdb_model.convert(filepath_checked, filepath_maxit_out, filepath_log):  # convert PDB to mmCIF
        logger.info("Converted PDB %s to CIF %s" % (filepath_checked, filepath_maxit_out))
        if validateCif(filepath_maxit_out):
            logger.info("Converted mmCIF is OK")
            logger.info("truncate maxit output")
            try:
                merger = Merger()
                merger.readModel(filepath_maxit_out)
                merger.truncateMaxitCat()
                merger.write(filepath_out)
                return True
            except Exception as e:
                logger.exception(e)
                generateErrorLog("Error: Failed to truncate/cleanup the file converted from PDB.")
                return False
        else:
            logger.error("Fail to validate mmCIF converted from PDB %s" % filepath_maxit_out)
            generateErrorLog("Error: Failed to validate mmCIF converted from PDB.")
            return False
    else:
        logger.error("Fail to convert PDB model file to mmCIF")
        generateErrorLog("Error: Failed to convert PDB to mmCIF.")
        return False

def processCifModel(filepath_in, filepath_out):
    try:
        if validateCif(filepath_in):
            logger.info("input mmCIF model file is OK")
            if fileCopy(filepath_in, filepath_out):
                return True
            else:
                generateErrorLog("Error: Failed to copy input cif file.")
                return False
        else:
            logger.error("fail to validate input mmCCIF model file")
            generateErrorLog("Error: uploaded file is not standard mmCIF format file. Please make sure the file complied with mmCIF format.")
            return False
    except UnicodeDecodeError:
        logger.error("CIF input file is not Unicode file, likely in wrong format or with special character. STOP!!!")
        generateErrorLog("Error: CIF input file is not Unicode file, likely in wrong format or with special character.")
        return False

def generateSummaryForCGI(filename_processed, filename_cgi_value="cgi_value"):
    try:
        # with open(filename_processed) as file:
        #     reader = PdbxReader(file)
        #     logger.debug("opening converted file at %s" % filename_processed)
        #     try:
        #         l_dc = []
        #         reader.read(l_dc)
        #         dc0 = l_dc[0]
        #         logger.info("Read converted model file")
        #     except Exception as e:
        #         logger.exception(e)
        #         generateErrorLog("Error: Failed to parse converted mmCIF.")
        io = IoAdapterCore()
        logger.debug("opening converted file at %s" % filename_processed)
        try:
            l_dc = io.readFile(filename_processed)
            dc0 = l_dc[0]
            logger.info("Read converted model file")
        except Exception as e:
            logger.exception(e)
            generateErrorLog("Error: Failed to parse converted mmCIF.")
    except IOError as e:
        logger.exception(e)
        generateErrorLog("Error: Failed to read converted mmCIF.")

    logger.info("Extract information from processed file %s" % filename_processed)
    if "symmetry" in dc0.getObjNameList():
        d_symmetry = convertCatObjToDict(dc0.getObj("symmetry"))
    else:
        d_symmetry = {}
    if "cell" in dc0.getObjNameList():
        d_cell = convertCatObjToDict(dc0.getObj("cell"))
    else:
        d_cell = {}
    if "entity_poly" in dc0.getObjNameList():
        d_entity_poly = convertCatObjToDict(dc0.getObj("entity_poly"))
    else:
        d_entity_poly = {}

    logger.info("Write to cgi_value")
    file = open(filename_cgi_value, 'w')
    if d_symmetry and d_cell:
        file.write('a=  %s :b=  %s :c= %s :alpha= %s :beta= %s :gamma= %s :space_group="%s"' %
                   (d_cell["_cell.length_a"][0], d_cell["_cell.length_b"][0], d_cell["_cell.length_c"][0],
                    d_cell["_cell.angle_alpha"][0], d_cell["_cell.angle_beta"][0], d_cell["_cell.angle_gamma"][0],
                    d_symmetry["_symmetry.space_group_name_H-M"][0]))
    else:
        file.write('a= ? :b= ? :c= ? :alpha= ? :beta= ? :gamma= ? :space_group= ?')
    file.write('\n')
    
    if d_entity_poly:
        l_entities = []
        for i in range(len(list(d_entity_poly.values())[0])):
            l_entity = []
            l_entity.append("ENTITY_ID=%s" % d_entity_poly["_entity_poly.entity_id"][i])
            l_entity.append("ENTITY_TYPE=%s" % d_entity_poly["_entity_poly.type"][i])
            l_entity.append("ONE_LETTER_SEQUNCE=%s" % ''.join(d_entity_poly["_entity_poly.pdbx_seq_one_letter_code"][i].split()))
            l_entity.append("CHAIN_ID=%s" % d_entity_poly["_entity_poly.pdbx_strand_id"][i])
            l_entity.append("TARGET_DB=%s" % d_entity_poly["_entity_poly.pdbx_target_identifier"][i])
            l_entity.append("SEQUNCE_DB= ")
            l_entity.append("SEQUNCE_NAME= ")
            l_entities.append(" : ".join(l_entity))
        file.write(" ; ".join(l_entities))
    else:
        l_entity = []
        l_entity.append("ENTITY_ID= ")
        l_entity.append("ENTITY_TYPE= ")
        l_entity.append("ONE_LETTER_SEQUNCE= ")
        l_entity.append("CHAIN_ID= ")
        l_entity.append("TARGET_DB= ")
        l_entity.append("SEQUNCE_DB= ")
        l_entity.append("SEQUNCE_NAME= ")
        file.write(" : ".join(l_entity))
    file.write('\n')
    file.close()


def main():
    parser = create_parser()
    args = parser.parse_args()
    (file_format, filepath_in) = parseArgs(args)

    logger.info("Current working directory: %s" % os.path.abspath(os.getcwd()))
    logger.debug("Input file format: %s" % file_format)
    logger.debug("Input file path: %s" % filepath_in)
    filepath_out = "pdb_extract_cgi_out.cif"
    if file_format == "PDB":
        filepath_maxit_out = "maxit_out.cif"
        if processPdbModel(filepath_in, filepath_maxit_out, filepath_out):
            generateSummaryForCGI(filepath_maxit_out)
    elif file_format == "CIF":
        if processCifModel(filepath_in, filepath_out):
            generateSummaryForCGI(filepath_out)
    else:
        logger.info("Error: Wrong file format selected %s" % file_format)
        generateErrorLog("Error: Wrong file format selected %s" % file_format)


if __name__ == "__main__":
    main()
