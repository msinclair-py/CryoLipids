#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2021-12-07
# Updates:
#   2022-06-07 CS refactor
# =============================================================================
"""
pdb_extract run controller
"""
import os
import sys
import warnings
import tempfile
import subprocess
import importlib
import re
import shutil

TOP_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, TOP_DIR)

from bin.parseArgs import create_parser, parseArgs
from extract.process_model.convertPdbModel import PdbModel
from extract.process_model.validateCifModel import validateCif
from extract.extract_template.extractTemplate import Template
from extract.merge.mergeLogTemplateToModel import Merger
from extract.util.cifDictCheck import Dict
from extract.util.exceptions import *

import logging
log_format = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(module)s - %(funcName)s:%(lineno)d - %(message)s')
f_handler = logging.FileHandler('pdb_extract.log')
f_handler.setLevel(logging.DEBUG)
f_handler.setFormatter(log_format)

c_handler = logging.StreamHandler()
c_handler.setLevel(logging.WARNING)
c_handler.setFormatter(log_format)

logger = logging.getLogger("PDB_EX")
logger.setLevel(logging.DEBUG)
logger.addHandler(f_handler)
logger.addHandler(c_handler)


class Process():
    """class for managing pdb_extract run
    class dictionary d_manager records all information including:
        locations of file input, output, and temporary ones in run folder
        run status
        file check status
    """

    def __init__(self):
        """
        create temporary run folder, all input files are copied to run folder
        under static name for subsequent process. Final results are copied 
        from run folder to user's working directory

        Parameters
        ----------
        folder_run : str, optional
            author-designated run folder. The default is "".
        tf_cleanup : bool, optional
            whether to cleanup run folder after process. The default is False.

        Returns
        -------
        None.

        """
        self.d_manager = {}  # one dict to rule them all
        self.d_manager["method"] = ''  # record method
        self.d_manager["folder"] = {}  # record folder structure
        self.d_manager["model"] = {}  # record model format and filepaths
        self.d_manager["software"] = {}  # record softwares and logs by process
        self.d_manager["template"] = {}  # record template filepaths
        self.d_manager["status"] = {}  # record run status of each step
        self.d_manager["log"] = {}  # record individual run log filepath
        self.d_software = {}  # record software categoory for Xray and EC
        self.d_em_software = {}  # record em_software category for EM
        self.d_nmr_software = {}  # record pdbx_nmr_software category for NMR

    def __fileCopy(self, filepath_source, filepath_destination):
        """
        file copy from A to B

        Parameters
        ----------
        filepath_source : str
            filepath from.
        filepath_destination : str
            filepath to.

        Returns
        -------
        bool
            whether the file copy process was successful.

        """
        fcopy = subprocess.run(["cp", filepath_source, filepath_destination])
        if fcopy.returncode == 0:
            if os.path.isfile(filepath_destination):
                return True
            else:
                return False
        else:
            return False

    def parseArgsInput(self, args_text_input=""):
        """
        parse user's arguments in running pdb_extract command

        Parameters
        ----------
        args_text_input : str, optional
            arguments after pdb_extract command. The default is "", and
            arguments will be taken from stdin.

        Returns
        -------
        None.

        """
        parser = create_parser()
        if args_text_input:
            args = parser.parse_args(args_text_input.split())
        else:
            args = parser.parse_args()
        d_arg = parseArgs(args)
        # print(d_arg)
        if d_arg:
            for key in d_arg:
                self.d_manager[key] = d_arg[key]  # args to start d_manager
            self.d_manager["status"]["parseArgsInput_OK"] = True
        else:
            self.d_manager["status"]["parseArgsInput_OK"] = False
        logger.debug("d_manager after parsing input args: %s" % self.d_manager)

    def createTempFolder(self, folder_run=""):
        """
        Create temp folder to store all files during pdb_extract run

        Parameters
        ----------
        folder_run : str, optional
            run folder, under which temp folder is created. The default is cur.
        tf_cleanup : bool , optional
            whether to remove the temp folder after run. The default is False.

        Returns
        -------
        None.

        """
        self.d_manager["folder"]["current"] = os.getcwd()

        if not folder_run:
            folder_run = os.getcwd()  # by default create temp in current dir
        try:
            self.d_manager["folder"]["temp"] = tempfile.mkdtemp(
                prefix="pdb_extract_", dir=folder_run)
            self.d_manager["status"]["createTempFolder_OK"] = True
        except OSError:
            self.d_manager["status"]["createTempFolder_OK"] = False

    def copyToTemp(self):
        """
        copy all input files to run folder:
            model file, mandatory
            template file, optional
            log files, optional

        Returns
        -------
        None.

        """
        tf_copyModelOK = True
        tf_copyTemplateOK = True
        tf_copyLogOK = True

        # copy model file
        if self.d_manager["model"]["filepath_author_in"]:
            if self.d_manager["model"]["format"] == "PDB":
                filepath_model_in_temp = os.path.join(
                    self.d_manager["folder"]["temp"], "in.pdb")
            elif self.d_manager["model"]["format"] == "CIF":
                filepath_model_in_temp = os.path.join(
                    self.d_manager["folder"]["temp"], "in.cif")

            if self.__fileCopy(self.d_manager["model"]["filepath_author_in"], filepath_model_in_temp):
                self.d_manager["model"]["filepath_temp_in"] = filepath_model_in_temp
            else:
                tf_copyModelOK = False

        # copy template file
        if self.d_manager["template"]:
            if "filepath_author_in" in self.d_manager["template"]:
                if self.__fileCopy(self.d_manager["template"]["filepath_author_in"],
                                   os.path.join(self.d_manager["folder"]["temp"], "template")):
                    self.d_manager["template"]["filepath_temp_in"] = \
                        os.path.join(
                            self.d_manager["folder"]["temp"], "template")
                else:
                    tf_copyTemplateOK = False

        # copy log files
        if self.d_manager["software"]:
            for process in self.d_manager["software"]:
                if "filepath_author_in" in self.d_manager["software"][process]:
                    filepath_log_in_temp = os.path.join(
                        self.d_manager["folder"]["temp"], process+".log")
                    if self.__fileCopy(self.d_manager["software"][process]["filepath_author_in"], filepath_log_in_temp):
                        self.d_manager["software"][process]["filepath_temp_in"] = filepath_log_in_temp
                    else:
                        tf_copyLogOK = False

        if tf_copyModelOK and tf_copyTemplateOK and tf_copyLogOK:
            self.d_manager["status"]["copyToTemp_OK"] = True
        else:
            self.d_manager["status"]["copyToTemp_OK"] = False

    def processModel(self):
        """
        convert PDB model by either pdb2cif in run folder.
        generate and check maxit output cif file.
        examine CIF model of either input or output

        Returns
        -------
        None.

        """
        if self.d_manager["model"]["format"] == "PDB":
            logger.info("input file in PDB format, run maxit to convert")
            filepath_in = self.d_manager["model"]["filepath_temp_in"]
            logger.debug("maxit input file: %s" % filepath_in)
            filepath_maxit_out = os.path.join(
                self.d_manager["folder"]["temp"], "maxit_out.cif")
            filepath_log = os.path.join(
                self.d_manager["folder"]["temp"], "maxit.log")

            pdb_model = PdbModel()  # create model conversion processs instance
            # convert PDB to mmCIF
            try:
                if pdb_model.convert(filepath_in, filepath_maxit_out, filepath_log):
                    logger.debug("maxit output file: %s" % filepath_maxit_out)
                    if validateCif(filepath_maxit_out):
                        self.d_manager["model"]["filepath_temp_process"] = filepath_maxit_out
                        self.d_manager["status"]["processModel_OK"] = True
                        self.d_manager["log"]["filepath_temp_maxit_log"] = filepath_log
                    else:
                        raise CifFormatError(
                            "Wrong format for converted mmCIF file")
                        self.d_manager["status"]["processModel_OK"] = False
                else:
                    self.d_manager["status"]["processModel_OK"] = False
            except UnicodeDecodeError:
                logger.error("PDB input file is not Unicode file, likely in wrong format or with special character. STOP!!!")
                sys.exit()
        elif self.d_manager["model"]["format"] == "CIF":
            logger.info("input file in CIF format, no conversion")
            filepath_in = self.d_manager["model"]["filepath_temp_in"]
            logger.debug("cif input file: %s" % filepath_in)
            try:
                if validateCif(filepath_in):
                    self.d_manager["model"]["filepath_temp_process"] = filepath_in
                    self.d_manager["status"]["processModel_OK"] = True
                else:
                    raise CifFormatError("Wrong format for input mmCIF file")
                    self.d_manager["status"]["processModel_OK"] = False
            except UnicodeDecodeError:
                logger.error("CIF input file is not Unicode file, likely in wrong format or with special character. STOP!!!")
                sys.exit()
        else:
            pass

    def parseSoftwareLog(self):
        """
        parse all log files by importing individual parser based on 
        method, process, and software.
        generate a dictionary for all log data

        Returns
        -------
        None.

        """
        if self.d_manager["method"] != "XRAY":
            self.d_software = {}
            self.d_log = {}
            logger.warning(
                "log parsing for method %s currently not supported" % self.d_manager["method"])
            return

        self.d_software = {"_software.classification": [],
                           "_software.name": []}  # record software in dict
        self.d_log = {}  # dictionary to store all parsed log data
        # record software name and log file location for each process
        for process in self.d_manager["software"]:
            software_name = self.d_manager["software"][process]["name"]
            logger.debug("process software %s in %s" %
                         (software_name, process))
            self.d_manager["software"][process]["parseOK"] = None
            # prepare software dict to add in the next merge process
            self.d_software["_software.name"].append(software_name)
            if process == "refinement":
                self.d_software["_software.classification"].append(
                    "refinement")
            elif process == "scaling":
                self.d_software["_software.classification"].append(
                    "data scaling")
            elif process == "indexing_integration":
                self.d_software["_software.classification"].append(
                    "data reduction")
            elif process in ("phasing", "molecular_replacement"):
                self.d_software["_software.classification"].append("phasing")

            # record log file location in the temp folder if any
            if "filepath_temp_in" not in self.d_manager["software"][process]:
                self.d_log[process] = {}
                logger.info("no software log to parse for %s in process %s" %
                            (software_name, process))

        # parse log files, 4.0 only handles one log file for scaling or
        # indexing_integration. If both log files exist, parse scaling log only
        try:
            log2parse = self.d_manager["software"]["scaling"]["filepath_temp_in"]
            logger.debug("process scaling log: %s" % log2parse)
            software2parse = self.d_manager["software"]["scaling"]["name"]
            process2parse = "scaling"
        except KeyError as e:
            logger.info("no log file for scaling")
            try:
                log2parse = self.d_manager["software"]["indexing_integration"]["filepath_temp_in"]
                logger.debug("process indexing log: %s" % log2parse)
                software2parse = self.d_manager["software"]["indexing_integration"]["name"]
                process2parse = "indexing_integration"
            except KeyError as e:
                logger.info("no log file for indexing_integration")
                log2parse = None
                software2parse = None
                process2parse = None
        if process2parse:
            # remove python-filename-incompatible chars in software name to
            # match the software's parser filename in extract/extract_log
            # folder, so that the corresponding parser can be imported
            # e.g. d*TREK, cctbx.xfel
            software_name_clean = re.sub('[-/ *\.]', '', software2parse)

            # locate folder of individual log parser
            l_import_folder = ["extract",
                               "extract_log",
                               self.d_manager["method"],
                               "scaling",
                               software_name_clean.lower()]
            import_path = '.'.join(l_import_folder)
            print(import_path)
            try:
                # import log parser from __init__.py of the specific folder
                # the __init__.py chooses version specific parser
                log_parser = importlib.import_module(import_path)
                logger.debug("software %s log parser import from %s" %
                             (software_name, import_path))
            except ModuleNotFoundError:
                logger.warning("not support parsing log for software %s in process %s" %
                               (software_name, process2parse))
                self.d_manager["software"][process2parse]["parseOK"] = False
            else:
                try:
                    self.d_log[process2parse] = log_parser.run(log2parse)
                    self.d_manager["software"][process2parse]["parseOK"] = True
                except Exception as e:
                    self.d_log[process2parse] = {}
                    self.d_manager["software"][process2parse]["parseOK"] = False
        logger.debug("created new _software category: %s" % self.d_software)
        logger.debug("stored parsed log data in self.d_log with keys: %s" %
                     self.d_log.keys())

    def parseTemplate(self):
        """
        parse the template file if provided in run folder.
        parsing new mmCIF template only.
        generate a dictionary data structure for all template  data

        Returns
        -------
        None.

        """
        self.d_template = {}  # dictionary to store parsed template data
        try:
            filepath_template = self.d_manager["template"]["filepath_temp_in"]
            logger.debug("process author's template: %s" % filepath_template)
        except KeyError:
            logger.info("no author-provided template")
            self.d_manager["status"]["parseTemplate_OK"] = None
        else:
            template = Template()  # create template parser instance
            if template.parse(filepath_template):
                self.d_manager["status"]["parseTemplate_OK"] = True
                self.d_template = template.d_template
            else:
                self.d_manager["status"]["parseTemplate_OK"] = False
        logger.debug("stored parsed template data self.d_template with keys:%s"
                     % self.d_template.keys())

    def mergeLogTemplate(self):
        """
        merge both log data and template data into maxit output cif in run folder.
        resolve conflict amongst all three sources, with for following priority:
            maxit output cif > log data > template data
            by assuming it's easier to make mistake in provide wrong log and
            more aasier to have typo in template.
        generate final combined model output cif file and run content check

        Returns
        -------
        None.

        """
        filepath_model = self.d_manager["model"]['filepath_temp_process']

        merger = Merger()

        #  read mmCIF file of converted model from temp folder
        try:
            merger.readModel(filepath_model)
        except Exception as e:
            logger.exception(e)

        #  truncate non-useful categories from maxit output of PDB2CIF
        #  do not touch if author's input is mmCIF file already
        if self.d_manager["model"]["format"] == "PDB":
            logger.info("truncate maxit output")
            merger.truncateMaxitCat()

        #  merge Log data from self.d_log to model
        try:
            self.d_manager["status"]["mergeLog_OK"] = merger.mergeLog(
                self.d_log)
        except Exception as e:
            self.d_manager["status"]["mergeLog_OK"] = False
            logger.exception(e)

        #  merge template from self.d_template to mosdel
        try:
            self.d_manager["status"]["mergeTemplate_OK"] = merger.mergeTemplate(
                self.d_template)
        except Exception as e:
            self.d_manager["status"]["mergeTemplate_OK"] = False
            logger.exception(e)

        #  process software
        #  If model file has no software category, add software from author's commandline arg or webpage input
        #  If model file has software category, for X-ray entries, merge author's software
        #  For EM entries, convert software to em_software and merge
        #  For NMR entries, convert software to pdbx_nmr_software and merge
        logger.info(
            "process software and merge software info from author command line")
        if self.d_manager["method"] == "XRAY":
            merger.processSoftwareXRAY(self.d_software)
        elif self.d_manager["method"] == "EM":
            merger.processSoftwareEM(self.d_em_software)
        elif self.d_manager["method"] == "NMR":
            merger.processSoftwareNMR(self.d_nmr_software)

        #  For X-ray, merge phasing method from author's commandline or webpage input into refine category
        #  If refine cat doesn't exist, add an empty refine category that is mandatory for X-ray OneDep DepUI
        if self.d_manager["method"] == "XRAY":
            if self.d_manager["phasing_method"]:
                phasing_method = self.d_manager["phasing_method"]
            elif "molecular_replacement" in self.d_manager["software"]:
                phasing_method = "MOLECULAR REPLACEMENT"
            else:
                phasing_method = ''
            merger.addRefine(phasing_method)

        # write merged file to temp folder
        filename_out = "out.cif"
        filepath_out = os.path.join(
            self.d_manager["folder"]["temp"], filename_out)
        logger.debug("write merged file to: %s" % filepath_out)
        try:
            merger.write(filepath_out)
            self.d_manager["model"]["filepath_temp_out"] = filepath_out
        except Exception as e:
            self.d_manager["model"]["filepath_temp_out"] = None
            logger.exception(e)

    def checkAgainstDict(self):
        """
        check final mmcif output in temp against dictionry

        Returns
        -------
        None.

        """
        try:
            pdb_extract_folder = TOP_DIR
            folder_dict = os.path.join(pdb_extract_folder, "data/dictionary")
            filename_dict = "mmcif_pdbx_v5_next.dic"
            filepath_dict = os.path.join(folder_dict, filename_dict)

            dict = Dict(filepath_dict)

            filepath = self.d_manager["model"]["filepath_temp_out"]
            dict.readModelCif(filepath)
            logger.debug("check output mmcif %s against dictionary %s" %
                         (filepath, filepath_dict))
            dict.check()  # check against mmCIF dictionary
            filepath_report = os.path.join(
                self.d_manager["folder"]["temp"], "dictionary_violation.log")
            dict.reportError(filepath_report)
            self.d_manager["log"]["filepath_temp_dictionary_check"] = filepath_report
            self.d_manager["status"]["checkAgainstDict_OK"] = True
        except Exception as e:
            self.d_manager["status"]["checkAgainstDict_OK"] = False
            logger.exception(e)

    def copyResultFromTemp(self):
        """
        copy final model file and necessary logging file for PDB_extract run
        to user's working directory

        Returns
        -------
        None
        """
        try:
            if self.__fileCopy(self.d_manager["model"]["filepath_temp_out"],
                               self.d_manager["model"]["filepath_author_out"]):
                self.d_manager["status"]["copyResultFromTemp_OK"] = True
            else:
                self.d_manager["status"]["copyResultFromTemp_OK"] = False
        except Exception as e:
            self.d_manager["status"]["copyResultFromTemp_OK"] = False
            logger.exception(e)

        if "filepath_temp_dictionary_check" in self.d_manager["log"]:
            filepath_author_dictionary_check = os.path.join(
                self.d_manager["folder"]["current"], "dictionary_violation.log")
            if self.__fileCopy(self.d_manager["log"]["filepath_temp_dictionary_check"], filepath_author_dictionary_check):
                self.d_manager["log"]["filepath_author_dictionary_check"] = filepath_author_dictionary_check

        if "filepath_temp_maxit_log" in self.d_manager["log"]:
            filepath_author_maxit_log = os.path.join(
                self.d_manager["folder"]["current"], "maxit.log")
            if self.__fileCopy(self.d_manager["log"]["filepath_temp_maxit_log"], filepath_author_maxit_log):
                self.d_manager["log"]["filepath_author_maxit_log"] = filepath_author_maxit_log


def reviewDict(d_):
    l_lines = []
    l_lines.append("######Review operational dictionary######")
    l_1 = ["method"]
    l_2 = ["status", "folder", "model", "template"]
    l_3 = ["software"]
    for key1 in l_1:
        l_lines.append("%s:\t%s" % (key1, d_[key1]))
        l_lines.append("##########")
    for key1 in l_2:
        l_lines.append("%s:" % key1)
        for key2 in d_[key1]:
            l_lines.append("\t%s:\t%s" % (key2, d_[key1][key2]))
        l_lines.append("##########")
    for key1 in l_3:
        l_lines.append("%s:" % key1)
        for key2 in d_[key1]:
            l_lines.append("\t%s:" % key2)
            for key3 in d_[key1][key2]:
                l_lines.append("\t\t%s:\t%s" % (key3, d_[key1][key2][key3]))
        l_lines.append("##########")
    return l_lines


def generateErrorLog(error, filename_error_log="pdb_extract_error_log"):
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


def runPdbExtract(args="", tf_cleanup=True):
    logger.info("Current working directory: %s" % os.path.abspath(os.getcwd()))
    p = Process()  # p.d_manager is the overall tracking dictionary
    logger.info("initiate pdb_extract processs")

    # Step 1: Parse input arguments
    logger.info("Parse input arguments starts")
    logger.info("Author's arguments %s" % args)
    logger.info("temp folder clean up set as %s" % tf_cleanup)
    p.parseArgsInput(args)
    # p.parseArgsInput()  # user input
    if not p.d_manager["status"]["parseArgsInput_OK"]:
        generateErrorLog("Fail to parse input arguments.")
        raise InputParameterError("Fail Step 1: Parse input arguments, stop")
    logger.info("Parse input arguments finishes")

    # Step 2: Create temp run folder to store all files in the process
    # p.createTempFolder()  # create temp folder in current folder
    logger.info("Create temp run folder starts")
    p.createTempFolder("/tmp")  # create temp folder in system /tmp
    if not p.d_manager["status"]["createTempFolder_OK"]:
        generateErrorLog("Fail to create temp run folder.")
        raise CreateFolderError(
            "Fail Step 2: create temporary run folder to store all files in the process, stop")
    logger.info("Create temp run folder finishes")

    # Step 3: Copy input files to temp folder
    logger.info("Copy input files to temp folder starts")
    p.copyToTemp()
    if not p.d_manager["status"]["copyToTemp_OK"]:
        generateErrorLog("Fail to copy input files to temp folder.")
        raise FileCopyError(
            "Fail Step 3: Copy input files to temp folder, stop")
    logger.info("Copy input files to temp folder finishes")

    # Step 4: Process input model file: PDB2mmCIF by Maxit; or validate mmCIF
    logger.info("Process input model file starts")
    p.processModel()  # Model.d_track tracks status of conversion substeps
    if not p.d_manager["status"]["processModel_OK"]:
        generateErrorLog("Fail to process input model file.")
        raise Exception(
            "Fail Step 4: Convert/Process input PDB/CIF model file, stop")
    logger.info("Process input model file finishes")

    # Step 5: Process software names and logs
    logger.info("Process software names and logs starts")
    p.parseSoftwareLog()
    for process in p.d_manager["software"]:
        software_name = p.d_manager["software"][process]["name"]
        if p.d_manager["software"][process]["parseOK"] == False:
            warnings.warn("Fail Step 5: Parse log file for process %s software %s, continue" %
                          (process, software_name))
    logger.info("Process software names and logs finishes")

    # Step 6: Process author's template
    logger.info("Process author's template starts")
    p.parseTemplate()
    if p.d_manager["status"]["parseTemplate_OK"] == False:
        warnings.warn("Fail Step 6: Process user's template, continue")
    logger.info("Process author's template finishes")

    # Step 7: Merge log and template data into model file, and implement other needed modifications
    logger.info("Merge log and template data into model file starts")
    p.mergeLogTemplate()
    if p.d_manager["status"]["mergeLog_OK"] == False:
        warnings.warn("Fail Step 7A: Merge log data into model file, continue")
    if p.d_manager["status"]["mergeTemplate_OK"] == False:
        warnings.warn(
            "Fail Step 7B: Merge template data into model file, continue")
    logger.info("Merge log and template data into model file finishes")

    # Step 8: Run comprehensive mmCIF check against mmCIF dictionary
    logger.info("Run comprehensive mmCIF check against mmCIF dictionary starts")
    p.checkAgainstDict()
    if not p.d_manager["status"]["checkAgainstDict_OK"]:
        warnings.warn(
            "Fail Step 8: Check merged file against mmCIF dictionary, continue")
    logger.info(
        "Run comprehensive mmCIF check against mmCIF dictionary finishes")

    # Step 9: Copy result files from temp folder to user's default folder
    logger.info("Copy result files from temp folder to user's folder starts")
    p.copyResultFromTemp()
    if not p.d_manager["status"]["copyResultFromTemp_OK"]:
        generateErrorLog("Fail to copy result file from temp folder.")
        raise FileCopyError(
            "Fail Step 9: Copy result files from temp folder to user's default folder, stop")
    logger.info("Copy result files from temp folder to user's folder finishes")

    # Review tracking dictionary
    logger.warning('\n'.join(reviewDict(p.d_manager)))

    if tf_cleanup:
        if os.path.isdir(p.d_manager["folder"]["temp"]):
            shutil.rmtree(p.d_manager["folder"]["temp"])
            logger.info("Remove temp folder finishes")


def main():
    runPdbExtract()


if __name__ == "__main__":
    main()
