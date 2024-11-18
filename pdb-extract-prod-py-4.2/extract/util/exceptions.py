#!/usr/bin/env python3
# =============================================================================
# Author:  Chenghua Shao
# Date:    2022-08-01
# Updates:
# =============================================================================
"""
Define PDB_extract specific exceptions
"""


class ExtractError(Exception):
    pass


class InputParameterError(ExtractError):
    pass

class CreateFolderError(ExtractError):
    pass

class FileCopyError(ExtractError):
    pass


class PdbFormatError(ExtractError):
    pass


class CifFormatError(ExtractError):
    pass


class MaxitError(ExtractError):
    pass


class LogParseError(ExtractError):
    pass


class TemplateParseError(ExtractError):
    pass
