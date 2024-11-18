#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 11:36:23 2021
pdb_extract parameter parser

@author: chenghua
"""
import os
import argparse

# updates from pdb_extract 3.x
# add -EM, -ECRYSTAL to differentiate experimenal method
# remove -iLOG, it can not be easily incorprated with argparse and it's really not necessary
# remove -d, -sp, these options only recorded software but do not offer log parsing, these are complimentary info that can be collected at DepUI 
# The new PDB_extract only handles file format conversion, log parse, and template incorporation

def create_parser():
    parser = argparse.ArgumentParser(description="PDB-extract is a program for \
                                     converting structure file format, \
                                         parsing log files, \
                                     and incorporating meta data of your \
                                         structure solution.",
                                     epilog="**Attention! DO NOT use '-iLOG' \
                                     before log file in PDB_extract 4.0 onward. \
                                         'README' file shows examples of usage.")
    # skip allow_abbrev option because it's not compatible with python 2

    # add mutually exclusive group for NMR and EM, default is XRAY
    group_method = parser.add_mutually_exclusive_group(required=False)
    group_method.add_argument('-NMR', action="store_true",
                              help="use for NMR method")
    group_method.add_argument('-EM', action="store_true",
                              help="use for Electron Microscopy method")
    group_method.add_argument('-ECRYSTAL', action="store_true",
                              help="use for Electron Crystallography method")

    # add mutually exclusive group for -iPDB and -iCIF
    group_input = parser.add_mutually_exclusive_group(required=True)
    group_input.add_argument('-iPDB', dest="pdb_input", metavar='PDB-input',
                       help="input filename, if in PDB format")
    group_input.add_argument('-iCIF', dest="cif_input", metavar='mmCIF-input',
                       help="input filename, if in mmCIF format")

    # add template and output
    parser.add_argument('-iENT', dest='template', 
                        help="optional author-provided template filename")
    parser.add_argument('-o', dest='cif_output', 
                        help="optional output filename, mmCIF format")

    # add phasing method
    parser.add_argument('-e', dest='phasing_method', 
                        help="optional phasing method")
	
    # add software and log
    parser.add_argument('-r', dest='refinement', nargs='+', 
                        help="optional software for refinement, followed by optional log")
    parser.add_argument('-s', dest='scaling', nargs='+',
                        help="optional software for reflection data scaling and merging, followed by optional log")
    parser.add_argument('-i', dest='indexing_integration', nargs='+',
                        help="optional software for reflection data indexing/integration, followed by optional log")
    parser.add_argument('-p', dest='phasing', nargs='+',
                        help="optional software for phasing, followed by optional log")
    parser.add_argument('-m', dest='molecular_replacement', nargs='+',
                        help="optional software for molecular replacement, followed by optional log")
    return parser


def parseArgs(args):
    d_arg = {}
    # process method args
    l_method = ["NMR", "EM", "ECRYSTAL"]
    for method in l_method:
        if vars(args)[method]:
            d_arg["method"] = method
            break
    else:
        d_arg["method"] = "XRAY"

    # process model args
    d_arg["model"] = {}
    if args.pdb_input:
        d_arg["model"]["format"] = "PDB"
        d_arg["model"]["filepath_author_in"] = args.pdb_input
    elif args.cif_input:
        d_arg["model"]["format"] = "CIF"
        d_arg["model"]["filepath_author_in"] = args.cif_input
    else:
        return
    if args.cif_output:
        d_arg["model"]["filepath_author_out"] = args.cif_output
    else:
        d_arg["model"]["filepath_author_out"] = "pdb_extract_out.cif"

    # process template args
    d_arg["template"] = {}
    if args.template:
        d_arg["template"]["filepath_author_in"] = args.template

    # process phasing method
    d_arg["phasing_method"] = ''
    if args.phasing_method:
        d_arg["phasing_method"] = args.phasing_method
		
    # process software args
    d_arg["software"] = {}
    l_process = ["refinement","scaling","indexing_integration","phasing","molecular_replacement"]
    for process in l_process:
        if vars(args)[process]:
            d_arg["software"][process] = {}
            # breakpoint()
            d_arg["software"][process]["name"] = vars(args)[process][0]
            try:
                d_arg["software"][process]["filepath_author_in"] = vars(args)[process][1]
            except IndexError:
                pass
    return d_arg

def main():
    # run command line test with 
    # python3 parseArgs.py -iPDB test.pdb -r CNS -s Scalepack scale.log
    parser = create_parser()
    # args = parser.parse_args()
    par = "-iPDB test.pdb -r CNS -s Scalepack scalepack.log -iENT template.cif \
        -e MR -o testout.cif"
    args = parser.parse_args(par.split())
    print(vars(args))
    d_arg = parseArgs(args)
    print()
    print(d_arg)

    par = "-ECRYSTAL -iPDB test.pdb -r CNS -s Scalepack scalepack.log \
        -iENT template.cif -o testout.cif"
    args = parser.parse_args(par.split())
    print(vars(args))
    d_arg = parseArgs(args)
    print()
    print(d_arg)


if __name__ == "__main__":
    main()
