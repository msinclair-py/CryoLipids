#!/usr/bin/env python
import os
from typing import List, Union

class PDB:
    """
    This class is intended to read in, process and write out pdb files using
    the proper formatting for which they are notoriously frustrating to work
    with.
    """
    def __init__(self, filename: str, resids: List[int], output_path: str = os.getcwd(),
                 resname: str = 'POV'):
        self.filename = filename
        self.resids = resids
        self.output_path = output_path
        self.resname = resname


    @staticmethod
    def parse_pdb(pdbcontents):
        _pdb_info = {'atom':     [0, 6],
                     'atomid':   [6, 11],
                     'atomname': [11, 17],
                     'resname':  [17, 20],
                     'chain':    [20, 23],
                     'resid':    [23, 26],
                     'x':        [26, 38],
                     'y':        [38, 46],
                     'z':        [46, 54],
                     'occ':      [54, 60],
                     'beta':     [60, 66],
                     'segname':  [66, -1]}

        parsed = []
        for line in pdbcontents:
            parsed.append([line[x:y] for x,y in _pdb_info.values()])

        return parsed


    @property
    def pdb_file(self):
        contents = open(self.filename, 'r').readlines()
        contents = [line for line in contents if line[:4] in ['ATOM', 'HETA']]
        return PDB.parse_pdb(contents)
    
    
    @staticmethod
    def format_line(line: Union[str, int, float]) -> str:
        _formatting = ['<6', '<5', '<6', '<3', '<3', '<3',
                      '<12', '<8', '<8', '<6' '<6']
        newline = [f'{l:{f}}' for l, f in zip(line, _formatting)]
        return ''.join(newline)


    def write_to_pdb_file(self):
        with open(self.filepath, 'w') as outfile:
            for line in self.lipid:
                outline = format_line(line)
                outfile.write(outline)
