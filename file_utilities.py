#!/usr/bin/env python
import os
from typing import List, Union

class PDB:
    """
    This class is intended to read in, process and write out pdb files using
    the proper formatting for which they are notoriously frustrating to work
    with.
    """
    def __init__(self, filename: str, resids: List[int] = [0], 
                 output_path: str = os.getcwd(), resname: str = 'POV'):

        if filename[-4:] != '.pdb':
            self.filename = f'{filename}.pdb'
        elif '.' in filename:
            self.filename = f'{filename.split(".")[0]}.pdb'
        else:
            self.filename = filename

        self.resids = resids
        self.output_path = output_path
        self.resname = resname
        


    @staticmethod
    def parse_pdb(pdbcontents: List[str], resids: List[int]) -> List[str]:
        _pdb_info = {'atom':     [0,   6],
                     'atomid':   [6,  11],
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

        _r0, _r1 = _pdb_info['resid']
        parsed = []
        for line in pdbcontents:
            if int(line[_r0:_r1].strip()) in resids or not resids[0]:
                parsed.append([line[x:y] for x,y in _pdb_info.values()])

        return parsed


    @property
    def contents(self) -> List[str]:
        lines = open(self.filename, 'r').readlines()
        lines = [line for line in lines if self.resname in line] 
        return PDB.parse_pdb(lines, self.resids)
    
    
    @staticmethod
    def format_line(line: Union[str, int, float]) -> str:
        _formatting = ['<6', '<5', '<6', '<3', '<3', '<3',
                      '<12', '<8', '<8', '<6', '<6']
        newline = [f'{l:{f}}' for l, f in zip(line, _formatting)]
        return f'{"".join(newline)}\n'


    def write_to_pdb_file(self, contents: List[str]) -> None:
        with open(f'{self.output_path}/processed_{self.filename}', 'w') as outfile:
            for line in contents:
                outline = PDB.format_line(line)
                outfile.write(outline)


class rtfParser:
    """
    This class parses CHARMM rtf parameter files to obtain lipid bonded
    parameters to be used in cryo-em lipid construction.
    """
    def __init__(self):
        pass
