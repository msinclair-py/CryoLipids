#!/usr/bin/env python
from copy import deepcopy
import os
from typing import Dict, List, Union

class PDB:
    """
    This class is intended to read in, process and write out pdb files using
    the proper formatting for which they are notoriously frustrating to work
    with.
    """
    def __init__(self, 
                 filename: str, 
                 outname: str, 
                 resids: List[int] = [0],
                 old_resnames: List[str] = ['POV'], 
                 output_path: str = os.getcwd()):

        if filename[-4:] != '.pdb':
            self.filename = f'{filename}.pdb'
        elif '.' in filename:
            self.filename = f'{filename.split(".")[0]}.pdb'
        else:
            self.filename = filename

        self.outname = outname
        self.resids = resids
        self.resnames = old_resnames
        self.output_path = output_path
        self.aa = self.amino_acids
        
    @staticmethod
    def parse_pdb(pdbcontents: List[str], resid: int, resname: str) -> List[str]:
        _pdb_info = {'atom':     [0,   6],
                     'atomid':   [6,  11],
                     'atomname': [11, 17],
                     'resname':  [17, 21],
                     'chain':    [21, 22],
                     'resid':    [22, 26],
                     'x':        [26, 38],
                     'y':        [38, 46],
                     'z':        [46, 54],
                     'occ':      [54, 60],
                     'beta':     [60, 66],} 
                     #'segname':  [66, -1]}

        _r0, _r1 = _pdb_info['resid']
        _r2, _r3 = _pdb_info['resname']
        _n0, _n1 = _pdb_info['atomname']
        parsed = []
        for line in pdbcontents:
            if int(line[_r0:_r1].strip()) == resid and line[_r2:_r3].strip() == resname:
                parsed.append([line[x:y].strip() for x, y in _pdb_info.values()])
                parsed[-1][0] = 'ATOM' # explicitly convert to ATOM here
                #if 'H' not in line[_n0:_n1]:
                #    parsed.append([line[x:y].strip() for x, y in _pdb_info.values()])
                #    parsed[-1][0] = 'ATOM' # explicitly convert to ATOM here

        return parsed

    @property
    def contents(self) -> List[str]:
        """Reads in input file and slices out selected lipids by both resid and
        resname filters as defined in the config.toml file.

        Returns:
            List[str]: list of all lines pertaining to unmodeled lipids to complete.
        """
        lines = [line for line in open(self.filename, 'r').readlines() 
                 if line[:6].strip() in ['ATOM', 'HETATM']]
        
        parsed_lines = []
        for resid, resname in zip(self.resids, self.resnames):
            parsed_lines += PDB.parse_pdb(lines, resid, resname)

        return parsed_lines
    
    @staticmethod
    def format_line(line: List[Union[str, int, float]]) -> str:
        new_line = f'{line[0]:<6}{line[1]:>5} '
        
        if len(line[2]) < 4:
            new_line += f' {line[2]:<4}'
        else:
            new_line += f'{line[2]:<5}'
            
        new_line += f'{line[3]:<4}{line[4]}{line[5]:>4}    '
        new_line += f'{float(line[6]):>8.3f}{float(line[7]):>8.3f}{float(line[8]):>8.3f}'
        new_line += f'{float(line[9]):>6.2f}{float(line[10]):>6.2f}\n'

        return new_line
    
    @property
    def protein(self):
        lines = [line for line in open(self.filename, 'r').readlines() 
                 if 'ATOM' == line[:4]]
        
        prot = []
        for line in lines:
            if any(aa in line for aa in self.aa):
                prot.append(line[:60] + '  1.00' + ' ' * 12 + '\n')
        
        if self.check_multichain(prot):
            return self.rename_chains(prot)

        return prot
    
    @staticmethod
    def check_multichain(pdb: List[str]) -> int:
        cur_resid = 0
        for line in pdb:
            resid = int(line[22:26].strip())
            if resid < cur_resid:
                return 1
            else:
                cur_resid = resid

        return 0

    @staticmethod
    def rename_chains(pdb: List[str]) -> List[str]:
        chains = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        cur_resid = 0
        chain = 0
        cur_chain = chains[chain]

        renamed = []
        for line in pdb:
            resid = int(line[22:26].strip())
            if resid < cur_resid:
                chain += 1
                cur_chain = chains[chain]
            
            cur_resid = resid
            renamed.append(f'{line[:21]}{cur_chain}{line[22:]}')

        return renamed

    @property
    def amino_acids(self) -> List[str]:
        return ['ARG', 'HIS', 'HSD', 'HSE', 'HSP', 'LYS', 'ASP', 'GLU', 
                'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA', 
                'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']

    def write_to_pdb_file(self, contents: List[str]) -> str:
        """
        Write to contents 
        Args:
            contents (List[str]): _description_
            hydrogens (bool, optional): _description_. Defaults to False.
        """
        fname = os.path.join(self.output_path, 
                             f'processed_{os.path.basename(self.filename)}')
        
        with open(fname, 'w') as outfile:
            for line in contents:
                if 'H' not in line[2]:
                    outline = self.format_line(line)
                    outfile.write(outline)
                
        return fname
                
    def merge_final_pdb(self, new_coordinates: Dict[str, List[str]]) -> None:
        """
        Merges protein coordinates with new lipid coordinates. Reindexes all atoms
        to ensure numerical continuity and renumbers lipid resids to avoid duplicates.

        Args:
            new_coordinates (Dict[str, List[str]]): Per lipid coordinates to be concatenated
                                                        with protein pdb representation
        """
        out_pdb, atom_index = self.renumber(self.protein, atom_idx=1, resid_idx=None)

        residue_index = 1
        for val in new_coordinates.values():
            formatted = [self.format_line(line) for line in val]

            new_lines, atom_index = self.renumber(formatted, 
                                                  atom_idx=atom_index, 
                                                  resid_idx=residue_index)
            
            out_pdb += new_lines
            residue_index += 1
        
        with open(f'{self.output_path}/{self.outname}', 'w') as outfile:
            outfile.write(''.join(out_pdb))
    
    @staticmethod
    def renumber(lines: List[str], atom_idx: int, 
                 resid_idx: Union[None, int] = None) -> tuple():
        new_lines = []
        for line in lines:
            new_line = line[:6]
            new_line += f'{atom_idx:>5}'
            if resid_idx is None:
                new_line += line[11:]
            else:
                new_line += line[11:22] + f'{resid_idx:>4}' + line[26:]

            new_lines.append(new_line)

            atom_idx += 1
        
        return new_lines, atom_idx

def unpack_lipids(cfg: Dict[str, Dict[str, str]]) -> tuple(List[str]):
    resids, restypes, resnames = list(), list(), list()
    for val in cfg['lipids'].values():
        resids.append(val['resid'])
        restypes.append(val['restype'])
        resnames.append(val['current_resname'])
        
    return resids, restypes, resnames
