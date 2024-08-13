#!/usr/bin/env python
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
                if 'H' not in line[_n0:_n1]:
                    parsed.append([line[x:y].strip() for x, y in _pdb_info.values()])
                    parsed[-1][0] = 'ATOM' # explicitly convert to ATOM here

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
        return f'{"".join(new_line)}\n'
    
    @property
    def protein(self):
        lines = [line for line in open(self.filename, 'r').readlines() 
                 if 'ATOM' == line[:4]]
        prot = [line for line in lines if any(aa in line for aa in self.aa)]
        return prot
    
    @property
    def amino_acids(self):
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
                outline = PDB.format_line(line)
                outfile.write(outline)
                
        return fname
                
    def merge_final_pdb(self, new_coordinates: Dict[str, List[str]]) -> None:
        """_summary_

        Args:
            new_coordinates (Dict[str, List[str]]): _description_

        Returns:
            _type_: _description_
        """
        with open(f'{self.output_path}/{self.outname}.pdb', 'w') as outfile:
            pass


class rtfParser:
    """
    This class parses CHARMM rtf parameter files to obtain lipid bonded
    parameters to be used in cryo-em lipid construction.
    """
    def __init__(self, lipids: List[str] = ['POPE','POPC','POPG','POPS',
                                            'POPI24','CHL1','PVCL2']):
        self.lipids = lipids
        self._rtfs = self.get_rtfs


    @property
    def get_rtfs(self) -> dict[str]:
        _files = {'top_all36_lipid.rtf': ['POPE','POPC','POPG','POPS'],
                 'toppar_all36_lipid_inositol.str': ['POPI24'],
                 'toppar_all36_lipid_cholesterol.str': ['CHL1'],
                 'toppar_all36_lipid_bacterial.str': ['PVCL2']}

        all_rtfs = {}
        for f in _files.keys():
            all_rtfs.update(rtfParser.unpack_rtf(f))
        
        rtfs = {lip: rtfParser.deconvolute_ic_table(all_rtfs[lip]) for lip in self.lipids}

        return rtfs


    @staticmethod
    def unpack_rtf(_file: str) -> Dict[str, List[str]]:
        with open(f'rtf_files/{_file}', 'r') as rfile:
            ics = {}
            for line in rfile:
                if line[:4] == 'RESI':
                    resi = line[5:11].strip()
                    ics.update({resi: []})
                elif line[:2] == 'IC':
                    ics[resi].append(line.strip())

        return ics


    @staticmethod
    def deconvolute_ic_table(ic_table: List[str]) -> Dict[str, Dict[str, List[float]]]:
        bonds = {}
        angles = {}
        dihedrals = {}
        impropers = {}
        for line in ic_table:
            _, a1, a2, a3, a4, *params = [ele.strip() for ele in line.split()]
            # handle case of ! defines S chirality comment
            if len(params) > 5:
                params = params[:5]

            params = [float(param) for param in params]
            # handle improper dihedral notation
            if '*' in a3:
                a3 = a3[1:]
                bonds.update({f'{a1}-{a3}': params[0], f'{a3}-{a4}': params[-1]})
                impropers.update({f'{a1}-{a2}-{a3}-{a4}': params[2]})

            else:
                bonds.update({f'{a1}-{a2}': params[0], f'{a3}-{a4}': params[-1]})
                dihedrals.update({f'{a1}-{a2}-{a3}-{a4}': params[2]})

            angles.update({f'{a1}-{a2}-{a3}': params[1], f'{a2}-{a3}-{a4}': params[-2]})

        return {'bond': bonds, 'angle': angles, 'dihedral': dihedrals, 'improper': impropers}

def unpack_lipids(cfg: Dict[str, Dict[str, str]]) -> tuple(List[str]):
    resids, restypes, resnames = list(), list(), list()
    for val in cfg['lipids'].values():
        resids.append(val['resid'])
        restypes.append(val['restype'])
        resnames.append(val['current_resname'])
        
    return resids, restypes, resnames