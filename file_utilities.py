#!/usr/bin/env python
import os
from typing import Dict, List, Union, Set

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
        fname = os.path.join(self.output_path, f'processed_{os.path.basename(self.filename)}')
        
        with open(fname, 'w') as outfile:
            for line in contents:
                outline = PDB.format_line(line)
                outfile.write(outline)


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


    @staticmethod
    def generate_graph(rtf_dict: Dict[str,float]) -> Dict[str, Set[str]]:
        graph = {}
        for key in rtfs['POPE']['bond'].keys():
            a1, a2 = key.split('-')
            try:
                graph[a1].add(a2)
            except KeyError:
                graph.update({a1: {a2}})

        return graph


    @staticmethod
    def generate_edges(graph: Dict[str, Set[str]]) -> List[str]:
        edges = []
        for node in graph:
            for neighbor in graph[node]:
                if not any(x in edges for x in [{node, neighbor}, 
                                                {neighbor, node}]):
                    edges.append({node, neighbor})

        return edges


    @property
    def connectivity_graphs(self) -> Dict[str, Dict[str, str]]:
        edges = {}
        for lipid in self._rtfs:
            graph = self.generate_graph(self._rtfs[lipid])
            edges.update({lipid: self.generate_edges(graph)})

        return edges
