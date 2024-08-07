#!/usr/bin/env python
from copy import deepcopy
import networkx as nx
import numpy as np
from typing import List

class MolecularGraph:
    """_summary_

    Returns:
        _type_: _description_
    """
    def __init__(self, lipid: str):
        edges = self.parse_rtf_file(lipid)
        self.graph = self.generate_graph(edges)
        
    @staticmethod
    def parse_rtf_file(lipid: str) -> List[str]:
        match lipid:
            case 'POPE':
                fname = 'top_all36_lipid.rtf'
            case 'POPC':
                fname = 'top_all36_lipid.rtf'
            case 'POPG':
                fname = 'top_all36_lipid.rtf'
            case _:
                raise NameError(f'Lipid type {lipid} NOT found in CHARMM forcefield!')
            
        rtf = open(f'rtf_files/{fname}').readlines()
        
        bond_lines = []
        dump = False
        for line in rtf:
            if dump and line[:4] == 'RESI':
                break
            
            if line[:9] == f'RESI {lipid}':
                dump = True
                
            if dump:
                if line[:4] == 'BOND' or line[:6] == 'DOUBLE':
                    bond_lines.append(line[6:].strip().split())
                    
        bonds = []
        for bond_line in bond_lines:
            n = len(bond_line) // 2
            i = 0
            while i < n:
                bond = bond_line[i*2:(i+1)*2]
                if not any(['H' in ele for ele in bond]):
                    bonds.append((bond))
                i += 1
                
        return np.array(bonds, ndmin=2)
    
    def generate_graph(self, bonds: np.ndarray) -> None:
        parent_nodes = ['C2']
        G = nx.DiGraph()
        G.add_node(parent_nodes[0])

        while True:
            new_nodes = []
            for parent_node in parent_nodes:
                
                idxs = np.where(bonds == parent_node)
                for (idx, j) in zip(idxs[0], idxs[1]):
                    jdx = abs(j - 1)
                    # identify atom that is not parent_node
                    new_atom = bonds[idx, jdx]
                    # add new edge
                    G.add_edge(parent_node, new_atom)
                    # populate new nodes list
                    new_nodes.append(new_atom)
                
                # remove bond(s) from bonds array
                bonds = np.delete(bonds, idxs[0], axis=0)

            if not new_nodes:
                break
            else:
                parent_nodes = deepcopy(new_nodes)
        
        return G