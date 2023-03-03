#!/usr/bin/env python
from copy import deepcopy
from itertools import combinations
import numpy as np
import persim
import ripser
from typing import Dict, List, Set

class MolecularGraph:
    """
    Collection of undirected graph methods for graph construction and traversal of atomic 
    information in rtf files. Intended for performing path matching and fragment library
    generation.
    """
    
    def __init__(self, rtf: Dict[str, Dict[str, float]], lipid: str):
        self._rtfs = rtf
        self._lipid = lipid


    @property
    def root(self) -> str:
        """
        Root node for given lipids to obtain the longest path through each
        molecule.
        """
        roots = {'POPE': 'N',
                 'POPC': 'C15',
                 'POPG': 'OC3',
                 'POPS': 'O13A',
                 'POPA': 'O12'}

        return roots[self._lipid]


    @property
    def branches(self) -> Dict[str, List[str]]:
        """
        Returns atoms where branching occurs in the longest continuous path
        through each given lipid molecule.
        """
        lip = {'C2': ['C3', 'O31'] + [f'C3{n}' for n in range(1, 17)],
               'P': ['O13', 'O14'],
               'C21': ['O22'],
               'C31': ['O32']}
        
        pc, pg, ps = [deepcopy(lip) for _ in range(3)]
        pc.update({'N': ['C13', 'C14']})
        pg.update({'C12': ['OC2']})
        ps.update({'C12': ['N']})

        branch_residues = {
                'POPE': lip,
                'POPC': pc,
                'POPG': pg,
                'POPS': ps,
                'POPA': lip
                }

        return branch_residues[self._lipid]


    @staticmethod
    def sliding_window(vec: List[str], n: int) -> List[List[str]]:
        """
        Create list of lists where each list corresponds to full sliding window sweep of
        size `k` where 3 <= k <= n.
        """
        windows = []
        for i in range(len(vec) - n + 1):
            windows.append(vec[i:i+n])

        return windows


    def complete_library(self,
            fragments: Dict[int, List[List[str]]]) -> Dict[int, List[List[str]]]:
        """
        Take fragments from longest pathway in structure and divergent paths (sn-2 tail, phosphate
        oxygens, etc.) to generate full fragment library.
        """
        ###### NOTE: This should be refactored for efficiency
        branch_map = self.branches
        for atom, frags in branch_map.items():
            matches = []
            for p in fragments.values():
                for v in p:
                    if atom in v:
                        matches.append(v)

            match atom:
                case 'C2':
                    new_paths = [frags[:n] for n in range(1, len(frags) + 1)]
                case 'P' | 'N':
                    new_paths = [[frag] for frag in frags] + [frags]
                case _:
                    new_paths = [frags] if isinstance(frags[0], str) else frags
           
            new_frags = dict()
            for m in matches:
                for p in new_paths:
                    temp = deepcopy(m)
                    temp.extend(p)
                    
                    try:
                        if temp not in new_frags[len(temp)]:
                            new_frags[len(temp)].append(temp)
                    except KeyError:
                        new_frags[len(temp)] = [temp]
            
            for length, path in new_frags.items():
                for p in path:
                    try:
                        fragments[length].append(p)
                    except KeyError:
                        fragments[length] = [p]

        return fragments


    def fragment_lipid(self, longest_path: List[str]) -> Dict[int, List[List[str]]]:
        """
        Take longest pathway through structure and obtain all fragments of size 3:N. Pass
        this initial fragment library to `complete_library` to obtain the full fragment
        library for the lipid type of instantiated class.
        """
        fragments = dict()
        for n in range(3, len(longest_path)):
            fragments[n] = self.sliding_window(longest_path, n)
        
        return self.complete_library(fragments)


    def dfs(self, node: str, path: List[str] = None, 
            visited: Set[str] = set()) -> List[str]:
        """
        Recursive depth-first search algorithm. Returns list of all identified
        paths in `graph`.
        """
        if not path: path = [node]
        visited.add(node)

        paths = []
        for tail in self._graph[node]:
            if tail not in visited:
                tail_path = path + [tail]
                paths.append(tail_path)
                paths.extend(self.dfs(tail, tail_path, visited))

        return paths

    
    @property
    def longest_path(self):
        """
        Performs dfs and then returns the longest pathway.
        """
        paths = self.dfs(self.root, visited=set())
        longest = tuple()
        for path in paths:
            if len(path) > len(longest):
                longest = path

        return longest


    @staticmethod
    def generate_graph(rtf_dict: Dict[str,float]) -> Dict[str, Set[str]]:
        graph = {}
        for key in rtf_dict['bond'].keys():
            a1, a2 = key.split('-')
            if not any(['H' in a for a in (a1, a2)]):
                try:
                    graph[a1].add(a2)
                except KeyError:
                    graph.update({a1: {a2}})

                try:
                    graph[a2].add(a1)
                except KeyError:
                    graph.update({a2: {a1}})

        return graph


    def generate_edges(self, graph: Dict[str, Set[str]]):
        edges = []
        for node in graph:
            for neighbor in graph[node]:
                if not any(x in edges for x in [{node, neighbor}, 
                                                {neighbor, node}]):
                    edges.append({node, neighbor})

        self._edges = edges


    @property
    def connectivity_graphs(self) -> Dict[str, Dict[str, str]]:
        graph = self.generate_graph(self._rtfs[self._lipid])
        self.generate_edges(graph)
        
        self._graph = graph
        return graph


    def get_coords(self, fragment: List[str]) -> np.ndarray:
        coords = np.array([])
        for line in self._rtfs[self._lipid]['bond']:
            print(line)
        return coords


class PersistentHomology:
    """
    Class for obtaining persistence diagrams and performing Wasserstein
    distance analysis on one point cloud versus a set of point clouds.
    Identifies the best possible match.
    """

    def __init__(self, cloud: List[List[float]], 
                    cloud_library: Dict[str, np.ndarray]):
        self.cloud = PersistentHomology.get_diagram(np.array(cloud))
        self.cloud_library = cloud_library


    @staticmethod
    def get_diagram(data: np.ndarray) -> List[np.ndarray]:
        """
        Using algebraic topology, obtain a persistence diagram
        of each fragment for comparison with experimentally-derived
        lipid fragments.
        """
        return ripser.ripser(data)['dgms']


    @property
    def best_fit(self) -> int:
        lowest = [0, 1e6]
        for idx, diagram in self.cloud_library.items():
            print(f'{self.cloud}\n\n{diagram}')
            distance = persim.wasserstein(self.cloud, diagram)
            if distance < lowest[1]:
                lowest = [idx, distance]

        return lowest[0]
