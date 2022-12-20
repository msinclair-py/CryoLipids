#!/usr/bin/env python
from copy import deepcopy
from itertools import combinations
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

        branch_residues = {
                'POPE': lip,
                'POPC': lip.copy().update({'N': ['C13', 'C14']}),
                'POPG': lip.copy().update({'C12': ['OC2']}),
                'POPS': lip.copy().update({'C12': ['N']}),
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


    def dfs(self, graph: Dict[str, Set[str]], node: str,
            path: List[str] = None, visited: Set[str] = set()) -> List[str]:
        """
        Recursive depth-first search algorithm. Returns list of all identified
        paths in `graph`.
        """
        if not path: path = [node]
        visited.add(node)

        paths = []
        for t in graph[node]:
            if t not in visited:
                t_path = path + [t]
                paths.append(t_path)
                paths.extend(self.dfs(graph, t, t_path, visited))

        return paths

    
    @property
    def longest_path(self):
        """
        Performs dfs and then returns the longest pathway.
        """
        paths = self.dfs(self._graph, self.root)
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

    
    def get_diagrams(fragments):
        return diagrams


    ### OBSOLETE ###
    @property
    def adjacency_list(self) -> Dict[str, Dict[str, str]]:
        adj = {}
        for (e1, e2) in self._edges:
            if any(['H' in e1, 'H' in e2]):
                continue

            try:
                adj[e1] = adj[e1] + [e2]
            except KeyError:
                adj[e1] = [e2]

            try:
                adj[e2] = adj[e2] + [e1]
            except KeyError:
                adj[e2] = [e1]

        return adj

    ### OBSOLETE ### BRUTE FORCE CITY
    @staticmethod
    def matt_algorithm(adj_list: Dict[str, List[str]], N: int) -> Dict[str, List[str]]:
        all_nodes = set(adj_list.keys())
        fragments = dict()
        for k in range(3, N):
            print(f'\nGenerating fragments of length {k}!')
            k_combinations = combinations(all_nodes, k)
            for combination in k_combinations:
                if MolecularGraph.path_exists(adj_list, combination):
                    try:
                        fragments[k].append(combination)
                    except KeyError:
                        fragments[k] = list(combination)

        return fragments


    @staticmethod
    def path_exists(adjacency: Dict[str, List[str]], path: List[str]) -> bool:
        neighbors = []
        for element in path:
            for val in adjacency[element]:
                if val in path:
                    neighbor = {element, val}
                    if neighbor not in neighbors:
                        neighbors.append(neighbor)

        if len(neighbors) >= len(path) - 1: # can be used to detect rings if there
            #print(neighbors, path)          # are `len(path)` number of connections
            return True                     # in neighbors

        return False

    
    ### OBSOLETE ###
    @staticmethod
    def has_path(graph: Dict[str, List[str]], src: str, 
                   dst: str, visited: Set[str] = set()) -> None:

        if src in visited:
            return False
        else:
            visited.add(src) # critical for not getting trapped
                             # in infinite recursion
        if src == dst:
            return True
       
        for neighbor in graph[src]:
            if MolecularGraph.has_path(graph, neighbor, dst, visited):
                return True

        return False


    ### OBSOLETE ###
    @staticmethod
    def depth_first_traversal(all_nodes: Set[str], 
                               adj_list: Dict[str, List[str]]) -> List[List[str]]:
        paths = []
        for s, d in combinations(all_nodes, 2):
            path = MolecularGraph.has_path(adj_list, s, d)
            print(path)

        return paths


    ### OBSOLETE ###
    def identify_subgraphs(all_nodes: Set[str]):
        potential = {N: combinations(all_nodes, N) for N in range(len(all_nodes))}
        for path in potential:
            dsf
        return 0
