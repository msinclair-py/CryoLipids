#!/usr/bin/env python
from copy import deepcopy
from molecular_graph import MolecularGraph
import networkx as nx
import numpy as np
from scipy.spatial.transform import Rotation
from typing import Dict, List, Union
from utilities import PDB

class Lipid(PDB):
    """
    Class for modeling a single lipid of resID `resid` and resType `restype`
    using an input CHARMM IC table and connectivity graph.
    """
    def __init__(self, pdbfile: str, resid: int,
                    restype: str,
                    current_resname: str):
        super().__init__(pdbfile, '', [resid], old_resnames=current_resname)
        self.lipid_type = restype
        self.pdb_contents = self.lipid_lines
        self.graph = MolecularGraph(restype).G
    
    @property
    def lipid_lines(self) -> List[str]:
        """
        Overloaded `contents` method to ensure we don't populate an empty list
        for `self.pdb_contents`.
        
        Returns:
            List[str]: Contents of incompletely modeled lipid for `self.resids[0]` resid
        """
        pdb = open(self.filename).readlines()
        parsed = self.parse_pdb(pdb, self.resids[0], self.resnames)
        
        for i in range(len(parsed)):
            parsed[i][3] = self.lipid_type
            
        return parsed
    
    def extract_coordinates(self) -> np.ndarray:
        """
        Extract just the coordinates from the `self.pdb_contents` list
        
        Returns:
            np.ndarray: Coordinate array of shape (N, 3)
        """
        coords = np.zeros((len(self.pdb_contents), 3))
        for i, line in enumerate(self.pdb_contents):
            coords[i, :] = line[6:9]
        return coords.astype(np.float64)
    
    def add_to_pdb(self, atom_names: List[str], coords: np.ndarray) -> None:
        """
        Adds new lines to internal representation of pdb stored in `self.pdb_contents`.
        
        Args:
            atom_names (List[str]): List of new atom names. Should correspond in index to coords
                                        found in the `coords` input array.
            coords (np.ndarray): Numpy array of XYZ coordinates for each atom in `atom_names`.
        """
        last_num = int(self.pdb_contents[-2][1])
        chain = self.pdb_contents[-2][4].strip()
        
        lines = []
        for name, coord in zip(atom_names, coords):
            last_num += 1
            line = ['ATOM', last_num, name, self.lipid_type, chain, self.resids[0],
                    *coord, 0.0, 0.0]
            lines.append(line)
            
        self.pdb_contents += lines
        
    def model(self, rotate_along_z: int=2) -> None:
        """
        Main lipid modeling function. 
        (i) Identifies terminally modeled atoms
        (ii) Generates list of what atoms are left to model
        (iii) Aligns an example of a complete lipid tail based on the 
            final two modeled atoms
        (iv) Rotates the newly modeled atoms such that a vector going from 
            new_atom_0 -> new_atom_n is aligned along the z axis ({0, 0, 1} 
            or {0, 0, -1} based on whether a majority of atoms are above or 
            below the phosphorous atom)
        (v) Store new system coordinates for later extraction
        
        Args:
            rotate_along_z (int): Numerical cutoff for the number of atoms in 
                                    a given chain to require alignment along
                                    the z-axis. Short fragments should not be
                                    perturbed in this manner (think phosphate
                                    oxygen for example). Default length of >2 atoms.
        """
        # first, find any missing nodes from the lipid graph
        missing_atoms = self.get_missing_atoms()
        
        # second, group missing atoms into continuous chains of atoms which are bonded
        missing_chains = self.get_missing_chains(missing_atoms)
        
        # third, identify vector to align/rotate and coordinates to perform rotation on
        for missing_chain in missing_chains:
            prev_atom1, *prev_atom2 = self.get_previous_atoms(missing_chain)
            vector_ref = np.vstack((self.get_coord(prev_atom2), self.get_coord(prev_atom1)))
            
            # fourth, staple on each group of missing atoms
            template_lipid = Template(f'lipids_from_rtf/{self.lipid_type}.pdb', 
                                      self.lipid_type)
            
            coords_to_rotate = template_lipid.atomic_coordinates(missing_chain)
            vector_align = template_lipid.atomic_coordinates(prev_atom2 + [prev_atom1])
            new_tail_coords = self.staple_fragment(vector_ref, vector_align, coords_to_rotate)
            
            # fifth, align placed group with z-axis
            if len(missing_chain) > rotate_along_z:
                new_tail_vec = np.vstack((new_tail_coords[0, :], new_tail_coords[-1, :]))
                
                phosphate_coord = self.get_coord('P')[-1]
                lipid_cog_z = np.mean(np.concatenate((self.get_coord(), new_tail_coords), axis=0), 
                                      axis=1)[-1]
                
                if lipid_cog_z > phosphate_coord:
                    z = np.array([0, 0, 1])
                else:
                    z = np.array([0, 0, -1])
                    
                align_vec = np.vstack((new_tail_vec[0], new_tail_vec[0] + z))
                new_tail_coords = self.staple_fragment(align_vec, new_tail_vec, new_tail_coords)
            
            self.add_to_pdb(missing_chain, new_tail_coords)
    
    def get_previous_atoms(self, chain: List[str]) -> List[str]:
        """
        Given a chain of continously connected atoms, return the preceding two
        atoms by atom name. If the first preceding atom is `C2` that means we have
        reached the top of the graph, and we need to take care to align the missing
        tail on by ensuring we align not only to the C1-C2 bond but also the other
        tail bond, either C2-O21 or C2-C3. If both tails are missing we will handle
        that edge case in the stapling method.
        
        Args:
            chain (List[str]): Atom names for given chain

        Returns:
            tuple[str]: (i-1 atom, i-2 atom)
        """
        prev1 = list(self.graph.predecessors(chain[0]))[0]
        if prev1 == 'C2':
            prev2 = ['C1']
            if chain[0] == 'C3':
                prev2 += ['O21']
            else:
                prev2 += ['C3']
        else:
            prev2 = list(self.graph.predecessors(prev1))
            
        return [prev1] + prev2

    def get_missing_chains(self, missing_atoms: List[str]) -> List[List[str]]:
        """
        Using the internal graph representation of lipid, identify continously
        connected chains of missing atoms. Topological sort allows us to then
        ensure these chains are sorted in order from the top of the graph 
        (e.g. closest to the `C2` atom) down.

        Args:
            missing_atoms (List[str]): List of missing atoms by atom name

        Returns:
            List[str]: List of lists of chains of connected missing atoms
        """
        connections = []
        sort_dict = {val: i for i, val in enumerate(nx.topological_sort(self.graph))}
        for i, node in enumerate(missing_atoms):
            if any([node in connect for connect in connections]):
                continue

            neighbors = list(self.graph.neighbors(node)) + list(self.graph.predecessors(node))
            
            other_nodes = deepcopy(missing_atoms)
            other_nodes.pop(i)
            current_connections = [node]
            while True:
                for j, n in enumerate(other_nodes):
                    if n in neighbors:
                        current_connections += [n]
                        other_nodes.pop(j)
                        neighbors += list(self.graph.neighbors(n)) + list(self.graph.predecessors(n))
                        break
                else:
                    break
                
            connections.append(sorted(current_connections, key=lambda k: sort_dict.get(k)))
        
        return connections
    
    @staticmethod
    def staple_fragment(v1: np.ndarray, v2: np.ndarray, 
                        fragment: np.ndarray) -> np.ndarray:
        """
        Given the coordinates of the terminal atoms in a given fragment, the
        template coordinates for the same two atoms, and the template coordinates
        for all missing atoms in the fragment being modeled we align the missing
        coordinates using the Kabsch algorithm.
        
        Args:
            v1 (np.ndarray): Array of the coordinates of two atoms comprising
                                the preceding bond we are stapling onto. Note
                                that these atoms are the current terminal atoms
                                in our incomplete lipid. 
            v2 (np.ndarray): Array of the coordinates of the same two atoms as
                                `v1` but coming from the complete, template lipid.
            fragment (np.ndarray): Coordinates of all missing atoms for a given fragment.
                                Originate from the template lipid.
        
        Returns:
            np.ndarray: Transformed coordinates of the missing atoms coming from `fragment`.
                            Because the atoms in `v1` should not have moved they are not
                            present in this array
        """
        to_origin = v2[0]
        from_origin = v1[0]
        
        rotmatrix = Rotation.align_vectors((v1 - v1[0]), (v2 - v2[0]))[0]
        
        return rotmatrix.apply(fragment - to_origin) + from_origin
        
    def get_missing_atoms(self) -> List[str]:
        """
        Get the atom names for all missing atoms. Uses an internal graph representation
        of a complete lipid to check this.

        Returns:
            List[str]: Names of the missing atoms
        """ 
        atom_names = [atom[2].strip() for atom in self.pdb_contents]
        missing_atoms = [atom for atom in self.graph.nodes if atom not in atom_names]
            
        return missing_atoms
        
    def get_coord(self, name: Union[List[str], str, None] = None) -> np.ndarray:
        """
        Gets coordinate(s) of input atom(s) from internal pdb representation.
        
        Args:
            name (Union[str, List[str]]): Either a single name (str) or a list of atom
                                            names (List[str]) to grab the coordinates of

        Raises:
            ValueError: If an atom of name `name` is not found in `self.pdb_contents`
                            throw an error. Either the code is misbehaving or we have
                            passed an illegal name, but either way we can't proceed

        Returns:
            np.ndarray: Coordinates stored in array of shape (N, 3)
        """
        if isinstance(name, list):
            coords = np.zeros((len(name), 3))
            for i, n in enumerate(name):
                for atom in self.pdb_contents:
                    if atom[2].strip() == n:
                        try:
                            coords[i] = [float(x.strip()) for x in atom[6:9]]
                        except AttributeError:
                            coords[i] = atom[6:9]
                            
                        break
            return coords
                
        elif isinstance(name, str):
            for atom in self.pdb_contents:
                if atom[2].strip() == name:
                    return np.array([float(x.strip()) for x in atom[6:9]])
            
        else:
            coords = np.zeros((len(self.pdb_contents), 3))
            for i, atom in enumerate(self.pdb_contents):
                try:
                    coords[i] = [float(x.strip()) for x in atom[6:9]]
                except AttributeError:
                    coords[i] = atom[6:9]
                    
            return coords
            
        raise ValueError(f'Atom {name} not found in partial lipid!')
    
    def update_coords(self, names: List[str], coords: np.ndarray) -> None:
        """
        Updates internal representation of coordinates found in `self.pdb_contents`

        Args:
            names (List[str]): Names of atoms for which to update the coordinates
            coords (np.ndarray): Array of coordinates
        """
        for (name, coord) in zip(names, coords):
            idx = np.where(np.array(self.pdb_contents, ndmin=2)[:, 2] == name)[0][0]
            self.pdb_contents[idx][6:9] = coord


class Template(PDB):
    """
    Child of `PDB` class used to access lipid coordinates of example
    rtf lipids and tie coordinate data to connected fragment data.
    """
    def __init__(self, pdbfile: str, restype: str):
        super().__init__(pdbfile, '', resids=[1], old_resnames=[restype])
        self.atoms = self.contents
        self.filter_heavy_atoms()

    
    def filter_heavy_atoms(self) -> None:
        """
        Filter out all hydrogens in model, ensuring we only have heavy atoms.
        """
        heavy = []
        for atom in self.atoms:
            heavy.append(atom)
            #if 'H' not in atom[2].strip():
            #    heavy.append(atom)
        
        self.heavy = Template.process(heavy)


    @staticmethod
    def process(atoms: List[str]) -> Dict[str, np.ndarray]:
        """
        Obtain coordinates for a set of atoms from the template system.

        Args:
            atoms (List[str]): Atoms for which to obtain coordinates

        Returns:
            Dict[str, np.ndarray]: Dictionary for atom:coords
        """
        processed = dict()
        for atom in atoms:
            coord = np.array([i.strip() for i in atom[6:9]], dtype=np.float64)
            processed[atom[2].strip()] = coord

        return processed

    def atomic_coordinates(self, names: List[str]) -> np.ndarray:
        """
        Extract atomic coordinates from internal heavy atom dictionary

        Args:
            names (List[str]): Names of atoms

        Returns:
            np.ndarray: Coordinates of atoms
        """
        return np.array([self.heavy[name] for name in names], dtype=np.float64)
    
    def missing_atoms(self, last_atom: str) -> np.ndarray:
        tail = last_atom[:2]
        num = int(last_atom[2:])
        all_atoms = [key for key in self.heavy.keys() if key[:2] == tail]
        
        remaining = []
        for atom in all_atoms:
            try:
                if int(atom[2:]) > num:
                    remaining.append(atom)
            except ValueError:
                continue

        return remaining, self.atomic_coordinates(remaining)
