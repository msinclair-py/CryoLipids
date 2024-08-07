#!/usr/bin/env python
from copy import deepcopy
import numpy as np
from scipy.spatial.transform import Rotation
from typing import Dict, List
from utilities import PDB

class Lipid(PDB):
    """
    Class for modeling a single lipid of resID `resid` and resType `restype`
    using an input CHARMM IC table and connectivity graph.
    """
    def __init__(self, pdbfile: str, resid: int,
                    lipid_type: str,
                    current_restype: str = 'POV',
                    collision: int = 0):

       super().__init__(pdbfile, [resid], resname=current_restype)
       self.pdb_contents = self.contents
       self.lipid_type = lipid_type
       self.unmodeled = deepcopy(self.pdb_contents)
       self.collision_detector = None
       self.collision = collision
    
    def extract_coordinates(self) -> np.ndarray:
        coords = np.zeros((len(self.pdb_contents), 3))
        for i, line in enumerate(self.pdb_contents):
            coords[i, :] = line[6:9]
        return coords.astype(np.float64)
    
    def add_to_pdb(self, atom_names: List[str], coords: np.ndarray) -> None:
        last_num = int(self.pdb_contents[-2][1])
        chain = self.pdb_contents[-2][4].strip()
        
        lines = []
        for name, coord in zip(atom_names, coords):
            last_num += 1
            line = ['ATOM', last_num, name, self.resname, chain, self.resids[0],
                    *coord, 0.0, 0.0]
            lines.append(line)
            
        self.pdb_contents += lines
        

    def model(self):
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
        (v) Checks for and repairs protein-lipid clashes
        """
        terminal_atoms = self.get_terminal_atoms()
        tail_map = {'sn1': 'C2', 'sn2': 'C3'}
        for (tail, terminus) in terminal_atoms.items():
            cur = f'{tail_map[tail]}{list(terminus.keys())[0]}'
            prev = f'{tail_map[tail]}{list(terminus.keys())[0] - 1}'
            prev_coords = self.get_coord(prev)
            vector_ref = np.array([list(terminus.values())[0], prev_coords])
            rtf_lipid = Template(f'lipids_from_rtf/{self.lipid_type}.pdb', 
                                 self.lipid_type)
            vector_comp = rtf_lipid.atomic_coordinates([cur, prev])
            new_tail_names, new_tail_coords = rtf_lipid.missing_atoms(cur)
            new_tail_coords = self.staple_tail(vector_ref, vector_comp, new_tail_coords)
            new_tail_vec = [new_tail_coords[0, :], new_tail_coords[-1, :]]
            
            # choose z based on whether a majority of tail atoms are above or below phosphate
            phosphate_coord = self.get_coord('P')
            tail_coords = self.get_coord([f'{tail_map[tail]}{i}' 
                                          for i in range(1, list(terminus.keys())[0])])
            
            num_above_phos = sum([1 for c in tail_coords if c[-1] > phosphate_coord[-1]])
            if num_above_phos / tail_coords.shape[0] >= 0.5:
                    z = np.array([0, 0, 1])
            else:
                    z = np.array([0, 0, -1])

            align_vec = np.vstack((new_tail_vec[0], new_tail_vec[0] + z))
            new_tail_coords = self.staple_tail(align_vec, new_tail_vec, new_tail_coords)
            
            self.add_to_pdb(new_tail_names, new_tail_coords)
            self.write_to_pdb_file(self.pdb_contents)
    
    @staticmethod
    def staple_tail(v1: np.ndarray, v2: np.ndarray, 
                    arr: np.ndarray) -> np.ndarray:
        """
        Given 2 arrays of 2 points each, which define our two vectors,
        generate the rotation matrix which aligns v2 onto v1 both in terms
        of angle and translation using the Kabsch algorithm. Apply this
        rotation onto the array `arr`.
        """ 
        center = v2[0]
        translation = v1[0]
        
        rotmatrix = Rotation.align_vectors((v1[1] - v1[0]).reshape(1, 3), 
                                           (v2[1] - v2[0]).reshape(1, 3))[0]
        
        return rotmatrix.apply(arr - center) + translation
        
    
    def get_terminal_atoms(self):
        """
        NOTE: Integrate networkx lipid model to identify missing atoms.
        
        Identify any incompleteness in lipid molecule. Should return
        the terminal headgroup atom, and terminal tail atoms which are
        where template lipids will be attached to complete.
        """
        highest_sn1 = [1, []]
        highest_sn2 = [1, []]
        
        for atom in self.contents:
            atom_name = atom[2].strip()
            try:
                atom_num = int(atom_name[2:])
                if 'C2' in atom_name and atom_num > highest_sn1[0]:
                    highest_sn1 = [atom_num, np.array([float(x.strip()) for x in atom[6:9]])]
                elif 'C3' in atom_name and atom_num > highest_sn2[0]:
                    highest_sn2 = [atom_num, np.array([float(x.strip()) for x in atom[6:9]])]
            except ValueError:
                continue
        
        sn1_dict = {highest_sn1[0]: highest_sn1[1]}
        sn2_dict = {highest_sn2[0]: highest_sn2[1]}
        return {'sn1': sn1_dict, 'sn2': sn2_dict}
    
        
    def get_coord(self, name):
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
                
        else:
            for atom in self.contents:
                if atom[2].strip() == name:
                    return np.array([float(x.strip()) for x in atom[6:9]])
            
        raise ValueError(f'Atom {name} not found in partial lipid!')
    
    def update_coords(self, names: List[str], coords: np.ndarray) -> None:
        """Does this even work??

        Args:
            names (List[str]): _description_
            coords (np.ndarray): _description_

        Raises:
            ValueError: _description_
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
        super().__init__(pdbfile, resname=restype)
        self.atoms = self.contents
        self.filter_heavy_atoms()

    
    def filter_heavy_atoms(self) -> None:
        heavy = []
        for atom in self.atoms:
            if 'H' not in atom[2].strip():
                heavy.append(atom)
        
        self.heavy = Template.process(heavy)


    @staticmethod
    def process(atoms: List[str]) -> Dict[str, np.ndarray]:
        processed = dict()
        for atom in atoms:
            coord = np.array([i.strip() for i in atom[6:9]], dtype=np.float64)
            processed[atom[2].strip()] = coord

        return processed


    def atomic_coordinates(self, names: List[str]) -> np.ndarray:
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
