#!/usr/bin/env python
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
                    ic_table: Dict[str, Dict[str, float]],
                    current_restype: str = 'POV'):

       super().__init__(pdbfile, [resid], resname=current_restype)
       self.ic_table = ic_table
       self.pdb_contents = self.contents
       self.unmodeled = self.pdb_contents
    
    def extract_coordinates(self) -> np.ndarray:
        pdb_contents = np.asarray(self.contents)
        return pdb_contents[:,6:9].astype(np.float64)


    def model(self):
        # identify terminated sequences

        # prioritize modeling order

        # model atoms one at a time

        # this method will perform first-pass modeling only
        return 1


    def _add_atom(self, atom_name: str) -> None:
        pass


    @staticmethod
    def vectorize(*args) -> List[np.ndarray]:
        assert len(args)%2 == 0,"ERROR: Must provide even number of points!"
        vecs = [args[2*x] - args[2*x+1] for x in range(int(len(args)/2))]
        return [vec/np.linalg.norm(vec) for vec in vecs]


    @staticmethod
    def measure_angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
        v1, v2 = Lipid.vectorize(a, b, c, b)
        return np.arccos(v1 @ v2)


    @staticmethod
    def measure_dihedral(a: np.ndarray, b: np.ndarray,
                         c: np.ndarray, d: np.ndarray) -> float:
        v1, v2, v3 = Lipid.vectorize(b, a, c, b, d, c)
        norm1 = np.cross(v1, v2)
        norm2 = np.cross(v2, v3)
        sign = np.sign(v1 @ v3)

        rad = np.arccos(norm1 @ norm2)

        if sign != 0:
            rad *= sign
    
        return rad * 180 / np.pi


    @staticmethod
    def measure_improper(a: np.ndarray, b: np.ndarray,
                         c: np.ndarray, d: np.ndarray) -> float:
        v1, v2, v3 = Lipid.vectorize(a, b, d, b, c, b)
        plane_normal = np.cross(v1, v2)
        return 90 - np.arccos(plane_normal @ v3) * 180 / np.pi


    @staticmethod
    def rotate_angle(a: np.ndarray, b: np.ndarray,
                     c: np.ndarray, angle: int) -> List[np.ndarray]:
        rotation_matrix = Rotation.from_euler('', angle, degrees=True)
        rotated_array = rotation_matrix.apply(array)
        return rotated_array


    def rotate_dihedral():
        pass


    def rotate_improper():
        pass
    
    @staticmethod
    def staple_tail(v1, v2, arr):
        """
        Given 2 arrays of 2 points each, which define our two vectors,
        generate the rotation matrix which aligns v2 onto v1 both in terms
        of angle and translation using the Kabsch algorithm. Apply this
        rotation onto the array `arr`.
        """ 
        rotmatrix = Rotation.align_vectors((v1[1] - v1[0]).reshape(1, 3), 
                                           (v2[1] - v2[0]).reshape(1, 3))[0]
        return rotmatrix.apply(arr)
        
    
    def get_terminal_atoms(self):
        """
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
        for atom in self.contents:
            if atom[2].strip() == name:
                return np.array([float(x.strip()) for x in atom[6:9]])
            
        raise ValueError(f'Atom {name} not found in partial lipid!')


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