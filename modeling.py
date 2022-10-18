#!/usr/bin/env python
from typing import Dict, List
from file_utilities import PDB
import numpy as np

class Lipid(PDB):
    """
    Class for modeling a single lipid of resID `resid` and resType `restype`
    using an input CHARMM IC table.
    """
    def __init__(self, pdbfile: str, resid: int, restype: str = 'POV', 
                        ic_table: Dict[str, Dict[str, float]]):
        super().__init__(pdbfile, [resid], resname=restype)
        self.ic_table = ic_table
        self.pdb_contents = self.contents


    def model(self):
        pass


    @staticmethod
    def vectorize(*args) -> List[np.ndarray]:
        assert len(args)%2 = 0,"ERROR: Must provide even number of points!"
        vecs = [args[2*x] - args[2*x+1] for x in range(int(len(args)/2))]
        return [vec/np.linalg.norm(vec) for vec in vecs]


    @staticmethod
    def measure_angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
        v1, v2 = Lipid.vectorize(a, b, c, b)
        return np.arccos(v1 @ v2)


    @staticmethod
    def measure_dihedral(a: np.ndarray, b: np.ndarray, 
                         c: np.ndarray, d: np.ndarray) -> float:
        v1, v2 = Lipid.vectorize(a, b, d, c)
        return np.arccos(v1 @ v2)


    @staticmethod
    def measure_improper(a: np.ndarray, b: np.ndarray,
                         c: np.ndarray, d: np.ndarray) -> float:
        v1, v2, v3 = Lipid.vectorize(a, b, d, b, c, b)
        plane_normal = np.cross(v1, v2)
        return 90 - np.arccos(plane_normal @ v3)

