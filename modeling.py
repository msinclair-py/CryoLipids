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


def calc_dihedral(u1, u2, u3, u4):
    """ Calculate dihedral angle method. From bioPython.PDB
    (adapted to np.array)
    Calculate the dihedral angle between 4 vectors
    representing 4 connected points. The angle is in
    [-pi, pi].
    """

    a1 = u2 - u1
    a2 = u3 - u2
    a3 = u4 - u3

    v1 = np.cross(a1, a2)
    v1 = v1 / (v1 * v1).sum(-1)**0.5
    v2 = np.cross(a2, a3)
    v2 = v2 / (v2 * v2).sum(-1)**0.5
    porm = np.sign((v1 * a3).sum(-1))
    rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
    if not porm == 0:
        rad = rad * porm

    return rad

a = np.array([1219.7, 4106.1, -7366.7])
b = np.array([1277.1, 4016.6, -7447.1])
c = np.array([1398.6, 3944.8, -7407.8])
d = np.array([1501.2, 3943.2, -7521.3])

alpha = calc_dihedral(a,b,c,d)
