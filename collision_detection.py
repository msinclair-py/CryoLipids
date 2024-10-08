import networkx as nx
import numpy as np
from scipy.optimize import linprog
from scipy.spatial import Delaunay
from scipy.spatial.transform import Rotation
from typing import Dict, List, Tuple
from modeling import Lipid

class CollisionDetector:
    """
    Perform collision detection of two point clouds. Available methods for
    computing detection include: linear programming, delaunay triangulation and 
    sparse voxel octree algorithms. These are more well described in the 
    corresponding docstrings for each class below.
    """
    def __init__(self, protein: np.ndarray, lipid: Lipid, 
                 method: int=0, **kwargs):
        if not isinstance(protein, np.ndarray):
            protein = np.array(protein)
        
        self.raw_lipid = np.array(lipid.pdb_contents)
        self.lipid_coordinates = lipid.extract_coordinates()
        methods = [OccupancyMap, LinearProgramming, DelaunayTriangulation, SVO]
        self.detector = methods[method](protein, **kwargs)
        
    def query_points(self, point: None=None) -> List[bool]:
        if point is None:
            clashes = self.detector.query(self.lipid_coordinates)
        else:
            clashes = self.detector.query(point)

        if any(clashes):
            return [line[2] for line in self.raw_lipid[clashes]]
        else:
            return False
        
class OccupancyMap:
    """
    Grids the protein space into a grid of shape (ceil(max) - floor(min) + 2*pad) 
    in each dimension. This allows for a clean, even gridspace to map atomic occupancies.
    Each voxel of the grid is assigned 1 if there are any protein atoms and 0 if there
    are no protein atoms. The grid can then be queried to see if any lipid atoms are in
    an occupied voxel, meaning a clash.
    ------
    Inputs
    ------
    points (np.ndarray): coordinates of protein atoms
    grid_spacing (float): spacing of grid in Angstroms
    pad (int): padding to apply in Angstroms
    """
    def __init__(self, points: np.ndarray, grid_spacing: float=1.5,
                 pad: int=1):
        self.points = points
        self.grid_spacing = grid_spacing
        self.pad = pad
        self.grid = self.fill_grid(points)
        
    def define_axis(self, points: np.ndarray) -> np.ndarray:
        """
        Defines a single axis in the grid space. Performs the padding as well
        as regularizing the grid space such that it is neatly divisible by
        the `self.grid_spacing` parameter.
        """
        min_ = np.floor(np.min(points)) - self.pad
        max_ = np.ceil(np.max(points)) + self.pad
        
        len_ = max_ - min_ + self.grid_spacing
        regularized_len = len_ - len_ % self.grid_spacing
        return min_, regularized_len
    
    def define_bounds(self) -> None:
        """
        Defines the full grid parameters in x, y, z.
        """
        minx, lenx = self.define_axis(self.points[:, 0])
        miny, leny = self.define_axis(self.points[:, 1])
        minz, lenz = self.define_axis(self.points[:, 2])

        self.mins = [minx, miny, minz]
        self.dims = [int(i / self.grid_spacing) for i in [lenx, leny, lenz]]
        
    def generate_grid(self) -> np.ndarray:
        """
        Initialize empty grid.
        """
        return np.zeros((self.dims))
    
    def fill_grid(self, points: np.ndarray) -> np.ndarray:
        """
        Fill grid with occupancies based on `points`.
        """
        self.define_bounds()
        grid = self.generate_grid()
        for point in points:
            xi, yi, zi = self.get_indices(point)
            grid[xi, yi, zi] = 1
        
        return grid
    
    def get_indices(self, point: np.ndarray) -> Tuple[float]:
        """
        Gets index corresponding to correct voxel in grid. Takes the
        `self.grid_spacing` into account.
        """
        return [int(np.floor(point[i] - self.mins[i]) / self.grid_spacing) 
                for i in range(len(self.mins))]
    
    def query(self, points: np.ndarray) -> List[int]:
        """
        Queries whether a given point or points are in any occupied voxels.
        """
        if not isinstance(points, np.ndarray):
            points = np.array(points, ndmin=2)
            
        clashes = []
        for i, point in enumerate(points):
            xi, yi, zi = self.get_indices(point)
            try:
                if self.grid[xi, yi, zi]:
                    clashes.append(i)
            except IndexError:
                continue
            
        return clashes
            
class LinearProgramming:
    """
    Perform collision detection using linear programming. This avoids the need
    for constructing an actual object such as in DelaunayTriangulation or 
    SparseVoxelOctree.
    """
    def __init__(self, points: np.ndarray):
        self.points = points
        self.n_points = len(points)
        
    def query(self, x: np.ndarray) -> bool:
        c = np.zeros(self.n_points)
        A = np.r_[self.points.T, np.ones((1, self.n_points))]
        b = np.r_[x, np.ones(1)]
        lp = linprog(c, A_eq=A, b_eq=b)
        return lp.success
    
    
class DelaunayTriangulation:
    """
    Compute the Delaunay triangulation of our reference set of points. Then
    query whether or not a new point is within our triangulation or not.
    """
    def __init__(self, points: np.ndarray):
        self.hull = Delaunay(points)
        
    def query(self, point: np.ndarray) -> bool:
        return self.hull.find_simplex(point) >= 0

class SparseVoxelOctree:
    """
    Recursively builds a Sparse Voxel Octree from a set of points.
    """
    def __init__(self, root, depth):
        self.root = root
        self.depth = depth
    
    def insert(self, point):
        node = self.root
        for i in range(self.depth):
            bit = (point[i] >> (self.depth - 1 - i)) & 1
            if node.children[bit] is None:
                node.children[bit] = SparseVoxelOctreeNode()
            node = node.children[bit]
        
    def query(self, point):
        node = self.root
        for i in range(self.depth):
            bit = (point[i] >> (self.depth - 1 - i)) & 1
            if node.children[bit] is None:
                return False
            node = node.children[bit]
        return True
  
    
class SparseVoxelOctreeNode:
    """
    Nodes in Sparse Voxel Octree.
    """
    def __init__(self):
        self.children = [None] * 8
  
        
class SVO:
    """
    Constructs Sparse Voxel Octree from point cloud data. Query method allows
    us to perform collision detection by identify whether a point is within
    one of the sparsely populated voxels. Preserves concave features unlike
    Delaunay methods, but suffers from the tradeoff of high resolution vs
    computational efficiency.
    """
    def __init__(self, points: np.ndarray, depth: int = 8):
        self.octree = SparseVoxelOctree(SparseVoxelOctreeNode(), depth)
        for point in points:
            self.octree.insert(point)
            
    def query(self, point: np.ndarray) -> bool:
        return self.octree.query(point)

class Repairer:
    """
    Class to orchestrate collision detection and repair. Takes in an instance of
    the `Lipid` class for both its coordinates and a few methods, protein coordinates
    as a numpy array, and an instance of CollisionDetector to perform the actual
    collision detection.
    """
    def __init__(self, lipid: Lipid, 
                 protein: np.ndarray, 
                 theta: float=15.,
                 clash_tolerance: int=24,
                 collision_detector: int=0,
                 **kwargs):
        self.lipid = lipid
        self.graph = lipid.graph
        self.top_sort = {val: i for i, val in enumerate(nx.topological_sort(lipid.graph))}
        self.detector = CollisionDetector(protein, lipid, method=collision_detector, **kwargs)
        self.theta = theta
        self.clash_tolerance = max(int(360. / theta), clash_tolerance)
    
    def check_collisions(self) -> None:
        """
        Main logic of class controlling the flow of collision detection and
        subsequent repair.

        Args:
            theta (float): Number of degrees to rotate clashes for repair.
            clash_tolerance (int): Min. number of clash repair attempts before
                            exiting due to lack of convergence.

        Returns:
            [int]: Success value, 1 means successful repair, 0 means repair has
                not converged indicating energy minimization is strictly required 
                to fix protein-lipid system.
        """

        clash_counter = 0
        is_repaired = True
        while True:
            clashes = self.detector.query_points()
            if not clashes:
                break
            
            clash_counter += 1

            if clash_counter > self.clash_tolerance:
                is_repaired = False
                break

            self.repair_tail_clashes(clashes)

        if is_repaired:
            return 1

        return 0

    def repair_tail_clashes(self, clashes: List[str]) -> None:
        """
        Repair clash which appears highest in the topological representation
        of the lipid (C2 is the top atom in the graph). While this may or may
        not fix clashes which appear lower in the topological heirarchy it is
        expected that subsequent collision repair will fix these.

        Args:
            clashes (List[str]): Atom names corresponding to atomic clashes
        """
        if len(clashes) == 1:
            clashes_to_resolve = clashes
        else:
            clashes_to_resolve = []
            for clash in clashes:
                ancestors = list(nx.ancestors(self.graph, clash))
                if not any([clash in ancestors for clash in clashes]):
                    clashes_to_resolve.append(clash)
        
        for clash in clashes_to_resolve:
            atoms_to_rotate = self.get_clash_rotation(clash)
            old_coords = self.lipid.get_coord(atoms_to_rotate)
            new_coords = self.rotate_tail(old_coords, self.theta)
            self.lipid.update_coords(atoms_to_rotate[2:], new_coords)
            self.detector.lipid_coordinates = self.lipid.extract_coordinates()
                
    def get_clash_rotation(self, clashing_atom: str) -> List[str]:
        """
        Determines what atoms comprise the bond we need to rotate about, and the
        rest of the remaining atoms which need to be rotated.
        Args:
            clashing_atoms (List[str]): List of atom names which are clashing with protein

        Returns:
            List[str]: List whose first two elements correspond to the bond to be rotated,
                and whose remaining elements are to be actually rotated in cartesian space.
        """            
        prev1 = list(self.graph.predecessors(clashing_atom))[0]
        prev2 = list(self.graph.predecessors(prev1))[0]
        
        return [prev2, prev1] + list(nx.descendants(self.graph, prev1))
    
    @staticmethod
    def rotate_tail(arr: np.ndarray, theta: float=15.) -> np.ndarray:
        """
        Rotates tail about arbitrary axis defined by first two atoms.

        Args:
            tail (np.ndarray): Array of atoms to rotate
            theta (float, optional): Angle to rotate about in degrees, 
                                        defaults to 15.0

        Returns:
            np.ndarray: Rotated atoms
        """
        center = arr[0]
        vector = arr[1] - center
        ux, uy, uz = vector / np.linalg.norm(vector)
        
        cos = np.cos(theta)
        sin = np.sin(theta)
        
        rot_matrix = np.array([
            [cos + ux**2*(1-cos), ux*uy*(1-cos) - uz*sin, ux*uz*(1-cos) + uy*sin],
            [uy*ux*(1-cos), cos + uy**2*(1-cos), uy*uz*(1-cos) - ux*sin],
            [uz*ux*(1-cos) - uy*sin, uz*uy*(1-cos) + ux*sin, cos + uz**2*(1-cos)]
            ])
        
        R = Rotation.from_matrix(rot_matrix)
        return R.apply(arr[2:] - center) + center

    def get_new_coords(self) -> List[str]:
        """
        Obtains both the new coordinates for `self.lipid`, and also a binary
        assignment for atoms to fix in energy minimization/implicit equilibration.

        Returns:
            np.ndarray: Array of shape (N, 4) where the first 3 columns are x, y, z and
                        the last column is a 1 or 0 which indicates it should be fixed
                        or not respectively.
        """
        pdb_lines = self.lipid.pdb_contents
        fixed = self.get_fixed(pdb_lines)
        for pdb_line, value in zip(pdb_lines, fixed):
            pdb_line[-1] = value
            
        return pdb_lines
    
    def get_fixed(self, lines: List[str]) -> List[int]:
        """
        Identify which atoms should be fixed by comparing to the original PDB.
        
        Args:
            coords (List[str]): The new lines of what will ultimately be our
                                output PDB.

        Returns:
            np.ndarray: A vector corresponding to which atom indices will be fixed.
        """
        original_atoms = [line[2] for line in self.lipid.lipid_lines]
        return [1 if line[2] in original_atoms else 0 for line in lines]

