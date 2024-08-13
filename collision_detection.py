import numpy as np
from scipy.optimize import linprog
from scipy.spatial import Delaunay
from scipy.spatial.transform import Rotation
from typing import List, Tuple
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
    NOTE: Incorporate networkx lipid graphs
    
    Class to orchestrate collision detection and repair. Takes in an instance of
    the `Lipid` class for both its coordinates and a few methods, protein coordinates
    as a numpy array, and an instance of CollisionDetector to perform the actual
    collision detection.
    """
    def __init__(self, lipid: Lipid, 
                 protein: np.ndarray, 
                 collision_detector: int=0,
                 **kwargs):
        self.lipid = lipid
        self.detector = CollisionDetector(protein, lipid, method=collision_detector, **kwargs)
    
    def check_collisions(self) -> None:
        """
        Main logic of class controlling the flow of collision detection and
        subsequent repair.
        """
        clash = True
        while clash:
            clashes = self.detector.query_points()
            if not clashes:
                break
            
            self.repair_tail_clashes(clashes)

    def repair_tail_clashes(self, clashes: List[str]) -> None:
        """_summary_

        Args:
            clashes (List[str]): _description_

        Raises:
            NotImplementedError: _description_
        """
        c2s, c3s = [], []

        
        for clash in clashes:
            if 'C2' in clash:
                c2s.append(clash)
            elif 'C3' in clash:
                c3s.append(clash)
            else:
                raise NotImplementedError('Clash on non-tail detected!')
        
        for tail, tail_length in zip([c2s, c3s], [18, 16]):
            if tail:
                atoms_to_rotate = self.get_clash_rotation(tail, tail_length)
                old_coords = self.lipid.get_coord(atoms_to_rotate)
                new_coords = self.rotate_tail(old_coords)
                self.lipid.update_coords(atoms_to_rotate[2:], new_coords)
                self.detector.lipid_coordinates = self.lipid.extract_coordinates()
                

    @staticmethod
    def get_clash_rotation(clashing_atoms: List[str], length: int) -> List[str]:
        """
        Determines what atoms comprise the bond we need to rotate about, and the
        rest of the tail atoms which need to be rotated.
        Args:
            clashing_atoms (List[str]): List of atom names which are clashing with protein
            length (int): Length of lipid tail (18 for SN-1, 16 for SN-2)

        Returns:
            List[str]: List whose first two elements correspond to the bond to be rotated,
                and whose remaining elements are to be actually rotated in cartesian space.
        """
        tail_type = clashing_atoms[0][:2]
        first_clash = min([int(name[2:]) for name in clashing_atoms])
 
        atoms_to_rotate = [f'{tail_type}{i}' for i in range(first_clash - 2, length + 1)]
        
        return atoms_to_rotate
    
    @staticmethod
    def rotate_tail(tail_atoms: np.ndarray, rotate_by: float=15.) -> np.ndarray:
        """
        Performs a rotation about the bond between the first two atoms of `tail_atoms`.
        Units of rotation are in degrees.
        Args:
            tail_atoms (np.ndarray): Numpy array of coordinates. First two particles comprise
                the arbitrary axis to be rotated about.
            rotate_by (float): Number of degrees to rotate about arbitrary axis.
        Returns:
            np.ndarray: Array of transformed coordinates. Shape is (N-2, 3) where N is the 
                number of rows in the input array `tail_atoms`.
        """
        a1, a2, *to_rotate = tail_atoms # unpack into first two atoms and rest of atoms
        vector = a2 - a1
        
        align, _ = Rotation.align_vectors(np.array([0, 0, 1]), vector)
        rotate = Rotation.from_euler('z', rotate_by, degrees=True)
        put_back = align.inv()
        
        tail_atoms = align.apply(tail_atoms - a1)
        tail_atoms = rotate.apply(tail_atoms)
        tail_atoms = put_back.apply(tail_atoms) + a1
        
        return put_back.apply(rotate.apply(align.apply(to_rotate - a1))) + a1

    def get_new_coords(self):
        return self.lipid.pdb_contents