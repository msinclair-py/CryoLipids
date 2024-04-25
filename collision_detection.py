import numpy as np
from scipy.optimize import linprog
from scipy.spatial import Delaunay
from typing import List

class CollisionDetector:
    """
    Perform collision detection of two point clouds. Available methods for
    computing detection include: linear programming, delaunay triangulation and 
    sparse voxel octree algorithms. These are more well described in the 
    corresponding docstrings for each class below.
    """
    def __init__(self, protein: np.ndarray, lipid: np.ndarray, 
                 method: int=0, **kwargs):
        self.lipid = lipid
        methods = [LinearProgramming, DelaunayTriangulation, SVO]
        self.detector = methods[method](protein, **kwargs)
        
    def query_points(self, point: None) -> List[bool]:
        if point is not None:
            self.detector.query(point)
        else:
            for point in self.lipid:
                self.detector.query(point)
    
    
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