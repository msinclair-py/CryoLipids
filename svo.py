#!/usr/bin/env python
import numpy as np

inp = 'misc/collision2.pdb'
lines = open(inp).readlines()[1:-1]
protein = np.array([[float(x) for x in line[31:54].strip().split()] 
                    for line in lines if 'POV' not in line])
lipid = np.array([[float(x) for x in line[31:54].strip().split()] 
                  for line in lines if 'POV' in line])

class Voxel:
    def __init__(self, center, size):
        self.center = center
        self.size = size
        self.children = None
        self.points = []

class SparseVoxelOctree:
    def __init__(self, points, depth=5, max_points_per_voxel=8):
        self.root = None
        self.depth = depth
        self.max_points_per_voxel = max_points_per_voxel
        self.build(points)

    def build(self, points):
        # Calculate the bounding box of the points
        min_coord = np.min(points, axis=0)
        max_coord = np.max(points, axis=0)
        bbox_center = (min_coord + max_coord) / 2
        bbox_size = max(max_coord - min_coord)

        # Create root voxel
        self.root = Voxel(bbox_center, bbox_size)

        # Recursively construct the octree
        self._subdivide(self.root, points, 0)

    def _subdivide(self, voxel, points, depth):
        if depth >= self.depth or len(points) <= self.max_points_per_voxel:
            voxel.points = points
            return

        # Subdivide voxel into 8 children
        voxel.children = []
        half_size = voxel.size / 2
        for i in range(8):
            child_center = voxel.center + half_size * np.array([
                (i & 1) * 2 - 1,
                ((i >> 1) & 1) * 2 - 1,
                ((i >> 2) & 1) * 2 - 1
            ])
            child_points = []
            for point in points:
                if np.all(np.abs(point - child_center) < half_size):
                    child_points.append(point)
            child_voxel = Voxel(child_center, half_size)
            self._subdivide(child_voxel, child_points, depth + 1)
            voxel.children.append(child_voxel)

    def query(self, query_point):
        return self._query_recursive(self.root, query_point)

    def _query_recursive(self, voxel, query_point):
        if voxel.children is None:
            return voxel.points
        
        points = []
        for child in voxel.children:
            if np.all(np.abs(query_point - child.center) < child.size / 2):
                points.extend(self._query_recursive(child, query_point))

        return points

octree = SparseVoxelOctree(protein, depth=2)
for pt in lipid:
    print(f'{octree.query(pt)=}')
