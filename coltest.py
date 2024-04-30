#!/usr/bin/env python
from collision_detection import CollisionDetector
import numpy as np

inp = 'misc/collision.pdb'
lines = open(inp).readlines()[1:-1]
protein = np.array([[float(x) for x in line[31:54].strip().split()] 
                    for line in lines if 'POV' not in line])
lipid = np.array([[float(x) for x in line[31:54].strip().split()] 
                  for line in lines if 'POV' in line])

linear_programming = CollisionDetector(protein, lipid)
delaunay_triangulation = CollisionDetector(protein, lipid, method=1)
#sparse_voxel_octree = CollisionDetector(protein, lipid, method=2)

print(f'{linear_programming.query_points()=}')
print(f'{delaunay_triangulation.query_points()=}')
#print(f'{sparse_voxel_octree.query_points()=}')
