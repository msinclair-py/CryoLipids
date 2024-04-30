#!/usr/bin/env python
import numpy as np

inp = 'misc/collision.pdb'
lines = open(inp).readlines()[1:-1]
protein = np.array([[float(x) for x in line[31:54].strip().split()]
                    for line in lines if 'POV' not in line])
lipid = np.array([[float(x) for x in line[31:54].strip().split()]
                  for line in lines if 'POV' in line])
llines = np.array([line for line in lines if 'POV' in line])

class OccupancyMap:
    def __init__(self, points, grid_spacing=1, pad=1):
        self.points = points
        self.grid_spacing = grid_spacing
        self.pad = pad
        self.grid = self.fill_grid(points)
        
    def define_bounds(self):
        self.minx = np.floor(np.min(self.points[:, 0])) - self.pad
        self.miny = np.floor(np.min(self.points[:, 1])) - self.pad
        self.minz = np.floor(np.min(self.points[:, 2])) - self.pad
        
        lenx = np.ceil(np.max(self.points[:, 0])) + self.pad - self.minx
        leny = np.ceil(np.max(self.points[:, 1])) + self.pad - self.miny
        lenz = np.ceil(np.max(self.points[:, 2])) + self.pad - self.minz

        self.dims = [int(i / self.grid_spacing) for i in [lenx, leny, lenz]]

    def generate_grid(self):
        return np.zeros((self.dims))

    def fill_grid(self, points):
        self.define_bounds()
        grid = self.generate_grid()
        for point in points:
            xi, yi, zi = self.get_indices(point)
            grid[xi, yi, zi] = 1

        return grid

    def get_indices(self, point):
        xi = np.floor(point[0]) - self.minx
        yi = np.floor(point[1]) - self.miny
        zi = np.floor(point[2]) - self.minz
        return (int(xi / self.grid_spacing), 
                int(yi / self.grid_spacing), 
                int(zi / self.grid_spacing))

    def query(self, points):
        if not isinstance(points, np.ndarray):
            points = np.array(points, ndmin=2)

        clashes = []
        for i, point in enumerate(points):
            xi, yi, zi = self.get_indices(point)
            if self.grid[xi, yi, zi]:
                clashes.append(i)

        return clashes

occ = OccupancyMap(protein, grid_spacing=1)
print(occ.query(lipid))
print(len(np.where(occ.grid == 0)[0]))
print(llines[occ.query(lipid)])
