#!/usr/bin/python
import sys, getopt
import meshio
import numpy as np

if __name__ == "__main__":

    points = np.empty((0,3), float)
    cells = []

    filename = sys.argv[-1] + ".3dm"

    print(filename)

    with open(filename, 'r') as infile:
        for line in infile:
            data = line.split()
            if data[0] == "ND":
                points = np.append(points,[[float(data[2]), float(data[3]), float(data[4])]], axis=0)
                #print(points);

            if data[0] == "TRI":
                nodes = np.array([[int(data[2])-1, int(data[3])-1, int(data[4])-1]])
                geo = ("triangle", nodes)
                #cells = np.append(cells,geo)
                cells.append(geo)
                #print(cells)

            if data[0] == "TET":
                nodes = np.array([[int(data[2])-1, int(data[3])-1, int(data[4])-1, int(data[5])-1]])
                #nodes = np.array([[data[2], data[3], data[4], data[5]]])
                geo = ("tetra", nodes)
                #cells = np.append(cells,geo)
                cells.append(geo)
                #print(cells)

    meshio.write_points_cells(
        sys.argv[-1] + ".vtk",
        points,
        cells,
        # Optionally provide extra data on points, cells, etc.
        # point_data=point_data,
        #point_data cell_data=cell_data,
        # field_data=field_data
        )

    print("DONE")
