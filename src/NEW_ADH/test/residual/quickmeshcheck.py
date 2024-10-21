#quick mesh check
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


#came from mesh builder
nds = np.array([[0.000000,0.000000],
    [0.000000,1.000000],
    [0.000000,2.000000],
    [0.000000,3.000000],
    [1.000000,0.000000],
    [1.000000,1.000000],
    [1.000000,2.000000],
    [1.000000,3.000000],
    [2.000000,0.000000],
    [2.000000,1.000000],
    [2.000000,2.000000],
    [2.000000,3.000000],
    [3.000000,0.000000],
    [3.000000,1.000000],
    [3.000000,2.000000],
    [3.000000,3.000000]])
connectivity = np.array([
[0,4,5],
[0,5,1],
[1,5,2],
[5,6,2],
[2,6,7],
[2,7,3],
[4,8,9],
[4,9,5],
[5,9,6],
[9,10,6],
[6,10,11],
[6,11,7],
[8,12,13],
[8,13,9],
[9,13,10],
[13,14,10],
[10,14,15],
[10,15,11]],dtype=np.int32)

# Create the triangulation object
triang = tri.Triangulation(nds[:,0], nds[:,1], connectivity)

# Plot the mesh
plt.triplot(triang, 'k-') 
plt.savefig("Rectangular_mesh.png")