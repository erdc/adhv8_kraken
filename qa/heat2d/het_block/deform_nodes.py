#! /usr/bin/env python

"""
quick script to map nodes of mesh to paraboloid
z = h_0 - \left(\frac{2(x-x_c)}{L_x}\right)^2 -\left(\frac{2(y-y_c)}{L_y}\right)^2   
"""
import os
import numpy as np

def para_z(x,y,z,h_0=10.0,Lx=50.,Ly=50.,x_c=25.0,y_c=-25.,c_x=10.,c_y=10.):
    return h_0 - (c_x*(x-x_c)/Lx)**2 - (c_y*(y-y_c)/Ly)**2
def map_mesh_file(meshbase,zfunc,ext='3dm'):
    assert os.path.isfile(meshbase+'.'+ext)
    fin = open(meshbase+'.'+ext,'r')
    fout=open(meshbase+'paraboloid'+'.'+ext,'w')
    
    line = fin.readline()
    while line:
        if line[0:2] == 'ND':
            terms = line.split()
            assert len(terms) == 5

            x,y,z=float(terms[2]),float(terms[3]),float(terms[4])
            
            zz=zfunc(x,y,z)
            new_line='{0} {1} {2} {3} {4:12.5e}\n'.format(terms[0],terms[1],
                                                   terms[2],terms[3],
                                                   zz)
            fout.write(new_line)
        else:
            fout.write(line)
        #
        line = fin.readline()
    #
    fin.close()
    fout.close()
    
    
            
