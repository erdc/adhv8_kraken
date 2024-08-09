import h5py as h5
import sys
import numpy as np


f1 = h5.File(sys.argv[1],'w')  

grp = f1.create_group('NodalAttributes', track_order=True)
dset = grp.create_dataset("Rain", (21), dtype='double')


dat = np.linspace(0,20,21,dtype=np.double)
print(dat)

dset[:] = dat
f1.close()

f1 = h5.File(sys.argv[1],'r') 
print(f1["NodalAttributes"]["Rain"][:].dtype)
#for a in f1.keys():
#    print(f1[a].keys())

