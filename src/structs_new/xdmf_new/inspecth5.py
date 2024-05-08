import h5py as h5
import sys


f1 = h5.File(sys.argv[1],'r')  

print(f1.keys())


for a in f1.keys():
    print(f1[a][:])