import h5py as h5
import sys


f1 = h5.File(sys.argv[1],'r')

print(f1.keys())


print(f1['Mesh'].keys())

print("Nodal PEs")
print(f1['Data']['NodalScalar']['0'][:])

print("Elemental PEs")
print(f1['Data']['ElementalScalar']['0'][:])
#for a in f1.keys():
#    print(f1[a].keys())
