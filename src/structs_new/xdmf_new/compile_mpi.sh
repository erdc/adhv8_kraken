export LIBRARY_PATH=/usr/local/opt/hdf5_parallel/lib/
mpicc standalone_mesh_mpi.c -o test_mpi -lhdf5 -I /usr/local/opt/hdf5_parallel/include/ -Wno-implicit-function-declaration
