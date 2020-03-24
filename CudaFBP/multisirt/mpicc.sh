mpicc -c main-fbp.c -o main-fbp.o
mpicc -c atom-fbp.c -o atom-fbp.o 
mpicc -c mrcfile_atom2.c -o mrcfile_atom.o

mpicc main-fbp.o atom-fbp.o mrcfile_atom.o -o main -lm  
rm *.o

