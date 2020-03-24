mpicc -c main-fbp.c -o main-fbp.o
mpicc -std=c99  -I /home/x5wan/software/OpenCV-2.2.0/include -c atom-fbp.c -o atom-fbp.o 
mpicc -c mrcfile_atom.c -o mrcfile_atom.o
mpicc -I /home/x5wan/software/OpenCV-2.2.0/include -c filter-prj3.c -o filter-prj.o
mpicc main-fbp.o atom-fbp.o mrcfile_atom.o filter-prj.o -o main -lm  -L /home/x5wan/software/OpenCV-2.2.0/lib -lopencv_core -lopencv_highgui -lopencv_imgproc
rm *.o

