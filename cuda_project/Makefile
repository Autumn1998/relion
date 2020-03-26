objects=main.o read_write_mrc.o
MAIN = main.cu
main:$(objects)
	nvcc -o main $(objects)
main.o:$(MAIN) read_write_mrc.h atom.h
	nvcc -c $(MAIN)
read_write_mrc.o:read_write_mrc.cpp read_write_mrc.h
	g++ -c read_write_mrc.cpp
.PHONY:clean
clean:
	rm main $(objects)
