stencil: stencil.c
	icc -std=c99 -qopenmp-stubs epic -xHost mpicc -g -pg -Ofast -Wall $^ -o $@ 



