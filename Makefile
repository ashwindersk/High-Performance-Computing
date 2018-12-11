stencil: stencil.c
	mpicc -std=c99  -g -pg -Ofast -Wall $^ -o $@ 



