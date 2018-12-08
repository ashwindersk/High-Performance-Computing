stencil: stencil.c
	icc -std=c99 -qopenmp-stubs epic -xHost mpiicc -g -pg -Ofast -Wall $^ -o $@ 



