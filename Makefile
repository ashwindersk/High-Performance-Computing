stencil: stencil.c
	icc -std=c99 -qopenmp-stubs epic -xHost -g -pg -Ofast -Wall $^ -o $@ 



