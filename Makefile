stencil: stencil.c
	icc -std=c99 -qopenmp-stubs -xHost mpiicc -g -pg -Ofast -Wall $^ -o $@ 



