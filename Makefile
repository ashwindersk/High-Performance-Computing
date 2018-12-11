stencil: stencil.c
	mpiicc -std=c99 -qopenmp-stubs  -xHost -g -pg -Ofast -Wall $^ -o $@



