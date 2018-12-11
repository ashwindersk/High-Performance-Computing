stencil: stencil.c
	icc -std=c99 -qopenmp-stubs  -xHost -g -pg -Ofast -Wall $^ -o $@



