stencil: stencil.c
	mpiicc -std=c99  -g -pg -Ofast -Wall $^ -o $@ 



