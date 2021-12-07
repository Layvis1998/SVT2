prog: flow.cpp counterflow.cpp
	g++  flow.cpp lib/libblas.so lib/libumfpack.a lib/libamd.a -O0 -g -o flow
	g++  counterflow.cpp lib/libblas.so lib/libumfpack.a lib/libamd.a -O0 -g -o counterflow
clean:
	rm prog
