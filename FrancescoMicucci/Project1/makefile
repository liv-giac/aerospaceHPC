.PHONY: all
all: project1

project1: project1.cpp parallelHeat.cpp parallelHeat.hpp
	mpic++ project1.cpp parallelHeat.cpp parallelHeat.hpp -o project1 -Wall -O3 -g
	mpiexec -n 6 ./project1

.PHONY: clear
clear:
	rm -f project1