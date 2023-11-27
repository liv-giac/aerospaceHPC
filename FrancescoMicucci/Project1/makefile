.PHONY: all
all: project1

project1: project1.cpp parallelHeat.cpp parallelHeat.hpp
	mpic++ project1.cpp parallelHeat.cpp parallelHeat.hpp -o project1 -Wall -O3 -g
	mpiexec -n 1 ./project1
	mpiexec -n 2 ./project1
	mpiexec -n 3 ./project1
	mpiexec -n 4 ./project1
	mpiexec -n 5 ./project1
	mpiexec -n 6 ./project1

.PHONY: clear
clear:
	rm -f project1