.PHONY: all
all: parallel_heat

parallel_heat: main.cpp parallel_heat.cpp parallel_heat.hpp
	mpic++ main.cpp parallel_heat.cpp -O3 -o parallel_heat
	mpirun -np 4 ./parallel_heat