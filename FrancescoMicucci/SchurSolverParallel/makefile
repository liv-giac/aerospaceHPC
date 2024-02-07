.PHONY: all
all: main1D # main3D

main1D: main1D.cpp heat1DSchur.cpp heat1DSchur.hpp
	mpic++ main1D.cpp heat1DSchur.cpp heat1DSchur.hpp -o main1D -Wall -O3 -march=native
	mpiexec -n 32 ./main1D 

# main3D: main3D.cpp heat3DSchur.cpp heat3DSchur.hpp
# 	g++ main3D.cpp heat3DSchur.cpp heat3DSchur.hpp -o main3D -Wall -O3 -g
# 	./main3D

.PHONY: clear
clear:
	rm -f main1D
	rm -f main3D