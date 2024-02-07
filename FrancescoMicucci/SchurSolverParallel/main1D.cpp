#include <chrono>
#include <iostream>
#include "heat1DSchur.hpp"

#include <fenv.h>

// Main function
int main(int argc, char *argv[])
{
#ifndef NDEBUG
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    MPI_Init(&argc, &argv);

    unsigned int spaces[] = {351, 3231, 32031, 320031, 3200031, 32000031, 320000031, 3200000031};

    constexpr size_t iterations = std::end(spaces) - std::begin(spaces);

    for (unsigned int idx = 0; idx < iterations; ++idx)
    {
        const unsigned int n = spaces[idx]; // Number of points
        // const int numDecomp = 4;               // Number of decomposition of the system
        double dx = 1.0 / (double)(n + 1);             // Step size
        //double *rhs = new double[n];           // Rhs of the problem
        //double *exactSolution = new double[n]; // Exact solution of the problem
        double diagElem = 2.0 / (dx * dx);     // Diagonal element of the matrix
        double upperElem = -1.0 / (dx * dx);   // Upper diagonal element of the matrix
        double lowerElem = -1.0 / (dx * dx);   // Lower diagonal element of the matrix

        auto start = std::chrono::high_resolution_clock::now();

        // Get MPI size
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // Get MPI rank
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        SchurSolver solver(MPI_COMM_WORLD, n, size, dx, diagElem, upperElem, lowerElem);
        solver.solve();

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

        if (rank == 0)
        {
            std::cout << "============================================" << std::endl;
            std::cout << "Error with " << n << " points: " << solver.getError() << std::endl;
            std::cout << "Time: " << (double)elapsed / 1e9 << " s " << std::endl;
            std::cout << "Seconds per point: " << ((double)elapsed / 1e9) / (double)n << std::endl;
            std::cout << "============================================" << std::endl;
        }
    }

    MPI_Finalize();

    return 0;
}