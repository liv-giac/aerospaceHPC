#include <chrono>
#include <iostream>
#include "heat1DSchur.hpp"

#include <fenv.h>

// Initialize rhs and exact solution of the pb.
void initialize(double rhs[], double exactSolution[], const int &n, const double &dx)
{
    for (int i = 0; i < n - 1; i++)
    {
        exactSolution[i] = std::sin((i + 1) * dx);
        rhs[i] = std::sin((i + 1) * dx);
    }
    exactSolution[n - 1] = std::sin(n * dx);
    rhs[n - 1] = std::sin(1.0) / (dx * dx) + std::sin(n * dx);
}

// Main function
int main(int argc, char *argv[])
{
#ifndef NDEBUG
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    MPI_Init(&argc, &argv);

    const int n = 27; // Number of points
    // const int numDecomp = 4;               // Number of decomposition of the system
    double dx = 1.0 / (n + 1);             // Step size
    double *rhs = new double[n];           // Rhs of the problem
    double *exactSolution = new double[n]; // Exact solution of the problem
    double diagElem = 2.0 / (dx * dx);     // Diagonal element of the matrix
    double upperElem = -1.0 / (dx * dx);   // Upper diagonal element of the matrix
    double lowerElem = -1.0 / (dx * dx);   // Lower diagonal element of the matrix

    auto start = std::chrono::high_resolution_clock::now();

    // Get MPI rank
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    initialize(rhs, exactSolution, n, dx);
    SchurSolver solver(MPI_COMM_WORLD, n, size, dx, rhs, exactSolution, diagElem, upperElem, lowerElem);
    solver.solve();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;

    // Get rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
        std::cout << "============================================" << std::endl;
        std::cout << "Error with " << n << " points: " << solver.getError() << std::endl;
        std::cout << "Time: " << time.count() * 1000 << " milliseconds " << std::endl;
        std::cout << "Seconds per point: " << time.count() / n << std::endl;
        std::cout << "============================================" << std::endl;
    }

    MPI_Finalize();

    return 0;
}