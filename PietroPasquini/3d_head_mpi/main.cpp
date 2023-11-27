#include <fstream>
#include <cstdlib>
#include <chrono>
#include "parallel_heat.hpp"

/**
 * Checks if a given integer is a perfect square.
 *
 * @param n The integer to check.
 * @return True if n is a perfect square, false otherwise.
 */
bool isPerfectSquare(int n)
{
    int root = sqrt(n);
    return n == root * root;
}

int main(int argc, char **argv)
{
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (isPerfectSquare(mpi_size) == false && mpi_rank == 0)
    {
        std::cout << "This application is meant to be run with a perfect square number of processes." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    if (mpi_rank == 0)
    {
        std::cout << "========================================================" << std::endl;
        std::cout << "Running with " << mpi_size << " processes." << std::endl;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    // Initialize the parallel heat equation solver
    ParallelHeat parallel_heat(mpi_rank, mpi_size);

    parallel_heat.initialize_u();

    for (unsigned int t = 1; t <= parallel_heat.TSteps; t++)
    {
        parallel_heat.time += parallel_heat.dt;

        if (mpi_rank == 0)
        {
            std::cout << "========================================================" << std::endl;
            std::cout << "Time step: " << t << std::endl;
        }

        if (mpi_rank == 0)
            std::cout << "Message passing..." << std::endl;

        // Prepare auxiliary vectors
        parallel_heat.message_passing();

        if (mpi_rank == 0)
            std::cout << "Building rhs..." << std::endl;

        // Build the right-hand side of the equation
        parallel_heat.build_rhs();

        if (mpi_rank == 0)
            std::cout << "Applying boundary conditions..." << std::endl;

        // Apply boundary conditions
        parallel_heat.bcs_x();

        if (mpi_rank == 0)
            std::cout << "Solving the linear systems..." << std::endl;

        // Solve the system of equations
        parallel_heat.thomas();

        MPI_Barrier(MPI_COMM_WORLD);

        if (mpi_rank == 0)
            std::cout << "Changing direction x -> y" << std::endl;

        // Change direction x -> y
        parallel_heat.change_direction_x_to_y();

        if (mpi_rank == 0)
            std::cout << "Applying boundary conditions..." << std::endl;

        // Apply boundary conditions
        parallel_heat.bcs_y();

        if (mpi_rank == 0)
            std::cout << "Solving the linear systems..." << std::endl;

        // Solve the system of equations
        parallel_heat.thomas();

        MPI_Barrier(MPI_COMM_WORLD);

        if (mpi_rank == 0)
            std::cout << "Changing direction y -> z" << std::endl;

        // Change direction y -> z
        parallel_heat.change_direction_y_to_z();

        if (mpi_rank == 0)
            std::cout << "Applying boundary conditions..." << std::endl;

        // Apply boundary conditions
        parallel_heat.bcs_z();

        if (mpi_rank == 0)
            std::cout << "Solving the linear systems..." << std::endl;

        // Solve the system of equations
        parallel_heat.thomas();

        MPI_Barrier(MPI_COMM_WORLD);

        if (mpi_rank == 0)
            std::cout << "Changing direction z -> x" << std::endl;

        // Change direction z -> x
        parallel_heat.change_direction_z_to_x();

        if (mpi_rank == 0)
            std::cout << "Adding rhs to solution..." << std::endl;

        // Sum the right-hand side to the solution
        parallel_heat.sum_rhs_to_u();

        if (mpi_rank == 0)
            std::cout << "========================================================" << std::endl;
    }

    if (mpi_rank == 0){
        std::cout << "========================================================" << std::endl;
        std::cout << "Calculating error..." << std::endl;
    }

    // Print the error
    parallel_heat.print_error();

    MPI_Finalize();
    auto end_time = std::chrono::high_resolution_clock::now();

    if (mpi_rank == 0)
    {
        std::cout << "========================================================" << std::endl;
        std::cout << "Time per element per time step: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() / 1.0E11 << " s" << std::endl;
        std::cout << "========================================================" << std::endl;
    }

    return EXIT_SUCCESS;
}