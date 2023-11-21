#include "parallel_heat.hpp"

#include <vector>

void ParallelHeat::printRhs()
{

    std::vector<double> rhs_gather;
    rhs_gather.resize(N * N * N);

    MPI_Gather(rhs, N * local_size * local_size, MPI_DOUBLE, rhs_gather.data(), N * local_size * local_size, MPI_DOUBLE, 0, comm2D);

    // Print rhs_local for every processor
    if (mpi_rank == 0)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    std::cout << rhs_gather[i * N * N + j * N + k] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
}


/**
 * @brief Calculates the index of a 3D array element given its 3D coordinates and the size of one dimension.
 *
 * @param i The index of the element along the first dimension.
 * @param j The index of the element along the second dimension.
 * @param k The index of the element along the third dimension.
 * @return The index of the element in the 1D array representation of the 3D array.
 */
inline unsigned int ParallelHeat::index(unsigned int i, unsigned int j, unsigned int k)
{
    return i + j * N + k * N * local_size;
}

/**
 * Calculates the force at a given point in space and time.
 *
 * @param x The x-coordinate of the point.
 * @param y The y-coordinate of the point.
 * @param z The z-coordinate of the point.
 * @param t The time at which to calculate the force.
 * @return The force at the given point and time.
 */
inline double ParallelHeat::force(double x, double y, double z, double t)
{
    using namespace std;
    return (cos(t) + 3.0 * sin(t)) * sin(x) * sin(y) * sin(z);
}

/**
 * Calculates the exact solution for a given set of input parameters.
 *
 * @param x The x-coordinate of the point to evaluate.
 * @param y The y-coordinate of the point to evaluate.
 * @param z The z-coordinate of the point to evaluate.
 * @param t The time at which to evaluate the solution.
 *
 * @return The exact solution for the given input parameters.
 */
inline double ParallelHeat::exact(double x, double y, double z, double t)
{
    using namespace std;
    return sin(x) * sin(y) * sin(z) * sin(t);
}

void ParallelHeat::initialize_u()
{
    double z = z_block_start_idx * dx;
    double y = y_block_start_idx * dx;
    double x = 0.0;

    for (unsigned int k = 0; k < local_size; k++)
    {
        z += dx;
        for (unsigned int j = 0; j < local_size; j++)
        {
            y += dx;
            for (unsigned int i = 0; i < N; i++)
            {
                x += dx;
                u[index(i, j, k)] = exact(x, y, z, 0.0);
            }
            x = 0.0;
        }
        y = y_block_start_idx * dx;
    }
}

void ParallelHeat::build_rhs()
{

    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                rhs[index(i, j, k)] = dt * (force((i + 1)*dx,(coord[0]*local_size + j + 1) * dx, (coord[1]*local_size + k + 1) * dx, time + 0.5 * dt) -
                                            (1.0 / dx2) * 6.0 * u[index(i, j, k)]);
                if (k > 0)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i, j, k - 1)];
                else if (coord[1] > 0)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * upper[index(i, j, 0)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * exact((i + 1) * dx, (coord[0]*local_size + j + 1) * dx, 0.0, time);

                if (k < local_size - 1)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i, j, k + 1)];
                else if(coord[1] < dims[1] - 1)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * lower[index(i, j, 0)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * exact((i + 1) * dx, (coord[0]*local_size + j + 1) * dx, 1.0, time);

                if (j > 0)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i, j - 1, k)];
                else if (coord[0] > 0)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * left[index(i, k, 0)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * exact((i + 1) * dx, 0.0, (coord[1]*local_size + k + 1) * dx, time);
                    
                if (j < local_size - 1)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i, j + 1, k)];
                else if(coord[0] < dims[0] - 1)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * right[index(i, k, 0)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * exact((i + 1) * dx, 1.0, (coord[1]*local_size + k + 1) * dx, time);

                if(i > 0)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i - 1, j, k)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * exact(0.0, (coord[0]*local_size + j + 1) * dx, (coord[1]*local_size + k + 1) * dx, time);
                if(i < N - 1)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i + 1, j, k)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * exact(1.0, (coord[0]*local_size + j + 1) * dx, (coord[1]*local_size + k + 1) * dx, time);
            }
        }
    }
}

void ParallelHeat::exchange_rhs_y()
{
    int source, process_right, process_left;
    int direction = 0;
    MPI_Datatype vertical_slice_type;
    MPI_Type_vector(local_size, N, N * local_size, MPI_DOUBLE, &vertical_slice_type);
    MPI_Type_commit(&vertical_slice_type);
    MPI_Cart_shift(comm2D, direction, 1, &source, &process_right);
    MPI_Cart_shift(comm2D, direction, -1, &source, &process_left);
    if (process_left >= 0)
    {
        // std::cout << "Process " << mpi_rank << " sending to " << process_left << std::endl;
        MPI_Send(u, 1, vertical_slice_type, process_left, mpi_rank, comm2D);
    }
    if (process_right >= 0)
    {
        // std::cout << "Process " << mpi_rank << " receiving from " << process_right << std::endl;
        MPI_Recv(right, N*local_size, MPI_DOUBLE, process_right, process_right, comm2D, MPI_STATUS_IGNORE);
    }
    if (process_right >= 0)
    {
        // std::cout << "Process " << mpi_rank << " sending to " << process_right << std::endl;
        MPI_Send(&u[N * (local_size - 1)], 1, vertical_slice_type, process_right, mpi_rank, comm2D);
    }
    if (process_left >= 0)
    {
        // std::cout << "Process " << mpi_rank << " receiving from " << process_left << std::endl;
        MPI_Recv(left,  N*local_size, MPI_DOUBLE, process_left, process_left, comm2D, MPI_STATUS_IGNORE);
    }
}

void ParallelHeat::exchange_rhs_z()
{
    int source, process_up, process_down;
    int direction = 1;
    MPI_Datatype horizontal_slice_type;
    MPI_Type_vector(local_size, N, N, MPI_DOUBLE, &horizontal_slice_type);
    MPI_Type_commit(&horizontal_slice_type);
    MPI_Cart_shift(comm2D, direction, 1, &source, &process_down);
    MPI_Cart_shift(comm2D, direction, -1, &source, &process_up);
    if (process_down >= 0)
    {
        // std::cout << "Process " << mpi_rank << " sending to " << process_down << std::endl;
        MPI_Send(&u[N * local_size * (local_size - 1)], 1, horizontal_slice_type, process_down, mpi_rank, comm2D);
    }
    if (process_up >= 0)
    {
        // std::cout << "Process " << mpi_rank << " receiving from " << process_up << std::endl;
        MPI_Recv(upper, N*local_size, MPI_DOUBLE, process_up, process_up, comm2D, MPI_STATUS_IGNORE);
    }
    if (process_up >= 0)
    {
        // std::cout << "Process " << mpi_rank << " sending to " << process_up << std::endl;
        MPI_Send(u, 1, horizontal_slice_type, process_up, mpi_rank, comm2D);
    }
    if (process_down >= 0)
    {
        // std::cout << "Process " << mpi_rank << " receiving from " << process_down << std::endl;
        MPI_Recv(lower, N*local_size, MPI_DOUBLE, process_down, process_down, comm2D, MPI_STATUS_IGNORE);
    }
}

void ParallelHeat::bcs_x()
{

    // Boundary conditions for x
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            rhs[index(0, j, k)] -= (exact(0.0, (coord[0]*local_size + j + 1)*dx, (coord[1]*local_size + k + 1)*dx, time + dt) -
                                    exact(0.0, (coord[0]*local_size + j + 1)*dx, (coord[1]*local_size + k + 1)*dx, time)) *
                                   offDiagCoeff;
            rhs[index(N - 1, j, k)] -= (exact(1.0, (coord[0]*local_size + j + 1)*dx, (coord[1]*local_size + k + 1)*dx, time + dt) -
                                        exact(1.0, (coord[0]*local_size + j + 1)*dx, (coord[1]*local_size + k + 1)*dx, time)) *
                                       offDiagCoeff;
        }
    }
}

void ParallelHeat::bcs_y()
{
    // Boundary conditions for x
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            rhs[index(0, j, k)] -= (exact((coord[0]*local_size + j + 1)*dx, 0.0, (coord[1]*local_size + k + 1)*dx, time + dt) -
                                    exact((coord[0]*local_size + j + 1)*dx, 0.0, (coord[1]*local_size + k + 1)*dx, time)) *
                                   offDiagCoeff;
            rhs[index(N - 1, j, k)] -= (exact((coord[0]*local_size + j + 1)*dx, 1.0, (coord[1]*local_size + k + 1)*dx, time + dt) -
                                        exact((coord[0]*local_size + j + 1)*dx, 1.0, (coord[1]*local_size + k + 1)*dx, time)) *
                                       offDiagCoeff;
        }
    }
}

void ParallelHeat::bcs_z()
{
    // Boundary conditions for x
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            rhs[index(0, j, k)] -= (exact((coord[0]*local_size + j + 1)*dx, (coord[1]*local_size + k + 1)*dx, 0.0, time + dt) -
                                    exact((coord[0]*local_size + j + 1)*dx, (coord[1]*local_size + k + 1)*dx, 0.0, time)) *
                                   offDiagCoeff;
            rhs[index(N - 1, j, k)] -= (exact((coord[0]*local_size + j + 1)*dx, (coord[1]*local_size + k + 1)*dx, 1.0, time + dt) -
                                        exact((coord[0]*local_size + j + 1)*dx, (coord[1]*local_size + k + 1)*dx, 1.0, time)) *
                                       offDiagCoeff;
        }
    }
}

void ParallelHeat::thomas()
{
    double *b_thomas;
    MPI_Alloc_mem(local_size * local_size * N * sizeof(double), MPI_INFO_NULL, &b_thomas);

    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                b_thomas[index(i, j, k)] = diagCoeff;
            }
        }
    }

    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 1; i < N; i++)
            {
                w = offDiagCoeff / b_thomas[index(i - 1, j, k)];
                b_thomas[index(i, j, k)] = b_thomas[index(i, j, k)] - w * offDiagCoeff;
                rhs[index(i, j, k)] -= rhs[index(i - 1, j, k)] * w;
            }
        }
    }

    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            rhs[index(N - 1, j, k)] /= b_thomas[index(N - 1, j, k)];
            for (int i = N - 2; i >= 0; i--)
            {
                rhs[index(i, j, k)] = (rhs[index(i, j, k)] - offDiagCoeff * rhs[index(i + 1, j, k)]) / b_thomas[index(i, j, k)];
            }
        }
    }
}


void ParallelHeat::rotate(int direction)
{
    int switcher = 2;
    int i = 0;
    int displacement;
    int source, dest;
    int coordswitch[dimensions];
    for (int stride = 1; stride < dims[direction]; stride++)
    {
        if (stride > switcher)
        {
            switcher *= 2;
            i = 0;
        }
        if (coord[direction] % switcher < switcher / 2)
        {
            if (coord[direction] % (switcher / 2) >= switcher / 2 - i)
                displacement = i;
            else
                displacement = stride;
        }
        else
        {
            if (coord[direction] % (switcher / 2) < i)
                displacement = -i;
            else
                displacement = -stride;
        }
        i++;
        MPI_Cart_shift(comm2D, direction, displacement, &source, &dest);
        MPI_Cart_coords(comm2D, dest, dimensions, coordswitch);
        for (unsigned int j = 0; j < local_size * local_size; j++)
        {
            MPI_Sendrecv_replace(rhs + local_size * coordswitch[direction] + j * N, local_size, MPI_DOUBLE,
                                 dest, dest, dest, mpi_rank,
                                 comm2D, MPI_STATUS_IGNORE);
        }
    }
}

void ParallelHeat::change_direction_x_to_y()
{
    
    rotate(direction);
    direction = (direction + 1) % 2;
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int p = 0; p < local_size; p++)
        {
            for (unsigned int j = 0; j < N / local_size; j++)
            {
                for (unsigned int i = 1; i + p < local_size; i++)
                {
                    std::swap(rhs[i + j * local_size + p * N + k * N * local_size + p], rhs[i + j * local_size + p * N + k * N * local_size + N * i - i + p]);
                }
            }
        }
    }
}

void ParallelHeat::change_direction_y_to_z()
{
    
    rotate(direction);
    direction = (direction + 1) % 2;
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int p = 0; p < local_size; p++)
        {
            for (unsigned int j = 0; j < N / local_size; j++)
            {
                for (unsigned int i = 1; i + k < local_size; i++)
                {
                    std::swap(rhs[i + j * local_size + p * N + k * N * local_size + k], rhs[i + j * local_size + p * N + k * N * local_size + local_size * N * i - i + k]);
                }
            }
        }
    }
}

void ParallelHeat::change_direction_z_to_x()
{
    
    rotate(direction);
    direction = (direction + 1) % 2;
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int p = 0; p < local_size; p++)
        {
            for (unsigned int j = 0; j < N / local_size; j++)
            {
                for (unsigned int i = 1; i + p < local_size; i++)
                {
                    std::swap(rhs[i + j * local_size + p * N + k * N * local_size + p], rhs[i + j * local_size + p * N + k * N * local_size + N * i - i + p]);
                }
            }
        }
    }
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int p = 1; p + k < local_size; p++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                std::swap(rhs[i + p * N + k * N * local_size + k * N], rhs[i + p * N + k * N * local_size + p * N * local_size - p * N + k * N]);
            }
        }
    }
}

void ParallelHeat::sum_rhs_to_u()
{
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                u[index(i, j, k)] += rhs[index(i, j, k)];
            }
        }
    }
}

void ParallelHeat::print_error()
{

    double error = 0.0;
    double error_tot = 0.0;
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                error += std::pow(u[index(i, j, k)] -
                                      exact((i + 1)*dx, (j + coord[0]*local_size + 1)*dx, (k + coord[1]*local_size + 1) * dx, 1.0),
                                  2.0);
            }
        }
    }

    MPI_Reduce(&error, &error_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    std::cout << "Error of process " << mpi_rank << " is: " << std::sqrt(error) / (double)(N * local_size * local_size) << std::endl;

    if (mpi_rank == 0)
    {
        std::cout << "Error with " << N << " x " << N << " x " << N
                  << " spatial elements and " << TSteps << " time steps is: " << std::sqrt(error_tot) / (double)(N * N * N) << std::endl;
    }
}