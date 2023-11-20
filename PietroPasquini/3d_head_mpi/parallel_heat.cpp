#include "parallel_heat.hpp"

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
                rhs[index(i, j, k)] = dt * (force(x, y, z, time + 0.5 * dt) -
                                            (1.0 / dx2) * 6.0 * u[index(i, j, k)]);
                if (k > 0)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i, j, k - 1)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * left[index(j, k, 0)];
                if (k < local_size - 1)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i, j, k + 1)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * right[index(j, k, 0)];
                if (j > 0)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i, j - 1, k)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * lower[index(i, k, 0)];
                if (j < local_size - 1)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i, j + 1, k)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * upper[index(i, k, 0)];
                if (x > dx)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i - 1, j, k)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * exact(x, y, 0.0, time);
                if (x < 1.0 - dx)
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * u[index(i + 1, j, k)];
                else
                    rhs[index(i, j, k)] += dt * 1.0 / dx2 * exact(x, y, 1.0, time);
            }
            x = 0.0;
        }
        y = y_block_start_idx * dx;
    }
}

void ParallelHeat::message_passing()
{
    // If there is only one process, fill vectors with exact values
    if (mpi_side_size == 1)
    {
        double x = z_block_start_idx * dx;
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                x += dx;
                lower[index(i, j, 0)] = exact((double)(i + 1) * dx, 0.0, x, time);
                upper[index(i, j, 0)] = exact((double)(i + 1) * dx, 1.0, x, time);
            }
            x = z_block_start_idx * dx;
        }
        x = y_block_start_idx * dx;
        for (unsigned int i = 0; i < N; i++)
        {
            for (unsigned int j = 0; j < local_size; j++)
            {
                x += dx;
                left[index(i, j, 0)] = exact((double)(i + 1) * dx, x, 0.0, time);
                right[index(i, j, 0)] = exact((double)(i + 1) * dx, x, 1.0, time);
            }
            x = y_block_start_idx * dx;
        }
        return;
    }
    double x = z_block_start_idx * dx;
    // Create MPI_Datatype for horizontal slice
    MPI_Datatype horizontal_slice_type;
    MPI_Type_vector(local_size, N, N * local_size, MPI_DOUBLE, &horizontal_slice_type);
    MPI_Type_commit(&horizontal_slice_type);
    // Send and recieve top
    if (mpi_rank < mpi_side_size) // bottom row (tag 0)
    {
        MPI_Sendrecv(&u[index(0, local_size - 1, 0)], 1, horizontal_slice_type, mpi_rank + mpi_side_size, 0,
                     upper, N * local_size, MPI_DOUBLE, mpi_rank + mpi_side_size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Fill other vectors with exact values
        x = z_block_start_idx * dx;
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                x += dx;
                lower[index(i, j, 0)] = exact((double)(i + 1) * dx, 0.0, x, time);
            }
            x = z_block_start_idx * dx;
        }
        if (mpi_rank == 0) // bottom left corner
        {
            x = y_block_start_idx * dx;
            for (unsigned int i = 0; i < N; i++)
            {
                for (unsigned int j = 0; j < local_size; j++)
                {
                    x += dx;
                    left[index(i, j, 0)] = exact((double)(i + 1) * dx, x, 0.0, time);
                }
                x = y_block_start_idx * dx;
            }
        }
        if (mpi_rank == mpi_side_size - 1) // bottom right corner
        {
            x = y_block_start_idx * dx;
            for (unsigned int j = 0; j < local_size; j++)
            {
                for (unsigned int i = 0; i < N; i++)
                {
                    x += dx;
                    right[index(i, j, 0)] = exact((double)(i + 1) * dx, x, 1.0, time);
                }
                x = y_block_start_idx * dx;
            }
        }
    }
    else if (mpi_rank >= mpi_side_size * (mpi_side_size - 1)) // top row (tag 1) // Send and recieve bottom
    {
        MPI_Sendrecv(u, 1, horizontal_slice_type, mpi_rank - mpi_side_size, 1,
                     lower, N * local_size, MPI_DOUBLE, mpi_rank - mpi_side_size, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Fill other vectors with exact values
        x = z_block_start_idx * dx;
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                x += dx;
                upper[index(i, j, 0)] = exact((double)(i + 1) * dx, 1.0, x, time);
            }
            x = z_block_start_idx * dx;
        }
        if (mpi_rank == mpi_side_size * mpi_side_size - 1) // top right corner
        {
            x = y_block_start_idx * dx;
            for (unsigned int j = 0; j < local_size; j++)
            {
                for (unsigned int i = 0; i < N; i++)
                {
                    x += dx;
                    right[index(i, j, 0)] = exact((double)(i + 1) * dx, x, 1.0, time);
                }
                x = y_block_start_idx * dx;
            }
        }
        if (mpi_rank == mpi_side_size * (mpi_side_size - 1)) // top left corner
        {
            x = y_block_start_idx * dx;
            for (unsigned int j = 0; j < local_size; j++)
            {
                for (unsigned int i = 0; i < N; i++)
                {
                    x += dx;
                    left[index(i, j, 0)] = exact((double)(i + 1) * dx, x, 0.0, time);
                }
                x = y_block_start_idx * dx;
            }
        }
    }
    else // center rows (tag 0 and 1)
    {
        MPI_Sendrecv(u, 1, horizontal_slice_type, mpi_rank - mpi_side_size, 1,
                     lower, N * local_size, MPI_DOUBLE, mpi_rank - mpi_side_size, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&u[index(0, local_size - 1, 0)], 1, horizontal_slice_type, mpi_rank + mpi_side_size, 0,
                     upper, N * local_size, MPI_DOUBLE, mpi_rank + mpi_side_size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // Send and recieve left
    if (mpi_rank % mpi_side_size == 0) // left column (tag 2)
    {
        MPI_Sendrecv(&u[index(0, 0, local_size - 1)], N * local_size, MPI_DOUBLE, mpi_rank + 1, 2,
                     right, N * local_size, MPI_DOUBLE, mpi_rank + 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Fill other vectors with exact values
        x = y_block_start_idx * dx;
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                x += dx;
                left[index(i, j, 0)] = exact((double)(i + 1) * dx, x, 0.0, time);
            }
            x = y_block_start_idx * dx;
        }
    }
    else if (mpi_rank % mpi_side_size == mpi_side_size - 1) // right column (tag 3) // Send and recieve right
    {
        MPI_Sendrecv(u, N * local_size, MPI_DOUBLE, mpi_rank - 1, 3,
                     left, N * local_size, MPI_DOUBLE, mpi_rank - 1, 2,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // Fill other vectors with exact values
        x = y_block_start_idx * dx;
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                x += dx;
                right[index(i, j, 0)] = exact((double)(i + 1) * dx, x, 1.0, time);
            }
            x = y_block_start_idx * dx;
        }
    }
    else // center columns
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < N; i++)
            {
                left[index(i, j, 0)] = u[index(i, j, z_block_start_idx)];
                right[index(i, j, 0)] = u[index(i, j, z_block_start_idx + local_size - 1)];
            }
        }
        MPI_Sendrecv(u, N * local_size, MPI_DOUBLE, mpi_rank - 1, 3,
                     left, N * local_size, MPI_DOUBLE, mpi_rank - 1, 2,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&u[index(0, 0, local_size - 1)], N * local_size, MPI_DOUBLE, mpi_rank + 1, 2,
                     right, N * local_size, MPI_DOUBLE, mpi_rank + 1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void ParallelHeat::bcs_x()
{
    double y = y_block_start_idx * dx;
    double z = z_block_start_idx * dx;

    // Boundary conditions for x
    for (unsigned int j = 0; j < local_size; j++)
    {
        y += dx;
        for (unsigned int k = 0; k < local_size; k++)
        {
            z += dx;
            rhs[index(0, j, k)] -= (exact(0.0, y, z, time + dt) -
                                    exact(0.0, y, z, time)) *
                                   offDiagCoeff;
            rhs[index(N - 1, j, k)] -= (exact(1.0, y, z, time + dt) -
                                        exact(1.0, y, z, time)) *
                                       offDiagCoeff;
        }
        z = z_block_start_idx * dx;
    }
}

void ParallelHeat::bcs_y()
{
    double y = z_block_start_idx * dx;
    double z = y_block_start_idx * dx;

    // Boundary conditions for y
    for (unsigned int j = 0; j < local_size; j++)
    {
        y += dx;
        for (unsigned int k = 0; k < local_size; k++)
        {
            z += dx;
            rhs[index(0, j, k)] -= (exact(0.0, y, z, time + dt) -
                                    exact(0.0, y, z, time)) *
                                   offDiagCoeff;
            rhs[index(N - 1, j, k)] -= (exact(1.0, y, z, time + dt) -
                                        exact(1.0, y, z, time)) *
                                       offDiagCoeff;
        }
        z = y_block_start_idx * dx;
    }
}

void ParallelHeat::bcs_z()
{
    double y = y_block_start_idx * dx;
    double z = z_block_start_idx * dx;

    // Boundary conditions for z
    for (unsigned int j = 0; j < local_size; j++)
    {
        y += dx;
        for (unsigned int k = 0; k < local_size; k++)
        {
            z += dx;
            rhs[index(0, j, k)] -= (exact(0.0, y, z, time + dt) -
                                    exact(0.0, y, z, time)) *
                                   offDiagCoeff;
            rhs[index(N - 1, j, k)]-= (exact(1.0, y, z, time + dt) -
                                        exact(1.0, y, z, time)) *
                                       offDiagCoeff;
        }
        z = z_block_start_idx * dx;
    }
}

void ParallelHeat::thomas()
{
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 1; i < N; i++)
            {
                rhs[index(i, j, k)] -= rhs[index(i - 1, j, k)] * w;
            }
            rhs[index(N - 1, j, k)] /= b;
            for (int i = N - 2; i >= 0; i--)
            {
                rhs[index(i, j, k)] = (rhs[index(i, j, k)] - offDiagCoeff * rhs[index(i + 1, j, k)]) / b;
            }
        }
    }
}

void ParallelHeat::rotate_block_x_to_y(double *block)
{
    double temp;
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < local_size; i++)
            {
                temp = block[i + j * local_size + k * local_size * local_size];
                block[i + j * local_size + k * local_size * local_size] =
                    block[j + i * local_size + k * local_size * local_size];
                block[j + i * local_size + k * local_size * local_size] = temp;
            }
        }
    }
}

void ParallelHeat::rotate_block_y_to_z(double *block)
{
    double temp;
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < local_size; i++)
            {
                temp = block[i + j * local_size + k * local_size * local_size];
                block[i + j * local_size + k * local_size * local_size] =
                    block[k + j * local_size + i * local_size * local_size];
                block[k + j * local_size + i * local_size * local_size] = temp;
            }
        }
    }
}

void ParallelHeat::rotate_block_z_to_x(double *block)
{
    double temp;
    for (unsigned int k = 0; k < local_size; k++)
    {
        for (unsigned int j = 0; j < local_size; j++)
        {
            for (unsigned int i = 0; i < local_size; i++)
            {
                temp = block[i + j * local_size + k * local_size * local_size];
                block[i + j * local_size + k * local_size * local_size] =
                    block[j + k * local_size + i * local_size * local_size];
                block[j + k * local_size + i * local_size * local_size] = temp;
            }
        }
    }
}

void ParallelHeat::change_direction_x_to_y()
{
    double *buffer;
    MPI_Alloc_mem(local_size * local_size * local_size * sizeof(double), MPI_INFO_NULL, &buffer);
    int column = mpi_rank % mpi_side_size;
    int row = mpi_rank / mpi_side_size;
    for (unsigned int pos = 0; pos < mpi_side_size; pos++)
    {

        for (unsigned int k = 0; k < local_size; k++)
        {
            for (unsigned int j = 0; j < local_size; j++)
            {
                for (unsigned int i = 0; i < local_size; i++)
                {
                    buffer[i + j * local_size + k * local_size * local_size] = rhs[index(i + pos * local_size, j, k)];
                }
            }
        }
        if (mpi_rank != column + pos * mpi_side_size) //  not owned blocks
        {
            MPI_Sendrecv_replace(buffer, local_size * local_size * local_size, MPI_DOUBLE, column + pos * mpi_side_size, 4,
                                 column + pos * mpi_side_size, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        rotate_block_x_to_y(buffer);
        for (unsigned int k = 0; k < local_size; k++)
        {
            for (unsigned int j = 0; j < local_size; j++)
            {
                for (unsigned int i = 0; i < local_size; i++)
                {
                    rhs[index(i + pos * local_size, j, k)] = buffer[i + j * local_size + k * local_size * local_size];
                }
            }
        }
    }
    MPI_Free_mem(buffer);
}

void ParallelHeat::change_direction_y_to_z()
{
    double *buffer;
    MPI_Alloc_mem(local_size * local_size * local_size * sizeof(double), MPI_INFO_NULL, &buffer);
    double column = mpi_rank % mpi_side_size;
    int row = mpi_rank / mpi_side_size;
    for (unsigned int pos = 0; pos < mpi_side_size; pos++)
    {
        for (unsigned int k = 0; k < local_size; k++)
        {
            for (unsigned int j = 0; j < local_size; j++)
            {
                for (unsigned int i = 0; i < local_size; i++)
                {
                    buffer[i + j * local_size + k * local_size * local_size] = rhs[index(i + pos * local_size, j, k)];
                }
            }
        }
        if (mpi_rank != column + pos * mpi_side_size) //  not owned blocks
        {
            MPI_Sendrecv_replace(buffer, local_size * local_size * local_size, MPI_DOUBLE, column + pos * mpi_side_size, 4,
                                 column + pos * mpi_side_size, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        rotate_block_y_to_z(buffer);
        for (unsigned int k = 0; k < local_size; k++)
        {
            for (unsigned int j = 0; j < local_size; j++)
            {
                for (unsigned int i = 0; i < local_size; i++)
                {
                    rhs[index(i + pos * local_size, j, k)] = buffer[i + j * local_size + k * local_size * local_size];
                }
            }
        }
    }
    MPI_Free_mem(buffer);
}

void ParallelHeat::change_direction_z_to_x()
{
    double *buffer;
    MPI_Alloc_mem(local_size * local_size * local_size * sizeof(double), MPI_INFO_NULL, &buffer);
    double column = mpi_rank % mpi_side_size;
    int row = mpi_rank / mpi_side_size;
    for (unsigned int pos = 0; pos < mpi_side_size; pos++)
    {

        for (unsigned int k = 0; k < local_size; k++)
        {
            for (unsigned int j = 0; j < local_size; j++)
            {
                for (unsigned int i = 0; i < local_size; i++)
                {
                    buffer[i + j * local_size + k * local_size * local_size] = rhs[index(i + pos * local_size, j, k)];
                }
            }
        }
        if (mpi_rank != column + pos * mpi_side_size) //  not owned blocks
        {
            MPI_Sendrecv_replace(buffer, local_size * local_size * local_size, MPI_DOUBLE, column + pos * mpi_side_size, 4,
                                 column + pos * mpi_side_size, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        rotate_block_z_to_x(buffer);
        for (unsigned int k = 0; k < local_size; k++)
        {
            for (unsigned int j = 0; j < local_size; j++)
            {
                for (unsigned int i = 0; i < local_size; i++)
                {
                    rhs[index(i + pos * local_size, j, k)] = buffer[i + j * local_size + k * local_size * local_size];
                }
            }
        }
    }
    MPI_Free_mem(buffer);
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
    double x = 0.0;
    double y = y_block_start_idx * dx;
    double z = z_block_start_idx * dx;
    double error = 0.0;
    double error_tot = 0.0;
    for (unsigned int k = 0; k < local_size; k++)
    {
        x += dx;
        for (unsigned int j = 0; j < local_size; j++)
        {
            y += dx;
            for (unsigned int i = 0; i < N; i++)
            {
                x += dx;
                error += std::pow(u[index(i, j, k)] -
                                      exact(x, y, z, 1.0),
                                  2.0);
            }
            x = 0.0;
        }
        y = y_block_start_idx * dx;
    }
    error = std::sqrt(error) / (double)(N * local_size * local_size);

    MPI_Reduce(&error, &error_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    std::cout << "Error of process " << mpi_rank << " is: " << error << std::endl;


    if (mpi_rank == 0)
    {
        std::cout << "Error with " << N << " x " << N << " x " << N
                  << " spatial elements and " << TSteps << " time steps is: " << error_tot << std::endl;
    }
}