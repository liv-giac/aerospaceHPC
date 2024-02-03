/* mpi-cart-2D.c -- test basic -cartesian functions
 * Written by Mary Thomas- Updated Mar, 2015
 * Based loosely on code from Pacheco’97,
 * Chap 7, pp. 121 & ff in PPMPI */
#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <bits/stdc++.h>

using vector = std::vector<double>;

void rotate(std::vector<double> &rhs_local, const int N, const int size, int *dims, int direction, int *coord, int dimensions, const int &my_rank, MPI_Comm comm2D)
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
        for (int j = 0; j < size * size; j++)
        {
            MPI_Sendrecv_replace(rhs_local.data() + size * coordswitch[direction] + j * N, size, MPI_DOUBLE,
                                 dest, 0, dest, 0,
                                 comm2D, MPI_STATUS_IGNORE);
        }
    }
}

void exchange_rhs_y(vector &rhs_local, vector &right, vector &left, const int N, const int size, const int &mpi_rank, MPI_Comm comm2D)
{
    int source, process_right, process_left;
    int direction = 0;
    MPI_Datatype vertical_slice_type;
    MPI_Type_vector(size, N, N * size, MPI_DOUBLE, &vertical_slice_type);
    MPI_Type_commit(&vertical_slice_type);
    MPI_Cart_shift(comm2D, direction, 1, &source, &process_right);
    MPI_Cart_shift(comm2D, direction, -1, &source, &process_left);
    if (process_left >= 0)
    {
        // std::cout << "Process " << mpi_rank << " sending to " << process_left << std::endl;
        MPI_Send(rhs_local.data(), 1, vertical_slice_type, process_left, mpi_rank, comm2D);
    }
    if (process_right >= 0)
    {
        // std::cout << "Process " << mpi_rank << " receiving from " << process_right << std::endl;
        MPI_Recv(right.data(), N * size, MPI_DOUBLE, process_right, process_right, comm2D, MPI_STATUS_IGNORE);
    }
    if (process_right >= 0)
    {
        // std::cout << "Process " << mpi_rank << " sending to " << process_right << std::endl;
        MPI_Send(rhs_local.data() + N * (size - 1), 1, vertical_slice_type, process_right, mpi_rank, comm2D);
    }
    if (process_left >= 0)
    {
        // std::cout << "Process " << mpi_rank << " receiving from " << process_left << std::endl;
        MPI_Recv(left.data(), N * size, MPI_DOUBLE, process_left, process_left, comm2D, MPI_STATUS_IGNORE);
    }
}

void exchange_rhs_z(vector &upper, vector &lower, vector &rhs_local, const int N, const int size, const int &mpi_rank, MPI_Comm comm2D)
{
    int source, process_up, process_down;
    int direction = 1;
    MPI_Datatype horizontal_slice_type;
    MPI_Type_vector(size, N, N, MPI_DOUBLE, &horizontal_slice_type);
    MPI_Type_commit(&horizontal_slice_type);
    MPI_Cart_shift(comm2D, direction, 1, &source, &process_down);
    MPI_Cart_shift(comm2D, direction, -1, &source, &process_up);
    if (process_down >= 0)
    {
        // std::cout << "Process " << mpi_rank << " sending to " << process_down << std::endl;
        MPI_Send(rhs_local.data() + N * size * (size - 1), 1, horizontal_slice_type, process_down, mpi_rank, comm2D);
    }
    if (process_up >= 0)
    {
        // std::cout << "Process " << mpi_rank << " receiving from " << process_up << std::endl;
        MPI_Recv(upper.data(), N * size, MPI_DOUBLE, process_up, process_up, comm2D, MPI_STATUS_IGNORE);
    }
    if (process_up >= 0)
    {
        // std::cout << "Process " << mpi_rank << " sending to " << process_up << std::endl;
        MPI_Send(rhs_local.data(), 1, horizontal_slice_type, process_up, mpi_rank, comm2D);
    }
    if (process_down >= 0)
    {
        // std::cout << "Process " << mpi_rank << " receiving from " << process_down << std::endl;
        MPI_Recv(lower.data(), N * size, MPI_DOUBLE, process_down, process_down, comm2D, MPI_STATUS_IGNORE);
    }
}

void printRhs(const vector &rhs_local, const int &N, const int &size, const int &my_rank, MPI_Comm comm2D)
{

    std::vector<double> rhs;
    rhs.resize(N * N * N);

    MPI_Gather(rhs_local.data(), N * size * size, MPI_DOUBLE, rhs.data(), N * size * size, MPI_DOUBLE, 0, comm2D);

    // Print rhs_local for every processor
    if (my_rank == 0)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    std::cout << rhs[i * N * N + j * N + k] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
}

int main(int argc, char *argv[])
{
    int N = 4;
    int size;
    vector rhs_local;
    vector upper, lower, right, left;

    int dimensions = 3,
        ierr;
    int p, my_rank, my_cart_rank, stride;
    MPI_Comm comm2D;
    int dims[dimensions], coord[dimensions];
    int wrap_around[dimensions];
    int reorder;

    /* start up initial MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /* process command line arguments*/
    dims[0] = dims[1] = dims[2] = 0;

    /* create cartesian topology for processes */
    MPI_Dims_create(p, dimensions, dims);
    if (my_rank == 0)
        std::cout << "PW[" << my_rank << "]/[" << p << "]: PEdims = [" << dims[0] << " x " << dims[1] << " x " << dims[2] << "]" << std::endl;

    /* create cartesian mapping */
    wrap_around[0] = wrap_around[1] = wrap_around[2] = 0; // set periodicity
    reorder = 1;
    ierr = 0;
    ierr = MPI_Cart_create(MPI_COMM_WORLD, dimensions, dims,
                           wrap_around, reorder, &comm2D);
    if (ierr != 0)
        std::cout << "ERROR[" << ierr << "] creating CART" << std::endl;
    size = N / dims[0];

    rhs_local.resize(size * size * N);
    upper.resize(size * N);
    lower.resize(size * N);
    right.resize(size * N);
    left.resize(size * N);

    if (my_rank == 0)
    {
        std::cout << "PW[" << my_rank << "]/[" << p << "]: PEdims = [" << dims[0] << " x " << dims[1] << " x " << dims[2] << "]" << std::endl;
        int source, right, left, up, down;
        MPI_Cart_shift(comm2D, 0, 1, &source, &right);
        MPI_Cart_shift(comm2D, 0, -1, &source, &left);
        MPI_Cart_shift(comm2D, 1, 1, &source, &down);
        MPI_Cart_shift(comm2D, 1, -1, &source, &up);
        std::cout << "PW[" << my_rank << "]: right = " << right << ", left = " << left << ", up = " << up << ", down = " << down << std::endl;
    }

    /* find my coordinates in the cartesian communicator group */
    MPI_Cart_coords(comm2D, my_rank, dimensions, coord);
    /* use my coords to find my rank in cartesian group*/
    MPI_Cart_rank(comm2D, coord, &my_cart_rank);
    std::cout << "PW[" << my_rank << "]: my_cart_rank PCM[" << my_cart_rank << "], my coords = (" << coord[0] << "," << coord[1] << "," << coord[2] << ")" << std::endl;

    // Initialize rhs_local
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < N; k++)
            {
                rhs_local[i * size * N + j * N + k] = i * N * N + j * N + k + coord[0] * size * N + coord[1] * size * N * N;
            }
        }
    }

    MPI_Barrier(comm2D);
    printRhs(rhs_local, N, size, my_rank, comm2D);
    MPI_Barrier(comm2D);

    // Rotate x to y
    if (my_rank == 0)
        std::cout << "Rotating x to y" << std::endl;
    rotate(rhs_local, N, size, dims, 0, coord, dimensions, my_rank, comm2D);
    // for(int k = 0; k < size; k++){
    //     j = 1;
    //     while(j < stride){
    //         std::swap(rhs_local[k * size * N + j], rhs_local[k * size * N + j + stride - 1]);
    //         j+=size;
    //     }
    // }
    for (int k = 0; k < size; k++)
    {
        for (int p = 0; p < size; p++)
        {
            for (int j = 0; j < N / size; j++)
            {
                for (int i = 1; i + p < size; i++)
                {
                    std::swap(rhs_local[i + j * size + p * N + k * N * size + p], rhs_local[i + j * size + p * N + k * N * size + N * i - i + p]);
                }
            }
        }
    }
    MPI_Barrier(comm2D);
    printRhs(rhs_local, N, size, my_rank, comm2D);
    MPI_Barrier(comm2D);

    // Rotate y to z
    if (my_rank == 0)
        std::cout << "Rotating y to z" << std::endl;
    rotate(rhs_local, N, size, dims, 1, coord, dimensions, my_rank, comm2D);
    for (int k = 0; k < size; k++)
    {
        for (int p = 0; p < size; p++)
        {
            for (int j = 0; j < N / size; j++)
            {
                for (int i = 1; i + k < size; i++)
                {
                    std::swap(rhs_local[i + j * size + p * N + k * N * size + k], rhs_local[i + j * size + p * N + k * N * size + size * N * i - i + k]);
                }
            }
        }
    }
    MPI_Barrier(comm2D);
    printRhs(rhs_local, N, size, my_rank, comm2D);
    MPI_Barrier(comm2D);

    // Rotate z to x
    if (my_rank == 0)
        std::cout << "Rotating z to x" << std::endl;
    rotate(rhs_local, N, size, dims, 1, coord, dimensions, my_rank, comm2D);
    // for (int k = 0; k < size; k++)
    // {
    //     for (int p = 0; p < size; p++)
    //     {
    //         for (int j = 0; j < N / size; j++)
    //         {
    //             for (int i = 1; i + p < size; i++)
    //             {
    //                 std::swap(rhs_local[i + j * size + p * N + k * N * size + p], rhs_local[i + j * size + p * N + k * N * size + N * i - i + p]);
    //             }
    //         }
    //     }
    // }
    // for (int k = 0; k < size; k++)
    // {
    //     for (int p = 1; p + k < size; p++)
    //     {

    //         for (int i = 0; i < N; i++)
    //         {
    //             std::swap(rhs_local[i + p * N + k * N * size + k * N], rhs_local[i + p * N + k * N * size + p * N * size - p * N + k * N]);
    //         }
    //     }
    // }
    MPI_Barrier(comm2D);
    printRhs(rhs_local, N, size, my_rank, comm2D);

    exchange_rhs_y(rhs_local, right, left, N, size, my_rank, comm2D);
    exchange_rhs_z(upper, lower, rhs_local, N, size, my_rank, comm2D);

    if (my_rank == 0)
    {
        std::cout << "Upper" << std::endl;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < N; j++)
            {
                std::cout << upper[i * N + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Lower" << std::endl;

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < N; j++)
            {
                std::cout << lower[i * N + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Right" << std::endl;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < N; j++)
            {
                std::cout << right[i * N + j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "Left" << std::endl;
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < N; j++)
            {
                std::cout << left[i * N + j] << " ";
            }
            std::cout << std::endl;
        }
    }

    MPI_Comm_free(&comm2D);
    MPI_Finalize();
    return 0;
} /* main */
