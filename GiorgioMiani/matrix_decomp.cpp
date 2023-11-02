/* mpi-cart-2D.c -- test basic -cartesian functions
 * Written by Mary Thomas- Updated Mar, 2015
 * Based loosely on code from Pacheco’97,
 * Chap 7, pp. 121 & ff in PPMPI */
#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include <vector>

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
        std::cout << "Processor: " << my_rank << " send to " << dest << std::endl;
        for (int j = 0; j < size * size; j++)
        {
            MPI_Sendrecv_replace(rhs_local.data() + size * coordswitch[direction] + j * N, size, MPI_DOUBLE,
                                 dest, 0, dest, 0,
                                 comm2D, MPI_STATUS_IGNORE);
        }
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

    int dimensions = 3, ierr;
    int p, my_rank, my_cart_rank;
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
    rotate(rhs_local, N, size, dims, 0, coord, dimensions, my_rank, comm2D);
    MPI_Barrier(comm2D);
    printRhs(rhs_local, N, size, my_rank, comm2D);
    MPI_Barrier(comm2D);

    // Rotate y to z
    rotate(rhs_local, N, size, dims, 1, coord, dimensions, my_rank, comm2D);
    MPI_Barrier(comm2D);
    printRhs(rhs_local, N, size, my_rank, comm2D);
    MPI_Barrier(comm2D);

    // Rotate z to x
    rotate(rhs_local, N, size, dims, 0, coord, dimensions, my_rank, comm2D);
    MPI_Barrier(comm2D);
    printRhs(rhs_local, N, size, my_rank, comm2D);
    MPI_Comm_free(&comm2D);
    MPI_Finalize();
    return 0;
} /* main */
