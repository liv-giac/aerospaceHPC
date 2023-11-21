#include <iostream>
#include <cmath>
#include <mpi.h>

class ParallelHeat
{
public:
    ParallelHeat(int mpi_rank, int mpi_size)
    {
        this->mpi_rank = mpi_rank;
        mpi_side_size = std::sqrt(mpi_size);
        N += N % mpi_side_size;
        local_size = N / mpi_side_size;
        TSteps = 100;
        dx2 = 1.0 / ((double)(N + 1) * (double)(N + 1));
        dx = 1.0 / (double)(N + 1);
        dt = 1.0 / (double)TSteps;

        offDiagCoeff = -dt / (2.0 * dx2);
        diagCoeff = 1.0 + dt / dx2;
        // w = offDiagCoeff / diagCoeff;
        b = diagCoeff; //- (offDiagCoeff * w);

        time = 0.0;

        z_block_start_idx = mpi_rank % mpi_side_size * local_size;
        y_block_start_idx = mpi_rank / mpi_side_size * local_size;

        MPI_Alloc_mem(local_size * local_size * N * sizeof(double), MPI_INFO_NULL, &u);
        MPI_Alloc_mem(local_size * local_size * N * sizeof(double), MPI_INFO_NULL, &rhs);

        // u = (double *)malloc(local_size * local_size * N * sizeof(double));
        // rhs = (double *)malloc(local_size * local_size * N * sizeof(double));

        MPI_Alloc_mem(local_size * N * sizeof(double), MPI_INFO_NULL, &right);
        MPI_Alloc_mem(local_size * N * sizeof(double), MPI_INFO_NULL, &left);
        MPI_Alloc_mem(local_size * N * sizeof(double), MPI_INFO_NULL, &upper);
        MPI_Alloc_mem(local_size * N * sizeof(double), MPI_INFO_NULL, &lower);

        // right = (double *)malloc(local_size * N * sizeof(double));
        // left = (double *)malloc(local_size * N * sizeof(double));
        // upper = (double *)malloc(local_size * N * sizeof(double));
        // lower = (double *)malloc(local_size * N * sizeof(double));

        if (mpi_rank == 0)
        {
            std::cout << "========================================================" << std::endl;
            std::cout << "Initializing parallel heat equation solver..." << std::endl;
            std::cout << "MPI size: " << mpi_size << std::endl;
            std::cout << "MPI side size: " << mpi_side_size << std::endl;
            std::cout << "N: " << N << std::endl;
            std::cout << "Time steps: " << TSteps << std::endl;
            std::cout << "dx: " << dx << std::endl;
            std::cout << "Local size: " << local_size << std::endl;
            std::cout << "========================================================" << std::endl;
        }

        int wrap_around[dimensions];
        int reorder;
        int my_cart_rank;

        // MPI_Alloc_mem(dimensions * sizeof(int), MPI_INFO_NULL, &dims);
        // MPI_Alloc_mem(dimensions * sizeof(int), MPI_INFO_NULL, &coord);

        dims[0] = dims[1] = dims[2] = 0;

        /* create cartesian topology for processes */
        MPI_Dims_create(mpi_size, dimensions, dims);
        std::cout << "PW[" << mpi_rank << "]/[" << mpi_size << "]: PEdims = [" << dims[0] << " x " << dims[1] << " x " << dims[2] << "]" << std::endl;
        /* create cartesian mapping */
        wrap_around[0] = wrap_around[1] = wrap_around[2] = 0; // set periodicity
        reorder = 1;
        int ierr = 0;
        ierr = MPI_Cart_create(MPI_COMM_WORLD, dimensions, dims,
                               wrap_around, reorder, &comm2D);
        if (ierr != 0)
            std::cout << "ERROR[" << ierr << "] creating CART" << std::endl;

        /* find my coordinates in the cartesian communicator group */
        MPI_Cart_coords(comm2D, mpi_rank, dimensions, coord);
        /* use my coords to find my rank in cartesian group*/
        MPI_Cart_rank(comm2D, coord, &my_cart_rank);
        std::cout << "PW[" << mpi_rank << "]: my_cart_rank PCM[" << my_cart_rank << "], my coords = (" << coord[0] << "," << coord[1] << "," << coord[2] << ")" << std::endl;
    };
    ~ParallelHeat()
    {
        MPI_Free_mem(right);
        MPI_Free_mem(left);
        MPI_Free_mem(upper);
        MPI_Free_mem(lower);
        MPI_Free_mem(u);
        MPI_Free_mem(rhs);
    };

    double time;
    unsigned int TSteps;
    double dt;

    void printRhs();

    void initialize_u();

    void build_rhs();

    void exchange_rhs_y();

    void exchange_rhs_z();

    void bcs_x();

    void bcs_y();

    void bcs_z();

    void thomas();

    void change_direction_x_to_y();

    void change_direction_y_to_z();

    void change_direction_z_to_x();

    void sum_rhs_to_u();

    void print_error();

    void rotate(int direction);

private:
    // Private member functions
    inline unsigned int index(unsigned int i, unsigned int j, unsigned int k);

    inline double force(double x, double y, double z, double t);

    inline double exact(double x, double y, double z, double t);

    // void rotate_block_x_to_y(double *block);

    // void rotate_block_y_to_z(double *block);

    // void rotate_block_z_to_x(double *block);

    // MPI-related variables
    unsigned int mpi_rank;
    unsigned int mpi_side_size;
    unsigned int local_size;
    int direction = 0;
    unsigned int z_block_start_idx;
    unsigned int y_block_start_idx;
    MPI_Comm comm2D;
    int dimensions = 3;
    int dims[3];
    int coord[3];

    // Parameters
    unsigned int N = 4;
    double dx2;
    double dx;

    double offDiagCoeff;
    double diagCoeff;
    double w;
    double b;

    // Arrays
    double *u;
    double *rhs;

    // Arrays for communication
    double *right;
    double *left;
    double *upper;
    double *lower;
};