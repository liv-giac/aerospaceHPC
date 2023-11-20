#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include <array>
#include <string>
#include <mpi.h>

//------------------------------------------VARIABLES-----------------------------------------//
#define X 1.0
#define T 1.0
#define N 4
#define nt 100

std::vector<double> solution;
std::array<double, N> my_sin;
std::array<double, N> diag;
std::array<double, nt + 1> my_sin_time;
std::vector<double> local_sol;

double dx;
double dx2;
double dt;

int local_size_x;
int local_size_y;
int local_size_z;

int local_start_x = 0;
int local_start_y = 0;
int local_start_z = 0;

int local_end_x;
int local_end_y;
int local_end_z;

// MPI cart variables
int dimensions = 3;
int wrap_around[3];
int reorder;
int dims[3];
int coord[3];
MPI_Comm comm2D;

//-------------------------------------PARALLEL PROGRAMMA-------------------------------------//

//--------------------------------------THOMAS ALGORITHM--------------------------------------//
void thomas_x(const double exdiag, const int rank, const int size, const int start, const int end, const int direction)
{

    double prev;
    double next;
    std::vector<double> diag;
    double w;

    int source, next_process, prev_process;

    diag.resize(local_size_z);
    diag[local_start_x] = 1.0 + dt / dx2;

    MPI_Cart_shift(comm2D, direction, 1, &source, &next_process);
    MPI_Cart_shift(comm2D, direction, -1, &source, &prev_process);

    for (int i = local_start_z; i <= local_end_z; i++)
    {
        for (int j = local_start_y; j <= local_end_y; j++)
        {

            if (coord[direction] != 0)
            {
                MPI_Recv(&prev, 1, MPI_DOUBLE, prev_process, 1, comm2D, MPI_STATUS_IGNORE);
                MPI_Recv(&diag[local_start_x], 1, MPI_DOUBLE, prev_process, 2, comm2D, MPI_STATUS_IGNORE);
            }

            for (int k = local_start_x; k <= local_end_x; k++)
            {
                if (k != local_start_x)
                    w = exdiag / diag[k - 1];
                diag[k] = diag[k] - w * exdiag;
                local_sol[k + j * local_size_x + i * local_size_x * local_size_y] -= w * prev;
                prev = local_sol[k + j * local_size_x + i * local_size_x * local_size_y];
            }

            if (coord[direction] != size - 1)
            {
                MPI_Send(&prev, 1, MPI_DOUBLE, next_process, 1, comm2D);
                MPI_Send(&diag[local_end_x], 1, MPI_DOUBLE, next_process, 2, comm2D);
            }
        }
    }

    MPI_Barrier(comm2D);

    for (int i = local_start_z; i < local_end_z; i++)
    {
        for (int j = local_start_y; j < local_end_y; j++)
        {

            if (coord[direction] != dims[direction] - 1)
                MPI_Recv(&next, 1, MPI_DOUBLE, next_process, 1, comm2D, MPI_STATUS_IGNORE);
            else
                local_sol[local_end_x + j * local_size_x + i * local_size_x * local_size_y] /= diag[local_end_x];

            for (int k = local_end_x; k >= local_start_x; k--)
            {
                if (k != local_end_x)
                    local_sol[k + j * local_size_x + i * local_size_x * local_size_y] -= (exdiag * next) / diag[k];
                next = local_sol[k + j * local_size_x + i * local_size_x * local_size_y];
            }

            if (coord[direction] != 0)
                MPI_Send(&next, 1, MPI_DOUBLE, prev_process, 1, comm2D);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void thomas_y(const double exdiag, const int rank, const int size, const int start, const int end, const int direction)
{

    double prev;
    double next;
    std::vector<double> diag;
    double w;

    int source, next_process, prev_process;

    diag.resize(local_size_y);
    diag[local_start_y] = 1.0 + dt / dx2;

    MPI_Cart_shift(comm2D, direction, 1, &source, &next_process);
    MPI_Cart_shift(comm2D, direction, -1, &source, &prev_process);

    for (int i = local_start_z; i <= local_end_z; i++)
    {
        for (int j = local_start_x; j <= local_end_x; j++)
        {

            if (coord[direction] != 0)
            {
                MPI_Recv(&prev, 1, MPI_DOUBLE, prev_process, 1, comm2D, MPI_STATUS_IGNORE);
                MPI_Recv(&diag[local_start_y], 1, MPI_DOUBLE, prev_process, 2, comm2D, MPI_STATUS_IGNORE);
            }

            for (int k = local_start_y; k <= local_end_y; k++)
            {
                if (k != local_start_y)
                    w = exdiag / diag[k - 1];
                diag[k] = diag[k] - w * exdiag;
                local_sol[j + k * local_size_x + i * local_size_x * local_size_y] -= w * prev;
                prev = local_sol[j + k * local_size_x + i * local_size_x * local_size_y];
            }

            if (coord[direction] != size - 1)
            {
                MPI_Send(&prev, 1, MPI_DOUBLE, next_process, 1, comm2D);
                MPI_Send(&diag[local_end_y], 1, MPI_DOUBLE, next_process, 2, comm2D);
            }
        }
    }

    MPI_Barrier(comm2D);

    for (int i = local_start_z; i < local_end_z; i++)
    {
        for (int j = local_start_x; j < local_end_x; j++)
        {

            if (coord[direction] != dims[direction] - 1)
                MPI_Recv(&next, 1, MPI_DOUBLE, next_process, 1, comm2D, MPI_STATUS_IGNORE);
            else
                local_sol[j + local_end_y * local_size_x + i * local_size_x * local_size_y] /= diag[local_end_y];

            for (int k = local_end_y; k >= local_start_y; k--)
            {
                if (k != local_end_y)
                    local_sol[j + k * local_size_x + i * local_size_x * local_size_y] -= (exdiag * next) / diag[k];
                next = local_sol[j + k * local_size_x + i * local_size_x * local_size_y];
            }

            if (coord[direction] != 0)
                MPI_Send(&next, 1, MPI_DOUBLE, prev_process, 1, comm2D);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void thomas_z(const double exdiag, const int rank, const int size, const int start, const int end, const int direction)
{

    double prev;
    double next;
    std::vector<double> diag;
    double w;

    int source, next_process, prev_process;

    diag.resize(local_size_z);
    diag[local_start_x] = 1.0 + dt / dx2;

    MPI_Cart_shift(comm2D, direction, 1, &source, &next_process);
    MPI_Cart_shift(comm2D, direction, -1, &source, &prev_process);

    for (int j = local_start_y; j < local_end_y; j++)
    {
        for (int i = local_start_x; i < local_end_x; i++)
        {

            if (coord[direction] != 0)
            {
                MPI_Recv(&prev, 1, MPI_DOUBLE, prev_process, 1, comm2D, MPI_STATUS_IGNORE);
                MPI_Recv(&diag[local_start_z], 1, MPI_DOUBLE, prev_process, 2, comm2D, MPI_STATUS_IGNORE);
            }

            for (int k = local_start_z; k <= local_end_z; k++)
            {
                if (k != local_start_z)
                    w = exdiag / diag[k - 1];
                diag[k] = diag[k] - w * exdiag;
                local_sol[i + j * local_size_x + k * local_size_x * local_size_y] -= w * prev;
                prev = local_sol[i + j * local_size_x + k * local_size_x * local_size_y];
            }

            if (coord[direction] != size - 1)
            {
                MPI_Send(&prev, 1, MPI_DOUBLE, next_process, 1, comm2D);
                MPI_Send(&diag[local_end_z], 1, MPI_DOUBLE, next_process, 2, comm2D);
            }
        }
    }

    MPI_Barrier(comm2D);

    for (int j = local_start_y; j < local_end_y; j++)
    {
        for (int i = local_start_x; i < local_end_x; i++)
        {

            if (coord[direction] != dims[direction] - 1)
                MPI_Recv(&next, 1, MPI_DOUBLE, next_process, 1, comm2D, MPI_STATUS_IGNORE);
            else
                local_sol[i + j * local_size_x + local_end_z * local_size_x * local_size_y] /= diag[local_end_z];

            for (int k = local_end_z; k >= local_start_z; k--)
            {
                if (k != local_end_z)
                    local_sol[i + j * local_size_x + k * local_size_x * local_size_y] -= (exdiag * next) / diag[k];
                next = local_sol[i + j * local_size_x + k * local_size_x * local_size_y];
            }

            if (coord[direction] != 0)
                MPI_Send(&next, 1, MPI_DOUBLE, prev_process, 1, comm2D);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

//----------------------------------------PARALLEL RHS----------------------------------------//
void RHS(const std::vector<int> &counts, const std::vector<int> &displs, const int start, const int end, const double dt, const double dx2, const int x)
{

    const double val_diag = 1.0 + dt / dx2;
    const double exdiag = -dt / 2.0 / dx2;
    const double coeff = dt / dx2;
    const double force = 3 * std::sin(dt * (x + 0.5)) + std::cos(dt * (x + 0.5));
    const double sinx = std::sin(X);

    for (int i = local_start_z; i <= local_end_z; i++)
        for (int j = local_start_y; j <= local_end_y; j++)
            for (int k = local_start_x; k <= local_end_x; k++)
            {

                const int local_pos = k + j * local_size_x + i * local_size_x * local_size_y;

                local_sol[local_pos] = dt * my_sin[coord[2] * N / dims[2] + i - 1] * my_sin[coord[1] * N / dims[1] + j - 1] * my_sin[coord[0] * N / dims[0] + k - 1] * force - 6 * coeff * solution[local_pos];

                if (k != local_start_x)
                    local_sol[local_pos] += coeff * solution[local_pos - 1];

                if (k == local_end_x)
                    local_sol[local_pos] += coeff * my_sin[coord[2] * N / dims[2] + i - 1] * my_sin[coord[1] * N / dims[1] + j - 1] * sinx * my_sin_time[x];
                else
                    local_sol[local_pos] += coeff * solution[local_pos + 1];

                if (j != local_start_y)
                    local_sol[local_pos] += coeff * solution[local_pos - local_size_x];

                if (j == local_end_y)
                    local_sol[local_pos] += coeff * my_sin[coord[2] * N / dims[2] + i - 1] * sinx * my_sin[coord[0] * N / dims[0] + k - 1] * my_sin_time[x];
                else
                    local_sol[local_pos] += coeff * solution[local_pos + local_size_x];

                if (i != local_start_z)
                    local_sol[local_pos] += coeff * solution[local_pos - local_size_x * local_size_y];

                if (i == local_end_z)
                    local_sol[local_pos] += coeff * sinx * my_sin[coord[1] * N / dims[1] + j - 1] * my_sin[coord[0] * N / dims[0] + k - 1] * my_sin_time[x];
                else
                    local_sol[local_pos] += coeff * solution[local_pos + local_size_x * local_size_y];
            }

    MPI_Barrier(MPI_COMM_WORLD);
}

//------------------------------------PARALLEL X DIRECTION-------------------------------------//
void parallel_x_direction(const std::vector<int> &counts, const std::vector<int> &displs, const int start, const int end, const double dt, const double dx2, const int x)
{

    const double val_diag = 1.0 + dt / dx2;
    const double exdiag = -dt / 2.0 / dx2;

    int size;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // boundary conditions
    if (coord[0] == dims[0] - 1)
    {
        for (int i = local_start_z; i <= local_end_z; i++)
            for (int j = local_start_y; j <= local_end_y; j++)
            {
                const int local_pos = (local_end_x) + j * local_size_x + i * local_size_x * local_size_y;
                local_sol[local_pos] -= exdiag * my_sin[coord[2] * N / dims[2] + i - 1] * std::sin(X) * my_sin[coord[1] * N / dims[1] + j - 1] * (my_sin_time[x + 1] - my_sin_time[x]);
            }
    }


    thomas_x(exdiag, rank, size, start, end, 0);
}

//------------------------------------PARALLEL Y DIRECTION-------------------------------------//
void parallel_y_direction(const std::vector<int> &counts, const std::vector<int> &displs, const int start, const int end, const double dt, const double dx2, const int x)
{

    const double val_diag = 1.0 + dt / dx2;
    const double exdiag = -dt / 2.0 / dx2;

    int size;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // boundary conditions
    if (coord[1] == dims[1] - 1)
    {
        for (int i = local_start_z; i <= local_end_z; i++)
            for (int k = local_start_x; k <= local_end_x; k++)
            {
                const int local_pos = k + local_end_y * local_size_x + i * local_size_x * local_size_y;
                local_sol[local_pos] -= exdiag * my_sin[coord[2] * N / dims[2] + i - 1] * std::sin(X) * my_sin[coord[0] * N / dims[0] + k - 1] * (my_sin_time[x + 1] - my_sin_time[x]);
            }
    }


    // solve Thomas algorithm
    // std::cout << "execute thomas" << std::endl;
    thomas_y(exdiag, rank, size, start, end, 1);
}

//------------------------------------PARALLEL Z DIRECTION-------------------------------------//
void parallel_z_direction(std::vector<int> &counts, const std::vector<int> &displs, const int start, const int end, const double dt, const double dx2, const int x)
{

    const double val_diag = 1.0 + dt / dx2;
    const double exdiag = -dt / 2.0 / dx2;

    int size;
    int rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // boundary conditions
    if (coord[2] == dims[2] - 1)
    {
        for (int j = local_start_y; j <= local_end_y; j++)
            for (int k = local_start_x; k <= local_end_x; k++)
            {
                const int local_pos = k + j * local_size_x + local_end_z * local_size_x * local_size_y;
                local_sol[local_pos] -= exdiag * my_sin[coord[1] * N / dims[1] + j - 1] * std::sin(X) * my_sin[coord[0] * N / dims[0] + k - 1] * (my_sin_time[x + 1] - my_sin_time[x]);
            }
    }

    // std::cout << "execute thomas on z" << std::endl;
    thomas_z(exdiag, rank, size, start, end, 2);
}

//------------------------------------------FINALIZE-------------------------------------------//
double finalize(const int x)
{

    for (int k = local_start_z; k <= local_end_z; k++)
        for (int j = local_start_y; j <= local_end_y; j++)
            for (int i = local_start_x; i <= local_end_x; i++)
                solution[i + j * local_size_x + k * local_size_x * local_size_y] += local_sol[i + j * local_size_x + k * local_size_x * local_size_y];

    MPI_Barrier(comm2D);

    // error
    if (x == nt - 1)
    {
        double error = 0.0;
        for (int k = local_start_z; k <= local_end_z; k++)
            for (int j = local_start_y; j <= local_end_y; j++)
                for (int i = local_start_x; i <= local_end_x; i++)
                    error += solution[i + j * local_size_x + k * local_size_x * local_size_y] - my_sin[coord[2] * N / dims[2] + k - 1] * my_sin[coord[1] * N / dims[1] + j - 1] * my_sin[coord[0] * N / dims[0] + i - 1] * my_sin_time[x + 1];
        
        // for (int k = local_start_z; k <= local_end_z; k++)
        //     for (int j = local_start_y; j <= local_end_y; j++)
        //         for (int i = local_start_x; i <= local_end_x; i++)
        //             if(solution[i + j * local_size_x + k * local_size_x * local_size_y] != solution[i + j * local_size_x + k * local_size_x * local_size_y])
        //                 std::cout << "Solution is NaN" << std::endl;
        
        // std::cout << "Error: " << error << std::endl;
        return error;
    }

    return 0.0;
}

//--------------------------------------------MAIN---------------------------------------------//
int main(int argc, char **argv)
{

    // i = asse z (piano)
    // j = asse y (riga)
    // k = asse x (colonna)

    auto start_time = std::chrono::high_resolution_clock::now();

    MPI_Init(&argc, &argv);

    int size, rank;
    int start_x, start_y, start_z;
    int end_x, end_y, end_z;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> starts(size);
    std::vector<int> ends(size);
    std::vector<int> counts(size);
    std::vector<int> displs(size);

    dx = X / (N + 1);
    dx2 = dx * dx;
    dt = T / nt;

    double error = 0.0;

    // Create cart with mpi
    dims[0] = dims[1] = dims[2] = 0;

    /* create cartesian topology for processes */
    MPI_Dims_create(size, dimensions, dims);
    std::cout << "PW[" << rank << "]/[" << size << "]: PEdims = [" << dims[0] << " x " << dims[1] << " x " << dims[2] << "]" << std::endl;
    /* create cartesian mapping */
    wrap_around[0] = wrap_around[1] = wrap_around[2] = 0; // set periodicity
    reorder = 1;
    int ierr = 0;
    ierr = MPI_Cart_create(MPI_COMM_WORLD, dimensions, dims,
                           wrap_around, reorder, &comm2D);
    if (ierr != 0)
        std::cout << "ERROR[" << ierr << "] creating CART" << std::endl;

    /* find my coordinates in the cartesian communicator group */
    MPI_Cart_coords(comm2D, rank, dimensions, coord);
    std::cout << "PW[" << rank << "], my coords = (" << coord[0] << "," << coord[1] << "," << coord[2] << ")" << std::endl;

    local_size_x = N / dims[0] + 2;
    local_size_y = N / dims[1] + 2;
    local_size_z = N / dims[2] + 2;

    local_sol.resize((local_size_x) * (local_size_y) * (local_size_z));
    solution.resize((local_size_x) * (local_size_y) * (local_size_z));

    local_end_x = local_size_x - 1;
    local_end_y = local_size_y - 1;
    local_end_z = local_size_z - 1;

    if (coord[0] == 0)
        local_start_x = 1;
    if (coord[0] == dims[0] - 1)
        local_end_x = local_size_x - 2;
    if (coord[1] == 0)
        local_start_y = 1;
    if (coord[1] == dims[1] - 1)
        local_end_y = local_size_y - 2;
    if (coord[2] == 0)
        local_start_z = 1;
    if (coord[2] == dims[2] - 1)
        local_end_z = local_size_z - 2;


    for (int i = 1; i <= N; i++)
        my_sin[i - 1] = std::sin(i * dx);

    for (int i = 0; i <= nt; i++)
        my_sin_time[i] = std::sin(i * dt);

    for (int i = 0; i < local_size_x * local_size_y * local_size_z; i++)
        solution[i] = 0.0;

    MPI_Barrier(comm2D);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    //--------------------------------TIME LOOP------------------------------------//

    for (int x = 0; x < nt; x++)
    {


        //-------------------------------RHS-----------------------------------//


        //  Bcast solution - not optimal but it works
        // std::cout << "RHS" << std::endl;
        RHS(counts, displs, start_x, end_x, dt, dx2, x);
        MPI_Barrier(comm2D);

        //---------------------------X DIRECTION-------------------------------//

        // std::cout << "Thomas on x direction" << std::endl;

        parallel_x_direction(counts, displs, start_x, end_x, dt, dx2, x);
        if(x == 0)
            for(int i = 0; i < local_size_x * local_size_y * local_size_z; i++)
                std::cout << local_sol[i] << std::endl;
        MPI_Barrier(comm2D);

        //---------------------------Y DIRECTION-------------------------------//

        // std::cout << "Thomas on y direction" << std::endl;

        parallel_y_direction(counts, displs, start_y, end_y, dt, dx2, x); // TO DO
        MPI_Barrier(comm2D);

        //---------------------------Z DIRECTION-------------------------------//

        // std::cout << "Thomas on z direction" << std::endl;

        parallel_z_direction(counts, displs, start_z, end_z, dt, dx2, x); // TO DO
        MPI_Barrier(comm2D);            

        //----------------------------FINALIZE---------------------------------//

        // std::cout << "Finalize" << std::endl;
        error = finalize(x);

        // if(solution[1] != solution[1])
        //     std::cout << "Solution is NaN" << std::endl;

    }

    double error_tot;
    MPI_Reduce(&error, &error_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    std::cout << "Error of process " << rank << " is: " << std::sqrt(error) / (double)(local_size_x * local_size_y * local_size_z) << std::endl;

    if (rank == 0)
    {
        std::cout << "Error with " << N << " x " << N << " x " << N
                  << " spatial elements is: " << std::sqrt(error_tot) / (double)(N * N * N) << std::endl;
        std::cout << "Time pre-loop = " << duration.count() / pow(10, 6) << std::endl;
    }

    MPI_Finalize();

    return 0;
}