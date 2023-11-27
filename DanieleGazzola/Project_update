#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include <array>
#include <string>
#include <mpi.h>

//------------------------------------------VARIABLES-----------------------------------------//
#define X  1.0
#define T  1.0
#define N  100
#define nt 100

std::vector<double>         solution;
std::vector<double>         delta_solution;

std::array<double, N>       diag;
std::array<double, N + 2>   my_sin;
std::array<double, nt + 1>  my_sin_time;

MPI_Comm comm3D;
int dims[3];
int coord[3];

//--------------------------------------THOMAS ALGORITHM--------------------------------------//
void thomas(const double exdiag, const int batch0, const int batch1, const int batch2, const int n){

    int batch;
    double prev;
    double next;

    int size0 = batch0 + 2;
    int size1 = batch1 + 2;
    int size2 = batch2 + 2;

    int rank_source, rank_prev, rank_next;

    MPI_Cart_shift(comm3D, n, -1, &rank_source, &rank_prev);
    MPI_Cart_shift(comm3D, n, 1, &rank_source, &rank_next);

    int start = 0;
    if(coord[n] == 0)
        start++;

    for (int i = 1; i <= batch2; i++){
        for (int j = 1; j <= batch1; j++){

            batch = i * size1 * size0 + j * size0;

            if (coord[n] != 0)
                MPI_Recv(&delta_solution[batch], 1, MPI_DOUBLE, rank_prev, batch, comm3D, MPI_STATUS_IGNORE);

            for (int k = start; k < batch0; k++)
                delta_solution[batch + k + 1] -= exdiag / diag[batch0 * coord[n] + k - 1] * delta_solution[batch + k];

            if (coord[n] != dims[n] - 1)
                MPI_Send(&delta_solution[batch + batch0], 1, MPI_DOUBLE, rank_next, batch, comm3D);
        }
    }

    MPI_Barrier(comm3D);

    int end = batch0 - 1;
    if(coord[n] == dims[n] - 1)
        end--;

    for (int i = 1; i <= batch2; i++){
        for (int j = 1; j <= batch1; j++){

            batch = i * size1 * size0 + j * size0;

            if (coord[n] != dims[n] - 1)
                MPI_Recv(&delta_solution[batch + batch0 + 1], 1, MPI_DOUBLE, rank_next, batch, comm3D, MPI_STATUS_IGNORE);
            else
                delta_solution[batch + batch0] = delta_solution[batch + batch0] / diag[N - 1];

            for (int k = end; k >= 0; k--)
                delta_solution[batch + k + 1] = (delta_solution[batch + k + 1] - exdiag * delta_solution[batch + k + 2]) / diag[batch0 * coord[n] + k];

            if (coord[n] != 0)
                MPI_Send(&delta_solution[batch + 1], 1, MPI_DOUBLE, rank_prev, batch, comm3D);
        }
    }

    MPI_Barrier(comm3D);
}

//----------------------------------------PARALLEL RHS----------------------------------------//
void RHS(const int batch_x, const int batch_y, const int batch_z,  const double dt, const double dx2, const int x){

    const double coeff = dt / dx2;
    const double force = 3 * std::sin(dt * (x + 0.5)) + std::cos(dt * (x + 0.5));

    const int size_x = batch_x + 2;
    const int size_y = batch_y + 2;
    const int size_z = batch_z + 2;

    if (coord[0] == dims[0] - 1)
        for (int i = 1; i <= batch_z; i++)
            for (int j = 1; j <= batch_y; j++)
                solution[i * size_y * size_x + j * size_x + (batch_x + 1)] = my_sin[batch_z * coord[2] + i] * my_sin[batch_y * coord[1] + j] * my_sin[N+1] * my_sin_time[x];

    if (coord[1] == dims[1] - 1)
        for (int i = 1; i <= batch_z; i++)
            for (int k = 1; k <= batch_x; k++)
                solution[i * size_y * size_x + (batch_y + 1) * size_x + k] = my_sin[batch_z * coord[2] + i] * my_sin[N+1] * my_sin[batch_x * coord[0] + k] * my_sin_time[x];

    if (coord[2] == dims[2] - 1)
        for (int j = 1; j <= batch_y; j++)
            for (int k = 1; k <= batch_x; k++)
                solution[(batch_z + 1) * size_y * size_x + j * size_x + k] = my_sin[N+1] * my_sin[batch_y * coord[1] + j] * my_sin[batch_x * coord[0] + k] * my_sin_time[x];

    if (coord[0] == 0)
        for (int i = 1; i <= batch_z; i++)
            for (int j = 1; j <= batch_y; j++)
                solution[i * size_y * size_x + j * size_x + 0] = 0;

    if (coord[1] == 0)
        for (int i = 1; i <= batch_z; i++)
            for (int k = 1; k <= batch_x; k++)
                solution[i * size_y * size_x + 0 * size_x + k] = 0;

    if (coord[2] == 0)
        for (int j = 1; j <= batch_y; j++)
            for (int k = 1; k <= batch_x; k++)
                solution[0 * size_y * size_x + j * size_x + k] = 0;

    // internal points
    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_y; j++)
            for (int k = 1; k <= batch_x; k++){

                int pos = i * size_y * size_x + j * size_x + k;

                delta_solution[pos]  = dt * my_sin[batch_z * coord[2] + i] * my_sin[batch_y * coord[1] + j] * my_sin[batch_x * coord[0] + k] * force;
                delta_solution[pos] -= 6 * coeff * solution[pos];
                delta_solution[pos] += coeff * (  solution[pos - 1] + solution[pos + 1]
                                                + solution[pos - size_x] + solution[pos + size_x]
                                                + solution[pos - size_y * size_x] + solution[pos + size_y * size_x]);
            }

    MPI_Barrier(comm3D);
}

//------------------------------------PARALLEL X DIRECTION------------------------------------//
void parallel_x_direction(const int batch_x, const int batch_y, const int batch_z, const double dt, const double dx2, const int x){

    const double exdiag = -dt / 2.0 / dx2;

    const int size_x = batch_x + 2;
    const int size_y = batch_y + 2;
    const int size_z = batch_z + 2;

    //boundary conditions
    if (coord[0] == dims[0] - 1)
        for (int i = 1; i <= batch_z; i++)
            for (int j = 1; j <= batch_y; j++)
                delta_solution[i * size_y * size_x + j * size_x + batch_x] -= exdiag * my_sin[batch_z * coord[2] + i] * my_sin[batch_y * coord[1] + j] * my_sin[N+1] * (my_sin_time[x + 1] - my_sin_time[x]);

    //solve Thomas algorithm
    thomas(exdiag, batch_x, batch_y, batch_z, 0);
}

//------------------------------------PARALLEL Y DIRECTION------------------------------------//
void parallel_y_direction(const int batch_x, const int batch_y, const int batch_z, const double dt, const double dx2, const int x){

    const double exdiag = -dt / 2.0 / dx2;

    const int size_x = batch_x + 2;
    const int size_y = batch_y + 2;
    const int size_z = batch_z + 2;

    //boundary conditions
    if (coord[1] == dims[1] - 1)
        for (int k = 1; k <= batch_x; k++)
            for (int i = 1; i <= batch_z; i++)
                delta_solution[k * size_z * size_y + i * size_y + batch_y] -= exdiag * my_sin[batch_z * coord[2] + i] * my_sin[N+1] * my_sin[batch_x * coord[0] + k] * (my_sin_time[x + 1] - my_sin_time[x]);

    //solve Thomas algorithm
    thomas(exdiag, batch_y, batch_z, batch_x, 1);
}

//------------------------------------PARALLEL Z DIRECTION------------------------------------//
void parallel_z_direction(const int batch_x, const int batch_y, const int batch_z, const double dt, const double dx2, const int x){

    const double exdiag = -dt / 2.0 / dx2;

    const int size_x = batch_x + 2;
    const int size_y = batch_y + 2;
    const int size_z = batch_z + 2;

    //boundary conditions
    if (coord[2] == dims[2] - 1)
        for (int j = 1; j <= batch_y; j++)
            for (int k = 1; k <= batch_x; k++)
                delta_solution[j * size_x * size_z + k * size_z + batch_z] -= exdiag * my_sin[N+1] * my_sin[batch_y * coord[1] + j] * my_sin[batch_x * coord[0] + k] * (my_sin_time[x + 1] - my_sin_time[x]);

    //solve Thomas algorithm
    thomas(exdiag, batch_z, batch_x, batch_y, 2);
}

//------------------------------------------FINALIZE------------------------------------------//
double finalize(const int batch_x, const int batch_y, const int batch_z, const int x){

    const int size_x = batch_x + 2;
    const int size_y = batch_y + 2;
    const int size_z = batch_z + 2;

    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_y; j++)
            for (int k = 1; k <= batch_x; k++)
                solution[i * size_y * size_x + j * size_x + k] += delta_solution[i * size_y * size_x + j * size_x + k];

    //error
    if (x == nt - 1){
        double error = 0.0;
        double errort = 0.0;
        for (int i = 1; i <= batch_z; i++)
            for (int j = 1; j <= batch_y; j++)
                for (int k = 1; k <= batch_x; k++)
                    error += std::pow(solution[i * size_y * size_x + j * size_x + k] - my_sin[batch_z * coord[2] + i] * my_sin[batch_y * coord[1] + j] * my_sin[batch_x * coord[0] + k] * my_sin_time[x + 1], 2);
        MPI_Reduce(&error, &errort, 1, MPI_DOUBLE, MPI_SUM, 0, comm3D);
        errort = std::sqrt(errort) / std::pow(N, 3);
        return errort;
    }

    return 0.0;
}

//-------------------------------------------ROTATE-------------------------------------------//
void rotate(const std::array<int, 3>& batch, const int n){

    const int size0 = batch[(n+0) % 3] + 2;
    const int size1 = batch[(n+1) % 3] + 2;
    const int size2 = batch[(n+2) % 3] + 2;

    const int local_size = size0 * size1 * size2;

    std::vector<double> temp(local_size);

    for (int i = 0; i < local_size; i++)
        temp[i] = delta_solution[i];

    for (int i = 1; i <= batch[(n+0) % 3]; i++)
        for (int j = 1; j <= batch[(n+2) % 3]; j++)
            for (int k = 1; k <= batch[(n+1) % 3]; k++)
                delta_solution[i * size2 * size1 + j * size1 + k] = temp[j * size1 * size0 + k * size0 + i];
}

//----------------------------------------COMUNICATOR-----------------------------------------//
void comunicator(const int batch_x, const int batch_y, const int batch_z){

    const int size_x = batch_x + 2;
    const int size_y = batch_y + 2;
    const int size_z = batch_z + 2;

    int rank_source, rank_prev, rank_next;

    MPI_Cart_shift(comm3D, 0, -1, &rank_source, &rank_prev);
    MPI_Cart_shift(comm3D, 0, 1, &rank_source, &rank_next);

    std::vector<double> send(batch_y * batch_z);
    std::vector<double> recv(batch_y * batch_z);

    //send to prev on x axis

    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_y; j++)
            send[(i - 1) * batch_y + (j - 1)] = solution[i * size_y * size_x + j * size_x + 1];

    MPI_Sendrecv(send.data(), batch_y * batch_z, MPI_DOUBLE, rank_prev, 0, recv.data(), batch_y * batch_z, MPI_DOUBLE, rank_next, 0, comm3D, MPI_STATUS_IGNORE);

    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_y; j++)
            solution[i * size_y * size_x + j * size_x + batch_x + 1] = recv[(i - 1) * batch_y + (j - 1)];

    MPI_Barrier(comm3D);

    // send to next on x axis

    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_y; j++)
            send[(i - 1) * batch_y + (j - 1)] = solution[i * size_y * size_x + j * size_x + batch_x];

    MPI_Sendrecv(send.data(), batch_y * batch_z, MPI_DOUBLE, rank_next, 0, recv.data(), batch_y * batch_z, MPI_DOUBLE, rank_prev, 0, comm3D, MPI_STATUS_IGNORE);

    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_y; j++)
            solution[i * size_y * size_x + j * size_x + 0] = recv[(i - 1) * batch_y + (j - 1)];

    MPI_Barrier(comm3D);

    MPI_Cart_shift(comm3D, 1, -1, &rank_source, &rank_prev);
    MPI_Cart_shift(comm3D, 1, 1, &rank_source, &rank_next);

    send.resize(batch_z * batch_x);
    recv.resize(batch_z * batch_x);

    //send to prev on y axis

    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_x; j++)
            send[(i - 1) * batch_x + (j - 1)] = solution[i * size_y * size_x + 1 * size_x + j];

    MPI_Sendrecv(send.data(), batch_z * batch_x, MPI_DOUBLE, rank_prev, 0, recv.data(), batch_z * batch_x, MPI_DOUBLE, rank_next, 0, comm3D, MPI_STATUS_IGNORE);

    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_x; j++)
            solution[i * size_y * size_x + (batch_y + 1) * size_x + j] = recv[(i - 1) * batch_x + (j - 1)];

    MPI_Barrier(comm3D);

    // send to next on y axis

    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_x; j++)
            send[(i - 1) * batch_x + (j - 1)] = solution[i * size_y * size_x + batch_y * size_x + j];

    MPI_Sendrecv(send.data(), batch_z * batch_x, MPI_DOUBLE, rank_next, 0, recv.data(), batch_z * batch_x, MPI_DOUBLE, rank_prev, 0, comm3D, MPI_STATUS_IGNORE);

    for (int i = 1; i <= batch_z; i++)
        for (int j = 1; j <= batch_x; j++)
            solution[i * size_y * size_x + 0 * size_x + j] = recv[(i - 1) * batch_x + (j - 1)];

    MPI_Barrier(comm3D);

    MPI_Cart_shift(comm3D, 2, -1, &rank_source, &rank_prev);
    MPI_Cart_shift(comm3D, 2, 1, &rank_source, &rank_next);

    send.resize(batch_x * batch_y);
    recv.resize(batch_x * batch_y);

    //send to prev on z axis

    for (int i = 1; i <= batch_y; i++)
        for (int j = 1; j <= batch_x; j++)
            send[(i - 1) * batch_x + (j - 1)] = solution[1 * size_y * size_x + i * size_x + j];

    MPI_Sendrecv(send.data(), batch_x * batch_y, MPI_DOUBLE, rank_prev, 0, recv.data(), batch_x * batch_y, MPI_DOUBLE, rank_next, 0, comm3D, MPI_STATUS_IGNORE);

    for (int i = 1; i <= batch_y; i++)
        for (int j = 1; j <= batch_x; j++)
            solution[(batch_z + 1) * size_y * size_x + i * size_x + j] = recv[(i - 1) * batch_x + (j - 1)];

    MPI_Barrier(comm3D);

    // send to next on z axis

    for (int i = 1; i <= batch_y; i++)
        for (int j = 1; j <= batch_x; j++)
            send[(i - 1) * batch_x + (j - 1)] = solution[batch_z * size_y * size_x + i * size_x + j];

    MPI_Sendrecv(send.data(), batch_x * batch_y, MPI_DOUBLE, rank_next, 0, recv.data(), batch_x * batch_y, MPI_DOUBLE, rank_prev, 0, comm3D, MPI_STATUS_IGNORE);

    for (int i = 1; i <= batch_y; i++)
        for (int j = 1; j <= batch_x; j++)
            solution[0 * size_y * size_x + i * size_x + j] = recv[(i - 1) * batch_x + (j - 1)];

    MPI_Barrier(comm3D);
}

//--------------------------------------------MAIN--------------------------------------------//
int main(int argc, char** argv){

    // i = asse z (piano)   (comm[2])
    // j = asse y (riga)    (comm[1])
    // k = asse x (colonna) (comm[0])

    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int start_x, start_y, start_z;
    int end_x, end_y, end_z;

    const double dx       = X / (N + 1);
    const double dx2      = dx * dx;
    const double dt       = T / nt;
    const double val_diag = 1.0 + dt / dx2;
    const double exdiag   = -dt / 2.0 / dx2;
          double error    = 0.0;

    const int dimensions  = 3;
    const int reorder     = 1;

    std::array<int, 3> wrap;

    dims[0] = dims[1] = dims[2] = 0;
    MPI_Dims_create(size, dimensions, dims);

    wrap[0] = wrap[1] = wrap[2] = 0;
    MPI_Cart_create(MPI_COMM_WORLD, dimensions, dims, wrap.data(), reorder, &comm3D);
    MPI_Cart_coords(comm3D, rank, dimensions, coord);

    const int batch_x = N / dims[0];
    const int batch_y = N / dims[1];
    const int batch_z = N / dims[2];

    std::array<int, 3> batchs;
    batchs[0] = batch_x;
    batchs[1] = batch_y;
    batchs[2] = batch_z;

    const int local_size = (batch_x + 2) * (batch_y + 2) * (batch_z + 2);

    solution.resize(local_size);
    delta_solution.resize(local_size);

    for (int i = 0; i < N + 2; i++)
        my_sin[i] = std::sin(i * dx);

    for (int i = 0; i <= nt; i++)
        my_sin_time[i] = std::sin(i * dt);

    for (int i = 0; i < local_size; i++)
        delta_solution[i] = 0.0;

    for (int i = 0; i < local_size; i++)
        solution[i] = 0.0;

    diag[0] = val_diag;
    for (int k = 1; k < N; k++)
        diag[k] = val_diag - exdiag / diag[k-1] * exdiag;

    MPI_Barrier(comm3D);

//--------------------------------TIME LOOP------------------------------------//

    for (int x = 0; x < nt; x++){

        //-------------------------------RHS-----------------------------------//

        RHS(batch_x, batch_y, batch_z, dt, dx2, x);

        //---------------------------X DIRECTION-------------------------------//

        parallel_x_direction(batch_x, batch_y, batch_z, dt, dx2, x);
        rotate(batchs, 0);

        //---------------------------Y DIRECTION-------------------------------//

        parallel_y_direction(batch_x, batch_y, batch_z, dt, dx2, x);
        rotate(batchs, 1);

        //---------------------------Z DIRECTION-------------------------------//

        parallel_z_direction(batch_x, batch_y, batch_z, dt, dx2, x);
        rotate(batchs, 2);

        //----------------------------FINALIZE---------------------------------//

        error = finalize(batch_x, batch_y, batch_z, x);
        //PENSO OK FIN QUA
        comunicator(batch_x, batch_y, batch_z);
    }

    if (rank == 0)
        std::cout << "Error: " << error << std::endl;

    MPI_Finalize();

    return 0;
}