#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include <array>
#include <string>
#include <mpi.h>
#include <omp.h>

//------------------------------------------VARIABLES-----------------------------------------//
#define X  1.0
#define T  1.0
#define N  100
#define nt 100

std::array<double, N * N * N>   delta_solution;
std::array<double, N * N * N>   solution;
std::array<double, N>           diag;
std::array<double, N>           my_sin;
std::array<double, nt + 1>      my_sin_time;

//------------------------------------SEQUENTIAL PROGRAMMA------------------------------------//

//----------------------------------------X  DIRECTION----------------------------------------//
void x_direction(const double dt, const double dx2, const int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    //boundary conditions
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
                delta_solution[i * N * N + j * N + N - 1] -= exdiag * my_sin[i] * my_sin[j] * std::sin(X) * (my_sin_time[x + 1] - my_sin_time[x]);

    diag[0] = val_diag;
    for (int k = 1; k < N; k++)
        diag[k] = val_diag - exdiag / diag[k-1] * exdiag;

    //solve Thomas algorithm
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            for (int k = 1; k < N; k++)
                delta_solution[i * N * N + j * N + k] -= exdiag / diag[k-1] * delta_solution[i * N * N + j * N + (k - 1)];

            delta_solution[i * N * N + j * N + (N - 1)] = delta_solution[i * N * N + j * N + (N - 1)] / diag[N - 1];
            for (int k = N-2; k >= 0; k--)
                delta_solution[i * N * N + j * N + k] = (delta_solution[i * N * N + j * N + k] - exdiag * delta_solution[i * N * N + j * N + (k + 1)]) / diag[k];
        }
    }
}

//----------------------------------------Y  DIRECTION----------------------------------------//
void y_direction(const double dt, const double dx2, const int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    //boundary conditions
    for (int i = 0; i < N; i++)
        for (int k = 0; k < N; k++)
                delta_solution[i * N * N + (N - 1) * N + k] -= exdiag * my_sin[i] * std::sin(X) * my_sin[k] * (my_sin_time[x + 1] - my_sin_time[x]);

    diag[0] = val_diag;
    for (int k = 1; k < N; k++)
        diag[k] = val_diag - exdiag / diag[k-1] * exdiag;

    //solve Thomas algorithm
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            for (int k = 1; k < N; k++)
                delta_solution[i * N * N + k * N + j] -= exdiag / diag[k-1] * delta_solution[i * N * N + (k - 1) * N + j];

            delta_solution[i * N * N + (N - 1) * N + j] = delta_solution[i * N * N + (N - 1) * N + j] / diag[N - 1];
            for (int k = N-2; k >= 0; k--)
                delta_solution[i * N * N + k * N + j] = (delta_solution[i * N * N + k * N + j] - exdiag * delta_solution[i * N * N + (k + 1) * N + j]) / diag[k];
        }
    }
}

//----------------------------------------Z  DIRECTION----------------------------------------//
void z_direction(const double dt, const double dx2, const int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    //boundary conditions
    for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++)
            delta_solution[(N -1) * N * N + j * N + k] -= exdiag * std::sin(X) * my_sin[j] * my_sin[k] * (my_sin_time[x + 1] - my_sin_time[x]);

    diag[0] = val_diag;
    for (int k = 1; k < N; k++)
        diag[k] = val_diag - exdiag / diag[k-1] * exdiag;

    //solve Thomas algorithm
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            for (int k = 1; k < N; k++)
                delta_solution[k * N * N + i * N + j] -= exdiag / diag[k-1] * delta_solution[(k - 1) * N * N + i * N + j];

            delta_solution[(N - 1) * N * N + i * N + j] = delta_solution[(N - 1) * N * N + i * N + j] / diag[N - 1];
            for (int k = N-2; k >= 0; k--)
                delta_solution[k * N * N + i * N + j] = (delta_solution[k * N * N + i * N + j] - exdiag * delta_solution[(k + 1) * N * N + i * N + j]) / diag[k];
        }
    }
}

//-------------------------------------PARALLEL PROGRAMMA-------------------------------------//

//--------------------------------------THOMAS ALGORITHM--------------------------------------//
void thomas(std::vector<double>& local_sol, const double exdiag, const int rank, const int size, const int start, const int end){

    int batch;
    double prev;
    double next;

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            batch = i * N + j;

            if (rank != 0)
                MPI_Recv(&prev, 1, MPI_DOUBLE, rank - 1, batch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int k = start; k < end; k++){
                if (k != 0)
                    local_sol[(k - start) * N * N + batch] -= exdiag / diag[k - 1] * prev;
                prev = local_sol[(k - start) * N * N + batch];
            }

            if (rank != size - 1)
                MPI_Send(&prev, 1, MPI_DOUBLE, rank + 1, batch, MPI_COMM_WORLD);
        }
    }

   // MPI_Barrier(MPI_COMM_WORLD);

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            batch = i * N + j;

            if (rank != size - 1)
                MPI_Recv(&next, 1, MPI_DOUBLE, rank + 1, batch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else
                local_sol[(end - start - 1) * N * N + batch] = local_sol[(end - start - 1) * N * N + batch] / diag[N - 1];

            for (int k = end - 1; k >= start; k--){
                if (k != N - 1)
                    local_sol[(k - start) * N * N + batch] = (local_sol[(k - start) * N * N + batch] - exdiag * next) / diag[k];
                next = local_sol[(k - start) * N * N + batch];
            }

            if (rank != 0)
                MPI_Send(&next, 1, MPI_DOUBLE, rank - 1, batch, MPI_COMM_WORLD);
        }
    }

   // MPI_Barrier(MPI_COMM_WORLD);
}

//----------------------------------------PARALLEL RHS----------------------------------------//
void RHS(const std::vector<int>& counts, const std::vector<int>& displs, const int start, const int end,  const double dt, const double dx2, const int x){

    const double val_diag = 1.0 + dt / dx2;
    const double exdiag =   -dt / 2.0 / dx2;
    const double coeff =    dt / dx2;
    const double force =    3 * std::sin(dt * (x + 0.5)) + std::cos(dt * (x + 0.5));
    const double sinx =     std::sin(X);

    std::vector<double> local_sol((end - start) * N * N);
    //#pragma omp parallel for num_threads(10)
    for (int i = start; i < end; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++){

                const int local_pos = (i - start) * N * N + j * N + k;
                const int global_pos = i * N * N + j * N + k;

                local_sol[local_pos] = dt * my_sin[i] * my_sin[j] * my_sin[k] * force - 6 * coeff * solution[global_pos];

                if (k != 0)
                    local_sol[local_pos] += coeff * solution[global_pos - 1];

                if (k == N - 1)
                    local_sol[local_pos] += coeff * my_sin[i] * my_sin[j] * sinx * my_sin_time[x];
                else
                    local_sol[local_pos] += coeff * solution[global_pos + 1];

                if (j != 0)
                    local_sol[local_pos] += coeff * solution[global_pos - N];

                if (j == N - 1)
                    local_sol[local_pos] += coeff * my_sin[i] * sinx * my_sin[k] * my_sin_time[x];
                else
                    local_sol[local_pos] += coeff * solution[global_pos + N];

                if (i != 0)
                    local_sol[local_pos] += coeff * solution[global_pos - N * N];

                if (i == N - 1)
                    local_sol[local_pos] += coeff * sinx * my_sin[j] * my_sin[k] * my_sin_time[x];
                else
                    local_sol[local_pos] += coeff * solution[global_pos + N * N];
            }

    MPI_Gatherv(local_sol.data(), ((end - start) * N * N), MPI_DOUBLE, delta_solution.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
}

//------------------------------------PARALLEL X DIRECTION-------------------------------------//
void parallel_x_direction(const std::vector<int>& counts, const std::vector<int>& displs, const int start, const int end, const double dt, const double dx2, const int x){

    const double val_diag = 1.0 + dt / dx2;
    const double exdiag = -dt / 2.0 / dx2;

    int size;
    int rank;

    std::array<double, N * N * N> global_sol;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //boundary conditions
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
                delta_solution[i * N * N + j * N + (N - 1)] -= exdiag * my_sin[i] * my_sin[j] * std::sin(X) * (my_sin_time[x + 1] - my_sin_time[x]);

    diag[0] = val_diag;
    for (int k = 1; k < N; k++)
        diag[k] = val_diag - exdiag / diag[k-1] * exdiag;

    //solve Thomas algorithm

    if (rank == 0)
    #pragma omp parallel for num_threads(10)
        for (int k = 0; k < N; k++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    global_sol[k * N * N + i * N + j] = delta_solution[i * N * N + j * N + k];

    std::vector<double> local_sol((end - start) * N * N);

    MPI_Scatterv(global_sol.data(), counts.data(), displs.data(), MPI_DOUBLE, local_sol.data(), ((end - start) * N * N), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    thomas(local_sol, exdiag, rank, size, start, end);

    MPI_Gatherv(local_sol.data(), ((end - start) * N * N), MPI_DOUBLE, global_sol.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    #pragma omp parallel for num_threads(10)
        for (int k = 0; k < N; k++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    delta_solution[i * N * N + j * N + k] = global_sol[k * N * N + i * N + j];
}

//------------------------------------PARALLEL Y DIRECTION-------------------------------------//
void parallel_y_direction(const std::vector<int>& counts, const std::vector<int>& displs, const int start, const int end, const double dt, const double dx2, const int x){

    const double val_diag = 1.0 + dt / dx2;
    const double exdiag = -dt / 2.0 / dx2;

    int size;
    int rank;

    std::array<double, N * N * N> global_sol;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //boundary conditions
    for (int i = 0; i < N; i++)
        for (int k = 0; k < N; k++)
                delta_solution[i * N * N + (N - 1) * N + k] -= exdiag * my_sin[i] * std::sin(X) * my_sin[k] * (my_sin_time[x + 1] - my_sin_time[x]);

    diag[0] = val_diag;
    for (int k = 1; k < N; k++)
        diag[k] = val_diag - exdiag / diag[k-1] * exdiag;

    //solve Thomas algorithm

    if (rank == 0)
    #pragma omp parallel for num_threads(10)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                for (int i = 0; i < N; i++)
                    global_sol[j * N * N + k * N + i] = delta_solution[i * N * N + j * N + k];

    std::vector<double> local_sol((end - start) * N * N);

    MPI_Scatterv(global_sol.data(), counts.data(), displs.data(), MPI_DOUBLE, local_sol.data(), ((end - start) * N * N), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    thomas(local_sol, exdiag, rank, size, start, end);

    MPI_Gatherv(local_sol.data(), ((end - start) * N * N), MPI_DOUBLE, global_sol.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    #pragma omp parallel for num_threads(10)
        for (int k = 0; k < N; k++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    delta_solution[i * N * N + j * N + k] = global_sol[j * N * N + k * N + i];
}

//------------------------------------PARALLEL Z DIRECTION-------------------------------------//
void parallel_z_direction(const std::vector<int>& counts, const std::vector<int>& displs, const int start, const int end, const double dt, const double dx2, const int x){

    const double val_diag = 1.0 + dt / dx2;
    const double exdiag = -dt / 2.0 / dx2;

    int size;
    int rank;

    std::array<double, N * N * N> global_sol;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //boundary conditions
    for (int j = 0; j < N; j++)
        for (int k = 0; k < N; k++)
                delta_solution[(N - 1) * N * N + j * N + k] -= exdiag * std::sin(X) * my_sin[j] * my_sin[k] * (my_sin_time[x + 1] - my_sin_time[x]);

    diag[0] = val_diag;
    for (int k = 1; k < N; k++)
        diag[k] = val_diag - exdiag / diag[k-1] * exdiag;

    //solve Thomas algorithm

    std::vector<double> local_sol((end - start) * N * N);

    MPI_Scatterv(delta_solution.data(), counts.data(), displs.data(), MPI_DOUBLE, local_sol.data(), ((end - start) * N * N), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    thomas(local_sol, exdiag, rank, size, start, end);

    MPI_Gatherv(local_sol.data(), ((end - start) * N * N), MPI_DOUBLE, delta_solution.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

//------------------------------------------FINALIZE-------------------------------------------//
double finalize(const int x){

    #pragma omp parallel for num_threads(10)
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                solution[i * N * N + j * N + k] += delta_solution[i * N * N + j * N + k];

    //error
    if (x == nt - 1){
        double error = 0.0;
        #pragma omp parallel for num_threads(10)
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++)
                    error += (solution[i * N * N + j * N + k] - my_sin[i] * my_sin[j] * my_sin[k] * my_sin_time[x + 1]) * (solution[i * N * N + j * N + k] - my_sin[i] * my_sin[j] * my_sin[k] * my_sin_time[x + 1]);
        error = std::sqrt(error) / std::pow(N, 3);
        return error;
    }

    return 0.0;
}

//--------------------------------------------MAIN---------------------------------------------//
int main(int argc, char** argv){

    // i = asse z (piano)
    // j = asse y (riga)
    // k = asse x (colonna)

auto start_time = std::chrono::high_resolution_clock::now();

    MPI_Init(&argc, &argv);

    int size, rank;
    int start, end;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> starts(size);
    std::vector<int> ends(size);
    std::vector<int> counts(size);
    std::vector<int> displs(size);

    double dx = X / (N + 1);
    double dx2 = dx * dx;
    double dt = T / nt;

    double error = 0.0;

    for (int i = 1; i <= N; i++)
        my_sin[i - 1] = std::sin(i * dx);

    for (int i = 0; i <= nt; i++)
        my_sin_time[i] = std::sin(i * dt);

    if (rank == 0){

        for (int i = 0; i < N * N * N; i++)
            solution[i] = 0.0;

        for (int i = 0; i < size; i++){
            starts[i] = (double)(i * N) / size;
            ends[i] = (double)((i + 1) * N) / size;
            counts[i] = (ends[i] - starts[i]) * N * N;
            displs[i] = starts[i] * N * N;
        }

        start = starts[0];
        end = ends[0];

        for (int i = 1; i < size; i++){
            MPI_Request request, request2;
            MPI_Isend(&starts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &request);
            MPI_Isend(&ends[i], 1, MPI_INT, i, 1, MPI_COMM_WORLD, &request);
        }
    } else {
        MPI_Request request, request2;
        MPI_Irecv(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
        MPI_Irecv(&end, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &request);
    }

   // MPI_Barrier(MPI_COMM_WORLD);

//--------------------------------TIME LOOP------------------------------------//

    for (int x = 0; x < nt; x++){

        //-------------------------------RHS-----------------------------------//

        //  Bcast solution - not optimal but it works
        MPI_Bcast(&solution, N * N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        RHS(counts, displs, start, end, dt, dx2, x);

        //---------------------------X DIRECTION-------------------------------//

        parallel_x_direction(counts, displs, start, end, dt, dx2, x);

        //---------------------------Y DIRECTION-------------------------------//

        parallel_y_direction(counts, displs, start, end, dt, dx2, x); // TO DO

        //---------------------------Z DIRECTION-------------------------------//

        parallel_z_direction(counts, displs, start, end, dt, dx2, x); // TO DO

        //----------------------------FINALIZE---------------------------------//

        if (rank == 0)
            error = finalize(x);

    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    if (rank == 0){
        std::cout << "Error = " << error << std::endl;
        std::cout << "Time = " << duration.count() / pow(10, 6) / (N*N*N*N) << std::endl;
    }

    MPI_Finalize();

    return 0;
}