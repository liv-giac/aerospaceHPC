#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include <array>
#include <string>
#include <mpi.h>

//----------------------------------VARIABLES----------------------------------//
#define X  1.0
#define T  1.0
#define N  5
#define nt 1

std::array<double, N * N * N> delta_solution;
std::array<double, N * N * N> solution;
std::array<double, N> diag;
std::array<double, N> my_sin;
std::array<double, nt + 1> my_sin_time;

double time1 = 0.0;
double time2 = 0.0;
double time3 = 0.0;
double time4 = 0.0;

//-------------------------------------RHS-------------------------------------//
void RHS(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;
    double coeff = dt / dx2;

    double force = 3 * std::sin(dt * (x + 0.5)) + std::cos(dt * (x + 0.5));

    double sinx = std::sin(X);

    auto start = std::chrono::high_resolution_clock::now();

    //RHS

    for (int i = 0; i < N; i++){
        int pos_i = i * N * N;
        for (int j = 0; j < N; j++){
            int pos_j = j * N;
            for (int k = 0; k < N; k++){
                delta_solution[pos_i + pos_j + k] = dt * my_sin[i] * my_sin[j] * my_sin[k] * force - 6 * coeff * solution[pos_i + pos_j + k];

                //k - 1
                if (k != 0)
                    delta_solution[pos_i + pos_j + k] += coeff * solution[pos_i + pos_j + k - 1];

                //k + 1
                if (k == N - 1)
                    delta_solution[pos_i + pos_j + k] += coeff * my_sin[i] * my_sin[j] * sinx * my_sin_time[x];
                else
                    delta_solution[pos_i + pos_j + k] += coeff * solution[pos_i + pos_j + k + 1];

                //j - 1
                if (j != 0)
                    delta_solution[pos_i + pos_j + k] += coeff * solution[i * N * N + (j - 1) * N + k];

                //j + 1
                if (j == N - 1)
                    delta_solution[pos_i + pos_j + k] += coeff * my_sin[i] * sinx * my_sin[k] * my_sin_time[x];
                else
                    delta_solution[pos_i + pos_j + k] += coeff * solution[i * N * N + (j + 1) * N + k];

                //i - 1
                if (i != 0)
                    delta_solution[pos_i + pos_j + k] += coeff * solution[(i - 1) * N * N + j * N + k];

                //i + 1
                if (i == N - 1)
                    delta_solution[pos_i + pos_j + k] += coeff * sinx * my_sin[j] * my_sin[k] * my_sin_time[x];
                else
                    delta_solution[pos_i + pos_j + k] += coeff * solution[(i + 1) * N * N + j * N + k];
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

    time1 += elapsed.count();
}

//-----------------------------PARALLEL X DIRECTION----------------------------//
void parallel_x_direction(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    int size, rank;

    auto starttime = std::chrono::high_resolution_clock::now();

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

    std::vector<double> global_sol(N * N * N);
    std::vector<int> starts(size);
    std::vector<int> ends(size);
    std::vector<int> counts(size);
    std::vector<int> displs(size);

    int start, end;

    if (rank == 0){

        for (int k = 0; k < N; k++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    global_sol[k * N * N + i * N + j] = delta_solution[i * N * N + j * N + k];

        for (int i = 0; i < size; i++){
            starts[i] = (double)(i * N) / size;
            ends[i] = (double)((i + 1) * N) / size;
            counts[i] = (ends[i] - starts[i]) * N * N;
            displs[i] = starts[i] * N * N;
        }

        start = starts[0];
        end = ends[0];

        for (int i = 1; i < size; i++){
            MPI_Send(&starts[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&ends[i], 1, MPI_INT, i, 1, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&end, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    std::vector<double> local_sol((end - start) * N * N);

    MPI_Scatterv(global_sol.data(), counts.data(), displs.data(), MPI_DOUBLE, local_sol.data(), ((end - start) * N * N), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int batch;

    double prev;

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            batch = i * N + j;

            if (rank != 0)
                MPI_Recv(&prev, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int k = start; k < end; k++){
                if (k != 0)
                    local_sol[(k - start) * N * N + batch] -= exdiag / diag[k - 1] * prev;
                prev = local_sol[(k - start) * N * N + batch];
            }

            if (rank != size - 1)
                MPI_Send(&prev, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
    }

    double next;

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            batch = i * N + j;

            if (rank != size - 1)
                MPI_Recv(&next, 1, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            else
                local_sol[(end - start - 1) * N * N + batch] = local_sol[(end - start - 1) * N * N + batch] / diag[N - 1];

            for (int k = end - 1; k >= start; k--){
                if (k != N - 1)
                    local_sol[(k - start) * N * N + batch] = (local_sol[(k - start) * N * N + batch] - exdiag * next) / diag[k];
                next = local_sol[(k - start) * N * N + batch];
            }

            if (rank != 0)
                MPI_Send(&next, 1, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
        }
    }

    MPI_Gatherv(local_sol.data(), ((end - start) * N * N), MPI_DOUBLE, global_sol.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
        for (int k = 0; k < N; k++)
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    delta_solution[i * N * N + j * N + k] = global_sol[k * N * N + i * N + j];

    auto finishtime = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finishtime - starttime);

    time2 += elapsed.count();
}

//---------------------------------X DIRECTION---------------------------------//
void x_direction(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    auto start = std::chrono::high_resolution_clock::now();

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

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

    time3 += elapsed.count();
}

//---------------------------------Y DIRECTION---------------------------------//
void y_direction(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    auto start = std::chrono::high_resolution_clock::now();

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

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

    time3 += elapsed.count();
}

//---------------------------------Z DIRECTION---------------------------------//
void z_direction(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    auto start = std::chrono::high_resolution_clock::now();

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

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

    time4 += elapsed.count();
}

//----------------------------------FINALIZE-----------------------------------//
double finalize(int x){

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                solution[i * N * N + j * N + k] += delta_solution[i * N * N + j * N + k];

    //error
    if (x == nt - 1){
        double error = 0.0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++)
                    error += (solution[i * N * N + j * N + k] - my_sin[i] * my_sin[j] * my_sin[k] * my_sin_time[x + 1]) * (solution[i * N * N + j * N + k] - my_sin[i] * my_sin[j] * my_sin[k] * my_sin_time[x + 1]);
        error = std::sqrt(error) / std::pow(N, 3);
        return error;
    }

    return 0.0;
}

//------------------------------------MAIN-------------------------------------//
int main(int argc, char** argv){

    // i = asse z (piano)
    // j = asse y (riga)
    // k = asse x (colonna)

    MPI_Init(&argc, &argv);

    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double dx = X / (N + 1);
    double dx2 = dx * dx;
    double dt = T / nt;

    double error = 0.0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                solution[i * N * N + j * N + k] = 0.0;

    for (int i = 1; i <= N; i++)
        my_sin[i - 1] = std::sin(i * dx);

    for (int i = 0; i <= nt; i++)
        my_sin_time[i] = std::sin(i * dt);

//--------------------------------TIME LOOP------------------------------------//

    for (int x = 0; x < nt; x++){

        //-------------------------------RHS-----------------------------------//
        if (rank == 0)
        RHS(dt, dx2, x);

        //---------------------------X DIRECTION-------------------------------//

        if (size == 1)
            x_direction(dt, dx2, x);
        else
            parallel_x_direction(dt, dx2, x);

        //---------------------------Y DIRECTION-------------------------------//
        if (rank == 0)
        y_direction(dt, dx2, x);

        //---------------------------Z DIRECTION-------------------------------//
        if (rank == 0)
        z_direction(dt, dx2, x);

        //----------------------------FINALIZE---------------------------------//
        if (rank == 0)
        error = finalize(x);

    }

    if (rank == 0){
        std::cout << "Error = " << error << std::endl;
        std::cout << "Time1 = " << time1 / 1e9 << std::endl;
        std::cout << "Time2 = " << time2 / 1e9 << std::endl;
        std::cout << "Time3 = " << time3 / 1e9 << std::endl;
        std::cout << "Time4 = " << time4 / 1e9 << std::endl;
    }

    MPI_Finalize();

    return 0;
}