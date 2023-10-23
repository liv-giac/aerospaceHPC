#include <mpi.h>
#include <omp.h>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include <array>
#include <string>

//  g++ -O3 Heat3DTD.cpp -o Heat3d -ftree-vectorize -march=native

//----------------------------------VARIABLES----------------------------------//
#define X  1.0
#define T  1.0
#define N  100
#define nt 100
//#define n_proc 2

std::array<std::array<std::array<double, N>, N>, N> delta_solution;
std::array<std::array<std::array<double, N>, N>, N> solution;
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

    auto start = std::chrono::high_resolution_clock::now();

    //rhs
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){

                //external force
                delta_solution[i][j][k] = dt * my_sin[i] * my_sin[j] * my_sin[k] * (3 * std::sin(dt * (x + 0.5)) + std::cos(dt * (x + 0.5)));

                //old solution
                delta_solution[i][j][k] -= 6 * dt / dx2 * solution[i][j][k];

                //k - 1
                if(k != 0)
                    delta_solution[i][j][k] += dt / dx2 * solution[i][j][k-1];

                //k + 1
                if(k == N - 1)
                    delta_solution[i][j][k] += dt / dx2 * my_sin[i] * my_sin[j] * std::sin(X) * my_sin_time[x];
                else
                    delta_solution[i][j][k] += dt / dx2 * solution[i][j][k+1];

                //j - 1
                if(j != 0)
                    delta_solution[i][j][k] += dt / dx2 * solution[i][j-1][k];

                //j + 1
                if(j == N - 1)
                    delta_solution[i][j][k] += dt / dx2 * my_sin[i] * std::sin(X) * my_sin[k] * my_sin_time[x];
                else
                    delta_solution[i][j][k] += dt / dx2 * solution[i][j+1][k];

                //i - 1
                if(i != 0)
                    delta_solution[i][j][k] += dt / dx2 * solution[i-1][j][k];

                //i + 1
                if(i == N - 1)
                    delta_solution[i][j][k] += dt / dx2 * std::sin(X) * my_sin[j] * my_sin[k] * my_sin_time[x];
                else
                    delta_solution[i][j][k] += dt / dx2 * solution[i+1][j][k];
            }
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

    time1 += elapsed.count();
}

//---------------------------------X DIRECTION---------------------------------//
void x_direction(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);  

    auto start = std::chrono::high_resolution_clock::now();

    //boundary conditions
    for(int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
                delta_solution[i][j][N - 1] -= exdiag * my_sin[i] * my_sin[j] * std::sin(X) * (my_sin_time[x + 1] - my_sin_time[x]);

    diag[0] = val_diag;
    for (int k = 1; k < N; k++){
          diag[k] = val_diag - exdiag / diag[k-1] * exdiag;
    }

    //solve Thomas algorithm
    int inizio = (double)N/mpi_size * mpi_rank;
    int fine = (double)N/mpi_size * (mpi_rank + 1);
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            //riceve
            if (mpi_rank != 0)
            {
                MPI_Recv(&delta_solution[i][j][inizio - 1], MPI_DOUBLE, mpi_rank - 1, 0, MPI_COMM_WORLD);
            }

            for (int k = inizio; k < fine; k++){
                if (k == 0) continue;
                // diag[k] -= exdiag / diag[k-1] * exdiag;
                delta_solution[i][j][k] -= exdiag / diag[k-1] * delta_solution[i][j][k - 1];
            }

            //manda
            if (mpi_rank != mpi_size - 1)
            {
                MPI_Send(&delta_solution[i][j][fine - 1], MPI_DOUBLE, mpi_rank + 1, 0, MPI_COMM_WORLD);
            }

        }
    }    


     for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            if (mpi_rank == mpi_size - 1)
                 delta_solution[i][j][N - 1] = delta_solution[i][j][N - 1] / diag[N - 1];

            if (mpi_rank != mpi_size - 1)
            {
                MPI_Recv(&delta_solution[i][j][fine], MPI_DOUBLE, mpi_rank + 1, 0, MPI_COMM_WORLD);
            }
     

            for (int k = /*N-2*/ fine - 1; k >= inzio; k--){
                if (k == N-1) continue;
                delta_solution[i][j][k] = (delta_solution[i][j][k] - exdiag * delta_solution[i][j][k + 1]) / diag[k];
            }   

            if (mpi_rank != 0)
            {
                MPI_Send(&delta_solution[i][j][inizio], MPI_DOUBLE, mpi_rank - 1, 0, MPI_COMM_WORLD);
            } 
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

    time2 += elapsed.count();
}

//---------------------------------Y DIRECTION---------------------------------//
void y_direction(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    auto start = std::chrono::high_resolution_clock::now();

    //boundary conditions
    for(int i = 0; i < N; i++)
        for (int k = 0; k < N; k++)
                delta_solution[i][N - 1][k] -= exdiag * my_sin[i] * std::sin(X) * my_sin[k] * (my_sin_time[x + 1] - my_sin_time[x]);

    //solve Thomas algorithm
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            int pos = i * N * N + j * N;

            for (int k = 0; k < N; k++)
                diag[k] = val_diag;

            for (int k = 1; k < N; k++){
                diag[k] -= exdiag / diag[k-1] * exdiag;
                delta_solution[i][k][j] -= exdiag / diag[k-1] * delta_solution[i][k - 1][j];
            }

            delta_solution[i][N - 1][j] = delta_solution[i][N - 1][j] / diag[N - 1];
            for (int k = N-2; k >= 0; k--)
                delta_solution[i][k][j] = (delta_solution[i][k][j] - exdiag * delta_solution[i][k + 1][j]) / diag[k];
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
            delta_solution[N - 1][j][k] -= exdiag * std::sin(X) * my_sin[j] * my_sin[k] * (my_sin_time[x + 1] - my_sin_time[x]);

    //solve Thomas algorithm
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            int pos = i * N * N + j * N;

            for (int k = 0; k < N; k++)
                diag[k] = val_diag;

            for (int k = 1; k < N; k++){
                diag[k] -= exdiag / diag[k-1] * exdiag;
                delta_solution[k][i][j] -= exdiag / diag[k-1] * delta_solution[k - 1][i][j];
            }

            delta_solution[N - 1][i][j] = delta_solution[N - 1][i][j] / diag[N - 1];
            for (int k = N-2; k >= 0; k--)
                delta_solution[k][i][j] = (delta_solution[k][i][j] - exdiag * delta_solution[k + 1][i][j]) / diag[k];
        }
    }

    auto finish = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);

    time4 += elapsed.count();
}

//----------------------------------FINALIZE-----------------------------------//
void finalize(int x){

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                solution[i][j][k] += delta_solution[i][j][k];

    //error
    if(x == nt - 1){
        double error = 0.0;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                for (int k = 0; k < N; k++)
                    error += (solution[i][j][k] - my_sin[i] * my_sin[j] * my_sin[k] * my_sin_time[x + 1]) * (solution[i][j][k] - my_sin[i] * my_sin[j] * my_sin[k] * my_sin_time[x + 1]);
        error = std::sqrt(error) / std::pow(N, 3);

        std::cout << "Error = " << error << std::endl;
    }
}

//------------------------------------MAIN-------------------------------------//
int main(int argc, char** argv){

    // i = asse z (piano)
    // j = asse y (riga)
    // k = asse x (colonna)

    MPI_Init(&argc, &argv);

    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);  

    double dx = X / (N + 1);
    double dx2 = dx * dx;
    double dt = T / nt;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                solution[i][j][k] = 0.0;

    for (int i = 1; i <= N; i++)
        my_sin[i - 1] = std::sin(i * dx);

    for (int i = 0; i <= nt; i++)
        my_sin_time[i] = std::sin(i * dt);

//--------------------------------TIME LOOP------------------------------------//

    for (int x = 0; x < nt; x++){

        //-------------------------------RHS-----------------------------------//

        RHS(dt, dx2, x);

        //---------------------------X DIRECTION-------------------------------//

        x_direction(dt, dx2, x);

        //---------------------------Y DIRECTION-------------------------------//

        y_direction(dt, dx2, x);

        //---------------------------Z DIRECTION-------------------------------//

        z_direction(dt, dx2, x);

        //----------------------------FINALIZE---------------------------------//

        finalize(x);

    }

    // std::cout << "Time1 = " << time1 / 1e9 << std::endl;
    // std::cout << "Time2 = " << time2 / 1e9 << std::endl;
    // std::cout << "Time3 = " << time3 / 1e9 << std::endl;
    // std::cout << "Time4 = " << time4 / 1e9 << std::endl;

    MPI_Finalize();

    return 0;
}