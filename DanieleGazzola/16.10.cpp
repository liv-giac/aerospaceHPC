#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>
#include <array>
#include <string>

//----------------------------------VARIABLES----------------------------------//
#define X  1.0
#define T  1.0
#define N  100
#define nt 100

std::array<std::array<std::array<double, N>, N>, N> delta_solution;
std::array<std::array<std::array<double, N>, N>, N> solution;
std::array<std::array<std::array<double, N>, N>, N> rhs;
std::array<double, N> diag;
std::array<double, N> my_sin;
std::array<double, nt + 1> my_sin_time;

double time_rhs = 0.0;
double time_delta = 0.0;
double time_thomas = 0.0;

//----------------------------------X DIRECTION---------------------------------//
void x_direction(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    double w;

    auto start = std::chrono::high_resolution_clock::now();

    //rhs
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){

                //external force
                rhs[i][j][k] = dt * my_sin[i] * my_sin[j] * my_sin[k] * (3 * std::sin(dt * (x + 1 / 2)) + std::cos(dt * (x + 1 / 2)));

                //old solution
                rhs[i][j][k] -= 6 * dt / dx2 * solution[i][j][k];

                //k - 1
                if(k != 0)
                    rhs[i][j][k] += dt / dx2 * solution[i][j][k-1];

                //k + 1
                if(k == N - 1)
                    rhs[i][j][k] += dt / dx2 * my_sin[i] * my_sin[j] * std::sin(X) * my_sin_time[x];
                else
                    rhs[i][j][k] += dt / dx2 * solution[i][j][k+1];

                //j - 1
                if(j != 0)
                    rhs[i][j][k] += dt / dx2 * solution[i][j-1][k];

                //j + 1
                if(j == N - 1)
                    rhs[i][j][k] += dt / dx2 * my_sin[i] * std::sin(X) * my_sin[k] * my_sin_time[x];
                else
                    rhs[i][j][k] += dt / dx2 * solution[i][j+1][k];

                //i - 1
                if(i != 0)
                    rhs[i][j][k] += dt / dx2 * solution[i-1][j][k];

                //i + 1
                if(i == N - 1)
                    rhs[i][j][k] += dt / dx2 * std::sin(X) * my_sin[j] * my_sin[k] * my_sin_time[x];
                else
                    rhs[i][j][k] += dt / dx2 * solution[i+1][j][k];
            }
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    auto start1 = std::chrono::high_resolution_clock::now();

    //boundary conditions
    for(int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                //k = N
                if(k == N - 1)
                    rhs[i][j][k] -= exdiag * my_sin[i] * my_sin[j] * std::sin(X) * (my_sin_time[x + 1] - my_sin_time[x]);

    auto stop1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);

    auto start2 = std::chrono::high_resolution_clock::now();

    //solve Thomas algorithm
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            int pos = i * N * N + j * N;

            for (int k = 0; k < N; k++)
                diag[k] = val_diag;

            for (int k = 1; k < N; k++){
                w = exdiag / diag[k-1];
                diag[k] -= w * exdiag;
                rhs[i][j][k] -= w * rhs[i][j][k - 1];
            }

            delta_solution[i][j][N - 1] = rhs[i][j][N - 1] / diag[N - 1];
            for (int k = N-2; k >= 0; k--)
                delta_solution[i][j][k] = (rhs[i][j][k] - exdiag * delta_solution[i][j][k + 1]) / diag[k];
        }
    }

    auto stop2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);

    time_rhs += duration.count();
    time_delta += duration1.count();
    time_thomas += duration2.count();

    if (x == nt - 1){
        std::cout << "RHS time = " << time_rhs / pow(10, 6) << " seconds" << std::endl;
        std::cout << "Delta time = " << time_delta / pow(10, 6) << " seconds" << std::endl;
        std::cout << "Thomas time = " << time_thomas / pow(10, 6) << " seconds" << std::endl;
    }
}

//----------------------------------Y DIRECTION---------------------------------//
void y_direction(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    double w;

    //boundary conditions
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){

                rhs[i][j][k] = delta_solution[i][j][k];

                //j = N
                if(j == N - 1)
                    rhs[i][j][k] -= exdiag * my_sin[i] * std::sin(X) * my_sin[k] * (my_sin_time[x + 1] - my_sin_time[x]);
            }
        }
    }

    //solve Thomas algorithm
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            int pos = i * N * N + j * N;

            for (int k = 0; k < N; k++)
                diag[k] = val_diag;

            for (int k = 1; k < N; k++){
                w = exdiag / diag[k-1];
                diag[k] -= w * exdiag;
                rhs[i][k][j] -= w * rhs[i][k - 1][j];
            }

            delta_solution[i][N - 1][j] = rhs[i][N - 1][j] / diag[N - 1];
            for (int k = N-2; k >= 0; k--)
                delta_solution[i][k][j] = (rhs[i][k][j] - exdiag * delta_solution[i][k + 1][j]) / diag[k];
        }
    }
}

//----------------------------------Z DIRECTION---------------------------------//
void z_direction(double dt, double dx2, int x){

    double val_diag = 1.0 + dt / dx2;
    double exdiag = -dt / 2.0 / dx2;

    double w;

    //boundary conditions
    for(int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){

                rhs[i][j][k] = delta_solution[i][j][k];

                //i = N
                if(i == N - 1)
                    rhs[i][j][k] -= exdiag * std::sin(X) * my_sin[j] * my_sin[k] * (my_sin_time[x + 1] - my_sin_time[x]);
            }
        }
    }

    //solve Thomas algorithm
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){

            int pos = i * N * N + j * N;

            for (int k = 0; k < N; k++)
                diag[k] = val_diag;

            for (int k = 1; k < N; k++){
                w = exdiag / diag[k-1];
                diag[k] -= w * exdiag;
                rhs[k][i][j] -= w * rhs[k - 1][i][j];
            }

            delta_solution[N - 1][i][j] = rhs[N - 1][i][j] / diag[N - 1];
            for (int k = N-2; k >= 0; k--)
                delta_solution[k][i][j] = (rhs[k][i][j] - exdiag * delta_solution[k + 1][i][j]) / diag[k];
        }
    }
}

//-----------------------------------FINALIZE-----------------------------------//
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

//----------------------------------MAIN---------------------------------------//
int main(int argc, char** argv){

    // i = asse z (piano)
    // j = asse y (riga)
    // k = asse x (colonna)

    auto start1 = std::chrono::high_resolution_clock::now();

    double dx = X / N;
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

    auto stop1 = std::chrono::high_resolution_clock::now();

    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);

    std::cout << "Initialization time = " << duration1.count() / pow(10, 6) << " seconds" << std::endl;

    double x_count = 0.0;
    double y_count = 0.0;
    double z_count = 0.0;
    double f_count = 0.0;

//--------------------------------TIME LOOP------------------------------------//

    for (int x = 0; x < nt; x++){

        //---------------------------X DIRECTION-------------------------------//

        auto start2 = std::chrono::high_resolution_clock::now();
        x_direction(dt, dx2, x);
        auto stop2 = std::chrono::high_resolution_clock::now();
        auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(stop2 - start2);
        x_count += duration2.count();

        //---------------------------Y DIRECTION-------------------------------//

        auto start3 = std::chrono::high_resolution_clock::now();
        y_direction(dt, dx2, x);
        auto stop3 = std::chrono::high_resolution_clock::now();
        auto duration3 = std::chrono::duration_cast<std::chrono::microseconds>(stop3 - start3);
        y_count += duration3.count();

        //---------------------------Z DIRECTION-------------------------------//

        auto start4 = std::chrono::high_resolution_clock::now();
        z_direction(dt, dx2, x);
        auto stop4 = std::chrono::high_resolution_clock::now();
        auto duration4 = std::chrono::duration_cast<std::chrono::microseconds>(stop4 - start4);
        z_count += duration4.count();

        //----------------------------FINALIZE---------------------------------//

        auto start5 = std::chrono::high_resolution_clock::now();
        finalize(x);
        auto stop5 = std::chrono::high_resolution_clock::now();
        auto duration5 = std::chrono::duration_cast<std::chrono::microseconds>(stop5 - start5);
        f_count += duration5.count();

    }

    std::cout << "X direction time = " << x_count / pow(10, 6) << " seconds" << std::endl;
    std::cout << "Y direction time = " << y_count / pow(10, 6) << " seconds" << std::endl;
    std::cout << "Z direction time = " << z_count / pow(10, 6) << " seconds" << std::endl;
    std::cout << "Finalize time = " << f_count / pow(10, 6) << " seconds" << std::endl;

    return 0;
}