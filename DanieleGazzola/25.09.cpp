#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>

#define N 100000000
#define X 1.0

int main(int argc, char** argv){

    double dx, error, a, w;
    double* T = new double[N];
    double* diag = new double[N];
    double* rhs = new double[N];
    double* sol = new double[N];
    std::chrono::duration<double> elapsed, elapsed2;

    for(int n = 10; n <= N; ){

    //SETUP
        auto time1 = std::chrono::system_clock::now();

        dx = X / (n + 1);
        error = 0;

        //diagonali non principali
        a = -1 / (dx * dx);

        //soluzione esatta
        for(int i = 0; i < n; i++)
            sol[i] = sin((i + 1) * dx);

        //right hand side
        for(int i = 0; i < n; i++)
            rhs[i] = sol[i];

        //diagonale principale
        for(int i = 0; i < n; i++)
            diag[i] = 2 / (dx * dx);

        //condizioni al bordo
        rhs[0] += 0;
        rhs[n-1] += sin(1) / (dx * dx);

    //SOLVE THOMAS ALGORITHM (tridiagonal matrix)
        auto time2 = std::chrono::system_clock::now();

        for(int i = 1; i < n; i++){
            w = a / diag[i-1];
            diag[i] -= w * a;
            rhs[i] -= w * rhs[i-1];
        }

        T[n-1] = rhs[n-1] / diag[n-1];
        for(int i = n-2; i >= 0; i--)
            T[i] = (rhs[i] - a * T[i+1]) / diag[i];

        auto time3 = std::chrono::system_clock::now();

        //errore
        for(int i = 0; i < n; i++)
            error += (T[i] - sol[i])*(T[i] - sol[i]);
        error = std::sqrt(error);

        elapsed2 = time3 - time2;
        elapsed = time2 - time1;

    //STAMPA
        std::cout << "N: " << n << std::endl;
        std::cout << "Error: " << error << std::endl;
        std::cout << "Setup time: " << elapsed.count() / n << std::endl;
        std::cout << "Solve time: " << elapsed2.count() / n << std::endl;
        std::cout << "Total time: " << elapsed.count() / n + elapsed2.count() / n << std::endl;
        std::cout << std::endl;

        n *= 10;
    }

    delete[] T;
    delete[] diag;
    delete[] rhs;
    delete[] sol;

    return 0;
}