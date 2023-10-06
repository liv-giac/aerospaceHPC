#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>
#include <vector>

#define Time 1.0
#define X 1.0
#define N 640

int main(int argc, char** argv){

    double a, b, dt, dx, T0, T1, nt, error, prec, t_dep;
    double* T = new double[N];
    double* sol = new double[N];

    for(double n = 10; n <= N; n *= 2){

        std::cout << "N: " << n << std::endl;

        nt = 4.0 * n * n;
        std::cout << "nt: " << nt << std::endl;

        error = 0.0;

        dx = 1 / (n + 1);
        dt = 1 / nt;

        a = dt / (dx * dx);
        b = -2 * a;

        for(int i = 0; i < n; i++)
            T[i] = 0;

        for(double t = 0; t <= Time; t += dt){

            T0 = 0;
            T1 = sin(X) * sin(t);
            t_dep = sin(t) + cos(t);

            for(int i = 0; i < n; i++){
                sol[i] = T[i] + b * T[i] + dt * sin((i + 1) * dx) * t_dep;
                if(i == 0)
                    sol[i] += a * (T0 + T[i+1]);
                if(i == n - 1)
                    sol[i] += a * (T[i-1] + T1);
                if(i > 0 && i < n - 1)
                    sol[i] += a * (T[i-1] + T[i+1]);
            }

            for(int i = 0; i < n; i++)
                T[i] = sol[i];
        }

        for(int i = 0; i < n; i++)
            sol[i] = sin((i + 1) * dx) * sin(Time);

        for(int i = 0; i < n; i++)
            error += std::pow(sol[i] - T[i], 2);
        error = std::sqrt(error) / n;

        std::cout << "Error: " << error << std::endl;
        std::cout << std::endl;
    }

    delete[] T;
    delete[] sol;

    return 0;
}