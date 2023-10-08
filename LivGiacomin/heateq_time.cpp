#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>


int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <integer>" << std::endl;
        return 1; // Return an error code
    }

    int n = std::atoi(argv[1]);
    double Time = 1.0;
    double X = 1.0;
    double dx=X/(n+1);
    double b,a,error,T0,T1,prec,t_dep;
    double* T = new double[n];
    double* sol = new double[n];
    std::chrono::duration<double> elapsed;
    double nt = 4.0 * n*n;
    double dt=1.0/nt;
    a = dt/ (dx*dx);
    b = -2.0 * a;
    error = 0.0;
 
    for(int i = 0; i < n; i++)
            T[i] = 0;
        
    auto time1= std::chrono::system_clock::now();
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

    auto time2 = std::chrono::system_clock::now();
    for(int i = 0; i < n; i++)
            error += std::pow(sol[i] - T[i], 2);
        error = std::sqrt(error)/n;
    // elapsed = time2 - time1;

        std::cout << "N: " << n << std::endl;
        std::cout << "Error: " << error << std::endl;
        std::cout << "Error scaled by n: " << error/n << std::endl;
         std::cout << "Iterations time: " <<elapsed.count()<<std::endl;
         std::cout << "Iterations time scaled by n: " << elapsed.count() / n << std::endl;

    delete[] T;
    delete[] sol;

    return 0;
}
    
    
    