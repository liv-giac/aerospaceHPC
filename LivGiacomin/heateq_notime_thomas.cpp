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
    double X = 1.0;
    double dx=X/(n+1);
    double w,a,error;
    double* T = new double[n];
    double* diag = new double[n];
    double* rhs = new double[n];
    double* sol = new double[n];
    std::chrono::duration<double> elapsed, elapsed2;
    error = 0;
    auto time1= std::chrono::system_clock::now();
    a = -1.0 / (dx * dx);
    for(int i = 0; i < n; i++){
            sol[i] = sin((i + 1) * dx);
            rhs[i] = sol[i];
            diag[i] = 2 / (dx * dx);
    }

    rhs[0] += 0;
    rhs[n-1] += sin(1) / (dx * dx);
    auto time2= std::chrono::system_clock::now();
    for(int i = 1; i < n; i++){
            w = a / diag[i-1];
            diag[i] -= w * a;
            rhs[i] -= w * rhs[i-1];
        }
    
    T[n-1] = rhs[n-1] / diag[n-1];
    for(int i = n-2; i >= 0; i--)
            T[i] = (rhs[i] - a * T[i+1]) / diag[i];
    auto time3 = std::chrono::system_clock::now();
    for(int i = 0; i < n; i++)
            error += (T[i] - sol[i])*(T[i] - sol[i]);
        error = std::sqrt(error);
    elapsed = time2 - time1;
    elapsed2 = time3- time2;
    auto tot=elapsed+elapsed2;
        std::cout << "N: " << n << std::endl;
        std::cout << "Error: " << error << std::endl;
        std::cout << "Error scaled by n: " << error << std::endl;
        std::cout << "Total time: " <<tot.count()<<" (of which "<<100*elapsed.count()/tot.count()<<"% setup and "<<100*elapsed2.count()/tot.count()<<"% solve) "<<std::endl;
        std::cout << "Setup time scaled by n: " << elapsed.count() / n << std::endl;
        std::cout << "Solve time scaled by n: " << elapsed2.count() / n << std::endl;
        
    delete[] T;
    delete[] diag;
    delete[] rhs;
    delete[] sol;

    return 0;
}
    
    
    