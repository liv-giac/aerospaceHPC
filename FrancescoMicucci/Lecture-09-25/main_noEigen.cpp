#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>

int main(int argc, char* argv[])
{
    const int N[6] = {100, 200, 300, 400, 500, 600};
    std::ofstream fout("data.csv");
    fout << "n," << "deltaX," << "SecondsPerPoint," << "err" << "\n";

    for(int n:N){
        double deltaX = 1.0/(n+1), w, error = 0;
        double* diag = new double[n];
        double* f = new double[n];
        double* exact_solution = new double[n];
        double* solution = new double[n];
        double upper_under_diag = -1.0/(deltaX*deltaX);
        std::fill_n(diag, n, 2.0/(deltaX*deltaX));

        auto start_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i <= n-1; i++){
            if(i < n - 1)
                f[i] = std::sin((i + 1) * deltaX);
            else
                f[n-1] = std::sin(1.0)/(deltaX*deltaX) + std::sin((i + 1) * deltaX);          
            
            exact_solution[i] = std::sin((i + 1) * deltaX); 
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> assembly_duration = end_time - start_time;

        start_time = std::chrono::high_resolution_clock::now();
        
        for(int i=1; i<n; i++){
            w = upper_under_diag/diag[i-1];
            diag[i] = diag[i] - w * upper_under_diag;
            f[i] = f[i] - w * f[i-1];
        }
        solution[n-1] = f[n-1]/diag[n-1];
        for(int i=n-2; i>=0; i--){
            solution[i] = (f[i] - upper_under_diag * solution[i+1])/diag[i];
        }   

        end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> solve_duration = end_time - start_time;

        for(int i=0; i<n; i++)
            error += (solution[i] - exact_solution[i]) * (solution[i] - exact_solution[i]);
        error = std::sqrt(error);
        error = error / n;

        std::cout << "Error with " << n << " points: " << error << std::endl; 
        std::cout << "Assembly duration: " << assembly_duration.count()*1000 << " milliseconds " << std::endl;
        std::cout << "Solving duration: " << solve_duration.count()*1000 << " milliseconds " << std::endl;
        std::cout << "Seconds per point: " << (assembly_duration.count() + solve_duration.count())/n << std::endl << std::endl;

        fout << n << "," << deltaX << "," << (assembly_duration.count() + solve_duration.count())/n << "," << error << "\n";

        delete[] diag;
        delete[] f;
        delete[] solution;
        delete[] exact_solution;
    }

    fout.close();
}