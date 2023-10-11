#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using SpMat = Eigen::SparseMatrix<double>;
using SpVec = Eigen::VectorXd;

int main(int argc, char* argv[]){
    const int N[5] = {100, 1000, 10000, 100000, 1000000};
    std::ofstream fout("data.csv");
    fout << "n," << "deltaX," << "SecondsPerPoint," << "err" << "\n";

    for(int n:N){
        double deltaX = 1.0/(n+1), w, error = 0;
        SpMat A(n,n);
        SpVec f(n), solution(n), exact_solution(n);

        auto start_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i <= n-1; i++){
            A.coeffRef(i,i) = 2/(deltaX*deltaX);
            if(i < n - 1){
                A.coeffRef(i+1,i) = -1.0/(deltaX*deltaX);
                A.coeffRef(i,i+1) = -1.0/(deltaX*deltaX); 
                f(i) = std::sin((i + 1) * deltaX);
            }
            else{
                f(n-1) = std::sin(1.0)/(deltaX*deltaX) + std::sin((i + 1) * deltaX);
            }           
            exact_solution(i) = std::sin((i + 1) * deltaX); 
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> assembly_duration = end_time - start_time;

        start_time = std::chrono::high_resolution_clock::now();

        for(int i=1; i<n; i++){
            w = A.coeffRef(i, i-1) / A.coeffRef(i-1, i-1);
            A.coeffRef(i, i) = A.coeffRef(i, i) - w * A.coeffRef(i-1, i);
            f(i) = f(i) - w * f(i-1);
        }
        solution(n-1) = f(n-1)/(A.coeffRef(n-1, n-1));
        for(int i=n-2; i>=0; i--){
            solution(i) = (f(i) - A.coeffRef(i, i+1) * solution(i+1))/A.coeffRef(i, i);
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
    }

    fout.close();
}