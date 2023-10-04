#include <chrono>
#include <cstdlib>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

    using SpMat = Eigen::SparseMatrix<double>;
    using SpVec = Eigen::VectorXd;


int main(int argc, char* argv[])
{
    std::ofstream fout("data.csv");
    fout << "deltaX,"  << "err" << "\n";

    int N[4] = {10, 20, 40, 80};
    for (int n:N) {
        const double deltaX = 1.0/(n-1);
        const double deltaT = 0.25*deltaX*deltaX;


        SpVec f(n);
        SpVec solution(n);
        SpVec b(n);
        double t = 0.0;

        auto start_time_solving
            = std::chrono::system_clock::now();

        // a = -1/(deltaX*deltaX);


        // for (int i = 0; i <= n - 1; i++){
        //     b(i) = 2/(deltaX*deltaX);
        //     f(i) = deltaT*std::sin((i + 1) * deltaX)*std::cos(t);
        // }
        // // f(n-1) += (std::sin(1.0)*std::cos(t))/(deltaX*deltaX);

        // double w;
        // for(int i = 1; i < n; i++){
        //     w = a/b(i - 1);
        //     b(i) -= w*a;
        //     f(i) -= w*f(i-1);
        // }

        // solution(n - 1) = f(n - 1)/b(n-1); 
        // for(int i = n - 2; i >= 0; i-- ){
        //     solution(i) = (f(i) - a*solution(i + 1))/b(i); 
        // }
        for(int i=0; i < n; i++){
            solution(i) = std::sin(i * deltaX)*std::sin(0); 
        }

        t += deltaT;
        SpVec new_solution(n);
        while(t <= 1.0){
            for (int i = 1; i < n - 1; i++){
                new_solution(i) = solution(i) + (deltaT/(deltaX*deltaX))*(solution(i-1) - 2*solution(i) + solution(i+1)) + deltaT*std::sin(i * deltaX)*(std::cos(t) + std::sin(t));
            }
            new_solution(0) = std::sin(0)*std::sin(t);
            new_solution(n - 1) = std::sin(1) * std::sin(t);
            solution = new_solution;
            t += deltaT;
        }

        auto end_time
            = std::chrono::system_clock::now();

        std::chrono::duration<double> solving_duration 
            = end_time - start_time_solving;

        SpVec real_solution(n);
        for(int i=0; i < n; i++){
            real_solution(i) = std::sin(i * deltaX)*std::sin(1); 
        }

        double error = 0;
        for(int i=1; i < n - 1; i++)
                error += (solution(i) - real_solution(i)) * (solution(i) - real_solution(i));
        error = std::sqrt(error);
        error = error / n;

        std::cout << "Nx: " << n << " Nt: " << 1/deltaT << std::endl;
        std::cout << "Error with " << n << " points: " << error << " delta x: " << deltaX << std::endl; 


        std::cout << "Total time: "
                    << (solving_duration.count())/n << " seconds/points " << std::endl << std::endl;

        fout << deltaX << "," << error << "\n";
    }

    fout.close();

}