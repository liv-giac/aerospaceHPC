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

    int N[5] = {11, 21, 41, 81, 161};
    for (int n:N) {
        const double deltaX = 1.0/(n-1);
        const double deltaT = 0.25*deltaX*deltaX;
        double t = 0.0;
        SpVec solution(n);
        SpVec new_solution(n);
        SpVec real_solution(n);
        double max_error = 0.0;
        double error = 0.0;

        auto start_time_solving
            = std::chrono::system_clock::now();

        //Initialization
        for(int i=0; i < n; i++){
            solution(i) = std::sin(i * deltaX)*std::sin(0); 
        }

        t += deltaT;
        while(t <= 1.0){
            for (int i = 1; i < n - 1; i++){
                new_solution(i) = solution(i) + (deltaT/(deltaX*deltaX))*(solution(i-1) - 2*solution(i) + solution(i+1)) + deltaT*std::sin(i * deltaX)*(std::cos(t) + std::sin(t));
            }
            new_solution(0) = std::sin(0)*std::sin(t);
            new_solution(n - 1) = std::sin(1) * std::sin(t);
            solution = new_solution;
            t += deltaT;

            for(int i=0; i < n; i++){
                real_solution(i) = std::sin(i * deltaX)*std::sin(t); 
            }

            //Calculate Error
            error = 0.0;
            for(int i=1; i < n - 1; i++)
                    error += (solution(i) - real_solution(i)) * (solution(i) - real_solution(i));
            error = std::sqrt(error);
            error = error / n;

            //Save onòy the max error for every temporal step
            if(error > max_error){
                max_error = error;
            }
            
        }

        auto end_time
            = std::chrono::system_clock::now();

        std::chrono::duration<double> solving_duration 
            = end_time - start_time_solving;

        std::cout << "Nx: " << n << " Nt: " << 1/deltaT << std::endl;
        std::cout << "Error with " << n << " points: " << max_error << " delta x: " << deltaX << std::endl; 
        std::cout << "Total time: "
                    << (solving_duration.count())/n << " seconds/points " << std::endl << std::endl;

        fout << deltaX << "," << max_error << "\n";
    }

    fout.close();

}