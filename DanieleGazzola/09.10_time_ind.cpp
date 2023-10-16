#include <cstdlib>
#include <iostream>
#include <math.h>
#include <chrono>
#include <fstream>
#include <unistd.h>
#include <string>
#include <sstream>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/unsupported/Eigen/SparseExtra>

#define X 1.0

//global variables
int N;

//-----------------------------SUPPORT FUNCTIONS-------------------------------//

void ram_usage(){
    pid_t pid = getpid();
    std::string statusFile = "/proc/" + std::to_string(pid) + "/status";
    std::ifstream file(statusFile.c_str());
    if (!file.is_open()) {
        std::cerr << "Impossibile aprire il file status del processo." << std::endl;
        return ;
    }
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("VmRSS") != std::string::npos) {
            break;
        }
    }
    file.close();

    //print ram usage
    std::istringstream iss(line);
    std::string trash;
    int ram;
    iss >> trash >> ram;

    std::cout << "RAM: " << ram << std::endl;
    /*
    std::ofstream outfile;
    outfile.open("RAM.csv", std::ios_base::app);
    outfile << N << "," << ram << std::endl;
    outfile.close();
    */
}

//-----------------------------------MAIN--------------------------------------//

int main(int argc, char** argv){

    //check input
    if(argc != 2){
        std::cout << "Usage: " << argv[0] << " N" << std::endl;
        return 1;
    }

    N = atoi(argv[1]);

    //variables declaration
    double dx = X / (N + 1);
    double dx2 = dx * dx;
    double error = 0.0;

    double diag = 6 / dx2;
    double exdiag = -1 / dx2;

    int size = N * N * N;

    double* sin = new double[N];

    Eigen::SparseMatrix<double> mat(size, size);
    Eigen::VectorXd rhs(size);
    Eigen::VectorXd sol(size);
    Eigen::VectorXd exact_sol(size);

//-------------------------CREATION MATRIX AND VECTORS-------------------------//

    //sin calculation
    for(int i = 1; i <= N; i++)
        sin[i - 1] = std::sin(i * dx);

    //exact_sol
    for(int i = 0; i < size; i++)
        exact_sol[i] = sin[i % N] * sin[(i % (N * N)) / N] * sin[i / (N * N)];

    //rhs
    for(int i = 0; i < size; i++){

        //forcing term
        rhs[i] = 3 * exact_sol[i];

        //boundary conditions

        //X = 1
        if(i % N == N - 1)
            rhs[i] -= exdiag * std::sin(X) * sin[(i % (N * N)) / N] * sin[i / (N * N)];

        //Y = 1
        if((i % (N * N)) / N == N - 1)
            rhs[i] -= exdiag * sin[i % N] * std::sin(X) * sin[i / (N * N)];

        //Z = 1
        if(i / (N * N) == N - 1)
            rhs[i] -= exdiag * sin[i % N] * sin[(i % (N * N)) / N] * std::sin(X);
    }

    auto start1 = std::chrono::high_resolution_clock::now();

    //mat initialization
    for (int i = 0; i < size; i++) {

        //diagonal - N * N
        if(i >= N * N)
            mat.insert(i, i - N * N) = exdiag;

        //diagonal - N
        if(i >= N)
            if ((i % (N * N)) / N != 0)
                mat.insert(i, i - N) = exdiag;

        //diagonal - 1
        if(i >= 1)
            if (i % N != 0)
                mat.insert(i, i - 1) = exdiag;

        //diagonal
        mat.insert(i, i) = diag;

        //diagnoal + 1
        if(i < size - 1)
            if (i % N != N - 1)
                mat.insert(i, i + 1) = exdiag;

        //diagonal + N
        if(i < size - N)
            if ((i % (N * N)) / N != N - 1)
                mat.insert(i, i + N) = exdiag;

        //diagonal + N * N
        if(i < size - N * N)
            mat.insert(i, i + N * N) = exdiag;
    }

    auto end1 = std::chrono::high_resolution_clock::now();

    mat.makeCompressed();

//-----------------------------COMPUTE SOLUTION---------------------------------//

    auto start2 = std::chrono::high_resolution_clock::now();

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(mat);
    sol = solver.solve(rhs);

    auto end2 = std::chrono::high_resolution_clock::now();

//----------------------------ERROR, TIME, OUTPUT-------------------------------//

    for(int i = 0; i < size; i++)
        error += (sol[i] - exact_sol[i]) * (sol[i] - exact_sol[i]);
    error = sqrt(error);

    auto elapsed1 = std::chrono::duration_cast<std::chrono::nanoseconds>(end1 - start1);
    auto elapsed2 = std::chrono::duration_cast<std::chrono::nanoseconds>(end2 - start2);

    std::cout << "Error: " << error / size << std::endl;
    std::cout << "Time for matrix initialization: " << elapsed1.count() / pow(10, 9) / size << " s / N^3" << std::endl;
    std::cout << "Time for solution: " << elapsed2.count() / pow(10, 9) / size << " s / N^3" << std::endl;

    // Controllo RAM
    ram_usage();

    delete[] sin;

    return 0;
}