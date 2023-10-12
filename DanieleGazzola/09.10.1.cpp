#include <cstdlib>
#include <iostream>
#include <math.h>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <eigen3/unsupported/Eigen/SparseExtra>

#define X 1.0

//global variables
int N;

//-----------------------------SUPPORT FUNCTIONS-------------------------------//

int posx(int n){
    return (n % N) + 1;
}

int posy(int n){
    return ((n % (N * N)) / N) + 1;
}

int posz(int n){
    return (n / (N * N)) + 1;
}

double solution(double x, double y, double z){
    return sin(x) * sin(y) * sin(z);
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

    Eigen::SparseMatrix<double> mat(size, size);
    Eigen::VectorXd rhs(size);
    Eigen::VectorXd sol(size);
    Eigen::VectorXd exact_sol(size);

//-------------------------CREATION MATRIX AND VECORS--------------------------//

    //exact_sol Time = 0
    for(int i = 0; i < size; i++)
        exact_sol[i] = solution(posx(i) * dx, posy(i) * dx, posz(i) * dx);

    //sol Time = 0
    for(int i = 0; i < size; i++)
        sol[i] = 0;

    //rhs Time = 0
    for(int i = 0; i < size; i++){

        //forcing term
        rhs[i] = 3 * exact_sol[i];

        //boundary conditions

        //X = 1
        if(posx(i) == N)
            rhs[i] -= exdiag * solution(X, posy(i) * dx, posz(i) * dx);

        //Y = 1
        if(posy(i) == N)
            rhs[i] -= exdiag * solution(posx(i)* dx, X, posz(i) * dx);

        //Z = 1
        if(posz(i) == N)
            rhs[i] -= exdiag * solution(posx(i) * dx, posy(i) * dx, X);
    }

    //mat initialization
    for (int i = 0; i < size; i++) {

        //diagonal
        mat.coeffRef(i, i) = diag;

        //diagonal - 1
        if(i > 0 && i % N != 0)
            mat.coeffRef(i, i - 1) = exdiag;

        //diagnoal + 1
        if(i < size - 1 && (i + 1) % N != 0)
            mat.coeffRef(i, i + 1) = exdiag;

        //diagonal - N
        if(i >= N && posy(i) != 1)
            mat.coeffRef(i, i - N) = exdiag;

        //diagonal + N
        if(i < size - N && posy(i) != N)
            mat.coeffRef(i, i + N) = exdiag;

        //diagonal - N * N
        if(i >= N * N)
            mat.coeffRef(i, i - N * N) = exdiag;

        //diagonal + N * N
        if(i < size - N * N)
            mat.coeffRef(i, i + N * N) = exdiag;

    }

    mat.makeCompressed();

//-----------------------------COMPUTE SOLUTION---------------------------------//

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
    solver.compute(mat);
    sol = solver.solve(rhs);

//-----------------------------ERROR AND OUTPUT---------------------------------//

    for(int i = 0; i < size; i++)
        error += (sol[i] - exact_sol[i]) * (sol[i] - exact_sol[i]);
    error = sqrt(error) / size;

    std::cout << "Error: " << error << std::endl;

    return 0;
}