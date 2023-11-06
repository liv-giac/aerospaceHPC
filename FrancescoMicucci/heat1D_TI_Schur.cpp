#include <cmath>
#include <chrono>
#include <iostream>

const int n = 1000000;      // Number of points
const int numDecomp = 5;    // Number of decomposition of the system

// Implementation of the Thomas Algorithm.
// Used to solve a linear system of the form: A * solution = rhs
// where A is a tridiagonal matrix.
void thomasAlgorithm(   
                    double* solution,
                    double* diag,
                    double* upperDiag,
                    double* lowerDiag,
                    double* rhs,
                    int dim){
    double w;
    
    for(int i = 1; i < dim; i++){
        w = lowerDiag[i - 1] / diag[i - 1];
        diag[i] -= w * upperDiag[i - 1];
        rhs[i] -= w * rhs[i - 1];
    }
    solution[dim - 1] = rhs[dim - 1] / diag[dim - 1];
    for(int i = dim - 2; i >= 0; i--){
        solution[i] = (rhs[i] - upperDiag[i] * solution[i + 1]) / diag[i];
    }
}

// Modified version of the Thomas algorithm used for solving simultaneously 
// 2 different linear systems having the same matrix A and rhs as follow:
//      |b c        |
//      |a b c      |
//  A = |  a b c    |
//      |    a b c  |
//      |      a b c|
// The rhs of the 2 systems are the following:
//          |a|             |0|
//          |0|             |0|
//  rhs1 =  |0|      rhs2 = |0|
//          |0|             |0|
//          |0|             |c|
void schurSubsystemsSolver(
                    double x[], 
                    double y[], 
                    double diag, 
                    double upperDiag, 
                    double lowerDiag, 
                    int dim){
    // Wikipedia notation
    double w;
    double* b = new double[dim];
    double* d = new double[dim];

    std::fill_n(b, dim, diag);
    d[0] = lowerDiag;

    for(int i = 1; i < dim; i++){
        w = lowerDiag / b[i - 1];
        b[i] = b[i] - w * upperDiag;
        d[i] = d[i] - w * d[i-1];
    }
    x[dim - 1] = d[dim - 1] / b[dim - 1];
    y[dim - 1] = upperDiag / b[dim - 1];
    for(int i = dim - 2; i >= 0; i--){
        x[i] = (d[i] - upperDiag * x[i + 1]) / b[i];
        y[i] = (- upperDiag * y[i + 1]) / b[i];
    }
}

// Extract the elements of the rhs associated with the matrix S
void extractRhsInterface(double rhs[], double rhsInterface[], int dim){
    int stride = (dim - numDecomp + 1) / numDecomp;
    
    for(int i = 1; i < numDecomp; i++){
        rhsInterface[i - 1] = rhs[i * (stride + 1) - 1];
    }
}

// Extract the elements of the rhs associated with the matrix Aj
void extractRhsAj(double rhs[], double rhsAi[], int dim, int j){
    int stride = (dim - numDecomp + 1) / numDecomp, lenghtRhsAi = stride;
    if(j == numDecomp - 1) lenghtRhsAi = dim - numDecomp + 1 - ((numDecomp - 1) * stride);
    
    for(int i = 0; i < lenghtRhsAi; i++){
        rhsAi[i] = rhs[i + j * (stride + 1)];
    }
}

// Compute the error of the solution w.r.t. the exact solution.
double computeError(double* solution, double* exactSolution, int dim){
    double error = 0;
    for(int i = 0; i < dim; i++){
        error += (solution[i] - exactSolution[i]) * (solution[i] - exactSolution[i]);
    }
    error = std::sqrt(error);
    error = error / dim;

    return error;
}

int main(int argc, char* argv[])
{
    int dimSubmatrix = (n - numDecomp + 1) / numDecomp;         // Dim of the first n-1 submatrices A0, A1, ..., An-2
    int dimLatestSubmatrix =                                    // Dim of the latest submatrix An-1
        n - numDecomp + 1 - ((numDecomp - 1) * dimSubmatrix);
    double dx               = 1.0 / (n + 1);                    // Step size
    double diagElem         = 2.0 / (dx * dx);                  // Value of the diagonal elements of the initial matrix 
    double upperElem        = -1.0 / (dx * dx);                 // Value of the elements of the 1st upper diagonal of the initial matrix
    double lowerElem        = -1.0 / (dx * dx);                 // Value of the elements of the 1st lower diagonal of the initial matrix
    double* rhs             = new double[n];                    // Rhs of the problem    
    double* exactSolution   = new double[n];                    // Exact solution of the problem
    double* solution        = new double[n];                    // Computed solution
    double* diagS           = new double[numDecomp - 1];        // Schur complement main diagonal
    double* upperDiagS      = new double[numDecomp - 2];        // Schur complement 1st upper diagonal
    double* lowerDiagS      = new double[numDecomp - 2];        // Schur complement 1st lower diagonal 
    double* xi              = new double[dimSubmatrix];         // xi corresponds to the column i-1 of (Ai^(-1)*Di)
    double* yi              = new double[dimSubmatrix];         // yi corresponds to the column i of (Ai^(-1)*Di)
    double* xm              = new double[dimLatestSubmatrix];   // xm corresponds to the column n-2 of (An-1^(-1)*Dn-1)
    double* ym              = new double[dimLatestSubmatrix];   // [in this case ym does not exist, we allocate memory just to use the function]
    double* rhsInterface    = new double[numDecomp - 1];        // Rhs elements associated with the matrix S
    double* xInt            = new double[numDecomp - 1];        // Solution of the linear system S * xInt = rhsInterface  

    // Vector storing the solutions of the linear systems: Ai * y = subrhs
    // The 1st element of the i-th solution is contained at the cell (i-1) * dimSubmatrix
    // The successive elements comes after the cell containing the 1st one.
    double* yStorage = new double[n - numDecomp + 1];

    auto start = std::chrono::high_resolution_clock::now();

    // Initialize rhs and exaxt solution of the pb.
    for (int i = 0; i <= n-1; i++){
        exactSolution[i] = std::sin((i + 1) * dx); 
        
        if(i < n - 1)
            rhs[i] = std::sin((i + 1) * dx);
        else
            rhs[n-1] = std::sin(1.0) / (dx * dx) + std::sin((i + 1) * dx); 
    }

    extractRhsInterface(rhs, rhsInterface, n);
    
    schurSubsystemsSolver(xi, yi, diagElem, upperElem, lowerElem, dimSubmatrix);

    // When we perform the following matrix multiplication Ei*Ai^(-1)*Di we get as result the following matrix:
    //                        i-1   i
    //                  |                   |
    // Ei*Ai^(-1)*Di =  |     s11  s12      |    i-1
    //                  |     s21  s22      |    i
    //                  |                   |  
    double s11 = upperElem * xi[0];
    double s12 = upperElem * yi[0];
    double s21 = lowerElem * xi[dimSubmatrix - 1];
    double s22 = lowerElem * yi[dimSubmatrix - 1];

    schurSubsystemsSolver(xm, ym, diagElem, upperElem, lowerElem, dimLatestSubmatrix);

    // When we perform the following matrix multiplication En-1*An-1^(-1)*Dn-1 we get as result the following matrix:
    //                                        n-2            
    //                        |                   |
    // En-1*An-1^(-1)*Dn-1 =  |                   |    
    //                        |                   |    
    //                        |               smm |    n-2
    double smm = upperElem * xm[0];

    // Assembling Schur complement matrix: S = F - E0*A0^(-1)*D0 - E1*A1^(-1)*D1 - ...
    // Obs: F = diag(diagElem)
    std::fill_n(diagS, numDecomp, diagElem);
    diagS[0] -= s22;
    diagS[numDecomp - 2] -= smm;
    for(int i = 1; i < numDecomp - 1; i++){
        diagS[i-1]      -= s11;
        diagS[i]        -= s22;
        upperDiagS[i-1] -= s12;
        lowerDiagS[i-1] -= s21;
    }

    double* subrhs = new double[dimSubmatrix];               // Rhs elements associated with the Matrix Ai
    double* y = new double[dimSubmatrix];                    // Solution of the linear system: Ai * y = subrhs 
    double* diag = new double[dimSubmatrix];                 // Main diagonal of the Matrix Ai
    double* upperDiag = new double[dimSubmatrix - 1];        // 1st upper diagonal of the Matrix Ai
    double* lowerDiag = new double[dimSubmatrix - 1];        // 1st lower diagonal of the Matrix Ai

    // Init diag, upperDiag, lowerDiag, subrhs of Matrix A0
    std::fill_n(diag, dimSubmatrix, diagElem);
    std::fill_n(upperDiag, dimSubmatrix - 1, upperElem);
    std::fill_n(lowerDiag, dimSubmatrix - 1, lowerElem);
    extractRhsAj(rhs, subrhs, n, 0);

    // Solve A0 * y = subrhs 
    thomasAlgorithm(y, diag, upperDiag, lowerDiag, subrhs, dimSubmatrix);

    // Store y
    for(int j = 0; j < dimSubmatrix; j++) yStorage[j] = y[j];
    
    // Update rhsInterface
    rhsInterface[0] -= lowerElem * y[dimSubmatrix - 1]; 

    for(int i = 1; i < numDecomp - 1; i++){
        // Init diag, upperDiag, lowerDiag, subrhs of Matrix Ai
        std::fill_n(diag, dimSubmatrix, diagElem);
        std::fill_n(upperDiag, dimSubmatrix - 1, upperElem);
        std::fill_n(lowerDiag, dimSubmatrix - 1, lowerElem); 
        extractRhsAj(rhs, subrhs, n, i);
        
        // Solve Ai * y = subrhs 
        thomasAlgorithm(y, diag, upperDiag, lowerDiag, subrhs, dimSubmatrix);

        // Store y
        for(int j = 0; j < dimSubmatrix; j++) yStorage[i * dimSubmatrix + j] = y[j];

        // Update rhsInterface
        rhsInterface[i - 1] -= upperElem * y[0];
        rhsInterface[i]     -= lowerElem * y[dimSubmatrix - 1]; 
    }

    subrhs = new double[dimLatestSubmatrix];
    y = new double[dimLatestSubmatrix];
    diag = new double[dimLatestSubmatrix];
    upperDiag = new double[dimLatestSubmatrix - 1];
    lowerDiag = new double[dimLatestSubmatrix - 1];

    // Init diag, upperDiag, lowerDiag, subrhs of Matrix An-1
    std::fill_n(diag, dimLatestSubmatrix, diagElem);
    std::fill_n(upperDiag, dimLatestSubmatrix - 1, upperElem);
    std::fill_n(lowerDiag, dimLatestSubmatrix - 1, lowerElem);
    extractRhsAj(rhs, subrhs, n, numDecomp - 1);

    // Solve An-1 * y = subrhs 
    thomasAlgorithm(y, diag, upperDiag, lowerDiag, subrhs, dimLatestSubmatrix);

    // Store y
    for(int j = 0; j < dimLatestSubmatrix; j++) yStorage[(numDecomp - 1) * dimSubmatrix + j] = y[j];

    // Update rhsInterface
    rhsInterface[numDecomp - 2] -= upperElem * y[0];     

    // Solve S * xInt = rhsInterface
    thomasAlgorithm(xInt, diagS, upperDiagS, lowerDiagS, rhsInterface, numDecomp - 1);

    // Compute the solution for all the points
    for(int j = 0; j < numDecomp - 1; j++) solution[(j + 1) * dimSubmatrix + j] = xInt[j];
    for(int j = 0; j < dimSubmatrix; j++){
        solution[j] = yStorage[j];
        solution[j] -= yi[j] * xInt[0];
    }
    for(int i = 1; i < numDecomp - 1; i++){
        for(int j = 0; j < dimSubmatrix; j++){
            solution[i * (dimSubmatrix + 1) + j] = yStorage[i * dimSubmatrix + j];
            solution[i * (dimSubmatrix + 1) + j] -= xi[j] * xInt[i - 1];
            solution[i * (dimSubmatrix + 1) + j] -= yi[j] * xInt[i];
        }     
    }
    for(int j = 0; j < dimLatestSubmatrix; j++){
        solution[(numDecomp - 1) * (dimSubmatrix + 1) + j] = yStorage[(numDecomp - 1) * dimSubmatrix + j];
        solution[(numDecomp - 1) * (dimSubmatrix + 1) + j] -= xm[j] * xInt[numDecomp - 2];
    }     

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = end - start;

    std::cout << "============================================" << std::endl;
    std::cout << "Error with " << n << " points: " << computeError(solution, exactSolution, n) << std::endl;
    std::cout << "Time: " << time.count() * 1000 << " milliseconds " << std::endl;
    std::cout << "Seconds per point: " << time.count() / n << std::endl;
    std::cout << "============================================" << std::endl;

    // Free memory
    delete[] rhs;
    delete[] exactSolution;
    delete[] diagS;
    delete[] upperDiagS;
    delete[] lowerDiagS;
    delete[] xi;
    delete[] yi;
    delete[] xm;
    delete[] ym;
    delete[] rhsInterface;
    delete[] xInt;
    delete[] yStorage;
    delete[] subrhs;
    delete[] y;
    delete[] diag;
    delete[] upperDiag;
    delete[] lowerDiag;
}