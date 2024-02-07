#ifndef SCHURSOLVER_HPP
#define SCHURSOLVER_HPP

#define NDEBUG

#define A_PLUS 0
#define A_MINUS 1

#define C_PLUS 0
#define C_MINUS 1

#define _OUT
#define _INOUT

#define _MPI_LOCAL
#define _MPI_GLOBAL

#include <cmath>
#include <iostream>
#include <mpi.h>

/*
 * Class implementing the Schur complement method for solving the heat equation.
 */
class SchurSolver
{
public:
    // Constructor.
    SchurSolver(const MPI_Comm &comm_, const unsigned int &n_, const int &numDecomp_,
                const double &dx_, double diagElem_,
                double upperElem_, double lowerElem_)
        : n(n_), local_n((n - numDecomp_ + 1) / numDecomp_), numDecomp(numDecomp_), schurSize(numDecomp_ - 1), dx(dx_), comm(comm_), diagElem(diagElem_), upperElem(upperElem_), lowerElem(lowerElem_)
    {
        // get MPI ID
        MPI_Comm_rank(comm, &rank);

        // get MPI size
        MPI_Comm_size(comm, &size);

        // print size and rank
        std::cout << "Rank: " << rank << " size: " << size << std::endl;

        error = 0.0;
        dimSubmatrix = (n - numDecomp + 1) / numDecomp;
        dimLatestSubmatrix = n - numDecomp + 1 - ((numDecomp - 1) * dimSubmatrix);

        setup();
    }

    // Destructor
    ~SchurSolver()
    {
        delete[] matrixDiag;
        delete[] matrixUpperDiag;
        delete[] matrixLowerDiag;
        delete[] localRhs;
        delete[] schurRhs;
    }

    // Solve the heat equation in 1D.
    void solve();

    // Return the error associated with the computed solution.
    double getError() const { return error; }

protected:
    const unsigned int n; // Number of points
    // TODO initialise
    const unsigned int local_n;            // Number of local points EXCLUDING the interface points
    const unsigned int numDecomp;          // Number of decomposition of the system
    const unsigned int schurSize; // Size of the Schur complement
    const double dx;              // Step size
    double *rhs;            // Rhs of the problem
    double *exactSolution;  // Exact solution of the problem
    _MPI_LOCAL double exactSolutionInterface[2]; // Solution at the interface points

    // MPI rank
    int rank;
    // MPI size
    int size;
    // MPI communicator
    MPI_Comm comm;

    _MPI_LOCAL double *matrixDiag;          // Diagonal elements of the initial matrix
    _MPI_LOCAL double *matrixUpperDiag;     // Upper diagonal elements of the initial matrix
    _MPI_LOCAL double *matrixLowerDiag;     // Lower diagonal elements of the initial matrix
    _MPI_LOCAL double *localRhs;            // Local rhs of the problem
    _MPI_LOCAL double lateralElements_D[2]; // First element is a+ and the second is c-
    _MPI_LOCAL double bottomElements_E[2];  // First element is c+ and the second is a-
    _MPI_GLOBAL double *schurRhs;           // Schur rhs; shared by all processors.
    _MPI_LOCAL double solutionInterfaceValues[2]; // First element is the value of the solution at the interface point i-1 and the second is the value of the solution at the interface point i

    double diagElem;                  // Value of the diagonal elements of the initial matrix
    double upperElem;                 // Value of the elements of the 1st upper diagonal of the initial matrix
    double lowerElem;                 // Value of the elements of the 1st lower diagonal of the initial matrix
    double error;                     // Error of the computed solution
    double *solution;                 // Solution of the problem
    _MPI_GLOBAL double *rhsInterface; // Rhs elements associated with the matrix S. These are GLOBAL
    _MPI_GLOBAL double *diagS;        // Schur complement main diagonal
    _MPI_GLOBAL double *upperDiagS;   // Schur complement 1st upper diagonal
    _MPI_GLOBAL double *lowerDiagS;   // Schur complement 1st lower diagonal
    unsigned int dimSubmatrix;                 // Dim of the first n-1 submatrices A0, A1, ..., An-2
    unsigned int dimLatestSubmatrix;           // Dim of the latest submatrix An-1

    double *xi; // xi corresponds to the column i-1 of (Ai^(-1)*Di)
    double *yi; // yi corresponds to the column i of (Ai^(-1)*Di)
    double *xm; // xm corresponds to the column n-2 of (An-1^(-1)*Dn-1)
    double *ym; // [in this case ym does not exist, we allocate memory just to use the function]

    // Vector storing the solutions of the linear systems: Ai * y = subrhs
    // The 1st element of the i-th solution is contained at the cell (i-1) * dimSubmatrix
    // The successive elements comes after the cell containing the 1st one.
    double *yStorage;

    // Initialize all the class components
    void setup();

    // Compute the elements componing the Schur complement
    void computeSchurComplement();

    /**
     * @brief Build the rhs of the Schur complement system. In the end, each processor will own a copy of schurRhs.
     *
     */
    void updateSchurRhs();

    // Compute the solution of the heat equation in 1D.
    void computeSolution();

    // Compute the error of the solution w.r.t. the exact solution.
    void computeError();

    // Build matrix
    void buildMatrix();

    // Build rhs
    void buildRhs();

    // Build lateral elements
    void buildLateralElements();

    // Build exact solution
    void buildExactSolution();

    // Implementation of the Thomas Algorithm.
    // Used to solve a linear system of the form: A * solution = rhs
    // where A is a tridiagonal matrix.
    /*[[deprecated]] void thomasAlgorithm(_OUT double *solution, _INOUT double *diag,
                                        double *upperDiag, double *lowerDiag,
                                        _INOUT double *rhs, int dim);*/

    void thomasAlgorithm(_OUT double *solution, const double *const diag,
                         const double *const upperDiag, const double *const lowerDiag,
                         const double *const rhs, const unsigned int &dim) const;

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
    void schurSubsystemsSolver(double x[], double y[],
                               double diag, double upperDiag,
                               double lowerDiag, int dim);

    // Extract the rhs elements associated with the matrix S from the rhs of the problem.
    void extractRhsInterface();

    // Extract the rhs elements associated with the matrix Ai from the rhs of the problem.
    void extractRhsAj(double rhsAi[], unsigned int j);
};

#endif