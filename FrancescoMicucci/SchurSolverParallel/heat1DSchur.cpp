#include "heat1DSchur.hpp"

#include <cstring>

void SchurSolver::solve()
{
    updateSchurRhs();
    computeSolution();
    computeError();
}

void SchurSolver::setup()
{
    rhsInterface = new double[schurSize]; // Allocate memory for rhsInterface
    solution = new double[n];                 // Allocate memory for solution
    yStorage = new double[n - numDecomp + 1]; // Allocate memory for yStorage
    diagS = new double[schurSize];        // Allocate memory for diagS
    upperDiagS = new double[schurSize - 1];   // Allocate memory for upperDiagS
    lowerDiagS = new double[schurSize - 1];   // Allocate memory for lowerDiagS

    xi = new double[dimSubmatrix];       // Allocate memory for xi
    yi = new double[dimSubmatrix];       // Allocate memory for yi
    xm = new double[dimLatestSubmatrix]; // Allocate memory for xm
    ym = new double[dimLatestSubmatrix]; // Allocate memory for ym

    matrixDiag = new double[local_n];      // Allocate memory for matrixDiag
    matrixUpperDiag = new double[local_n]; // Allocate memory for matrixUpperDiag
    matrixLowerDiag = new double[local_n]; // Allocate memory for matrixLowerDiag

    localRhs = new double[local_n];   // Allocate memory for localRhs
    schurRhs = new double[numDecomp]; // There are (numProcessors - 1) elements of the Schur rhs

    extractRhsInterface(); // Extract rhsInterface from rhs

    buildMatrix();
    buildRhs();
    extractRhsInterface();
    buildLateralElements();
    computeSchurComplement();
}

void SchurSolver::computeSchurComplement()
{
    schurSubsystemsSolver(xi, yi, diagElem, upperElem, lowerElem, local_n);

    // When we perform the following matrix multiplication Ei*Ai^(-1)*Di we get as result the following matrix:
    //                        i-1   i
    //                  |                   |
    // Ei*Ai^(-1)*Di =  |     s11  s12      |    i-1
    //                  |     s21  s22      |    i
    //                  |                   |
    // Assemble local buffer
    double s_elems[4];
    s_elems[0] = bottomElements_E[C_PLUS] * xi[0];            // s11
    s_elems[1] = bottomElements_E[C_PLUS] * yi[0];            // s12
    s_elems[2] = bottomElements_E[A_MINUS] * xi[local_n - 1]; // s21
    s_elems[3] = bottomElements_E[A_MINUS] * yi[local_n - 1]; // s22

    // Every processor will assemble the Schur complement locally; we need every process to gather all the elements.
    double s_buffer[4 * size];

    MPI_Allgather(s_elems, 4, MPI_DOUBLE, s_buffer, 4, MPI_DOUBLE, comm);

    // Assembling Schur complement matrix: S = F - E0*A0^(-1)*D0 - E1*A1^(-1)*D1 - ...
    // Obs: F = diag(diagElem)
    std::fill_n(diagS, numDecomp, diagElem);

    // First processor only has the last value that is non-zero
    diagS[0] -= s_buffer[3];
    for (int i = 1; i < schurSize; i++)
    {
        diagS[i - 1] -= s_buffer[i * 4];
        upperDiagS[i - 1] -= s_buffer[i * 4 + 1];
        lowerDiagS[i - 1] -= s_buffer[i * 4 + 2];
        diagS[i] -= s_buffer[i * 4 + 3];
    }
    // Last processor only has the first value that is non-zero
    diagS[numDecomp - 2] -= s_buffer[(size - 1) * 4];

    if (rank == 0)
    {
        std::cout << "Schur's complement matrix:" << std::endl;
        for (int i = 0; i < schurSize; i++)
        {
            for (int j = 0; j < schurSize; j++)
            {
                if (i == j)
                    std::cout << diagS[i] << " ";
                else if (i == j - 1)
                    std::cout << lowerDiagS[i] << " ";
                else if (i == j + 1)
                    std::cout << upperDiagS[i] << " ";
                else
                    std::cout << 0.0 << " ";
            }
            std::cout << std::endl;
        }
    }
}

void SchurSolver::updateSchurRhs()
{
    /*
        Every processor will perform its calculations; then we will reduce them and subtract the interface rhs to obtain schurRhs.

        Here we will solve the linear system Ai * yStorage = localRhs, where subrhs is the rhs of the submatrix Ai,
            and then we will multiply the result by the lateral elements of the matrix Ei.
        In the end, we will reduce all these intermediate products and subtract the interface rhs to obtain the Schur rhs.

        We store the intermediate result (A^-1 f) into yStorage since we will need it to compute the solution.

        TODO: rhsInterface and localRhs should already be filled and ready.
            rhsInterface should contain the interface elements of the rhs of the problem.
    */
    /*   REMOVEME FIXME just as a reminder
             lateralElements_D[A_PLUS] = lowerElem;
             lateralElements_D[C_MINUS] = upperElem;
             bottomElements_E[C_PLUS] = upperElem;
             bottomElements_E[A_MINUS] = lowerElem;
    */
    // Store the local intermediate products
    // double *y = new double[dimSubmatrix];             // Solution of the linear system: Ai * y = subrhs
    // double *diag = new double[dimSubmatrix];          // Main diagonal of the Matrix Ai
    // double *upperDiag = new double[dimSubmatrix - 1]; // 1st upper diagonal of the Matrix Ai
    // double *lowerDiag = new double[dimSubmatrix - 1]; // 1st lower diagonal of the Matrix Ai

    /*
    We will allGather() these and sum them together to obtain the Schur rhs.
        Each processor contributes to (at most) 2 values of the Schur rhs.
    */
    double localSchurRhsBuffer[2];

    if (rank == 0)
    {
        /*
            The first processor only needs to fill the first value
        */
        /*std::fill_n(diag, dimSubmatrix, diagElem);
        std::fill_n(upperDiag, dimSubmatrix - 1, upperElem);
        std::fill_n(lowerDiag, dimSubmatrix - 1, lowerElem);*/
        // extractRhsAj(subrhs, 0);

        /*  Solve A0 * yStorage = subrhs
            Store the intermediate result (A^-1 f), we will use it later to compute the solution
        */
        thomasAlgorithm(yStorage, matrixDiag, matrixUpperDiag, matrixLowerDiag, localRhs, dimSubmatrix);

        // Update rhsInterface
        // rhsSchur = rhsInterface - sum(E0 * (A0^-1 * D0))
        // rhsSchur = rhsInterface - E0 * (A0^-1 * D0)

        // Final step: compute E0 * (A0^-1 * D0)
        localSchurRhsBuffer[0] = bottomElements_E[A_MINUS] * yStorage[dimSubmatrix - 1];
        localSchurRhsBuffer[1] = 0.0;
    }
    if (rank == size - 1)
    {
        /*
            The first processor only needs to fill the last value
        */
        // Solve An-1 * y = subrhs
        thomasAlgorithm(yStorage, matrixDiag, matrixUpperDiag, matrixLowerDiag, localRhs, dimSubmatrix);

        localSchurRhsBuffer[0] = 0.0;
        localSchurRhsBuffer[1] = bottomElements_E[C_PLUS] * yStorage[0];
    }
    else
    {
        // Solve Ai * y = subrhs
        thomasAlgorithm(yStorage, matrixDiag, matrixUpperDiag, matrixLowerDiag, localRhs, dimSubmatrix);

        // Update rhsInterface
        localSchurRhsBuffer[0] = bottomElements_E[C_PLUS] * yStorage[0];
        localSchurRhsBuffer[1] = bottomElements_E[A_MINUS] * yStorage[dimSubmatrix - 1];
    }

    /*
    AllGather step

        At this point each processor has its two values. We need to allGather them into the appropriate positions and then reduce them.
        We need a buffer that accomodates 2 elements for each processor; the sum will overlap in one element per operation, eg:
            schurRhs_0   =   b_0 + a_1  == buffer[1] + buffer[2]
            schurRhs_1   =   b_1 + a_2  == buffer[3] + buffer[4]
            schurRhs_2   =   b_2 + a_3  == buffer[5] + buffer[6]
            ...

        where each processor will have
            a_i, b_i as the two elements of the buffer.
        In the end, each processor will own a copy of schurRhs.

    */
    double *intermediateProductsBuffer = new double[numDecomp * 2];

    MPI_Allgather(localSchurRhsBuffer, 2, MPI_DOUBLE, intermediateProductsBuffer, 2, MPI_DOUBLE, comm);

    // Reduce step
    schurRhs[0] = -(intermediateProductsBuffer[1] + intermediateProductsBuffer[2]); // b_0 + a_1
    for (unsigned int i = 1; i < schurSize - 1; ++i)
    {
        schurRhs[i] = -(intermediateProductsBuffer[i * 2 + 1] + intermediateProductsBuffer[i * 2 + 2]);
    }
    schurRhs[schurSize - 1] = -(intermediateProductsBuffer[(schurSize - 1) * 2 + 1] + intermediateProductsBuffer[(schurSize - 1) * 2 + 2]);

    // Final step: compute rhsInterface - E0 * (A0^-1 * D0)
    for (int i = 0; i < schurSize; i++)
    {
        schurRhs[i] += rhsInterface[i];
    }

    // Free memory
    delete[] intermediateProductsBuffer;

    // Print the Schur rhs
    if (rank == 0)
    {
        std::cout << "Schur's rhs:" << std::endl;
        for (int i = 0; i < schurSize; i++)
        {
            std::cout << schurRhs[i] << " ";
        }
        std::cout << std::endl;
    }
}

void SchurSolver::computeSolution()
{
    double *xInt = new double[schurSize]; // Solution of the linear system S * xInt = rhsInterface

    // Solve S * xInt = rhsInterface
    thomasAlgorithm(xInt, diagS, upperDiagS, lowerDiagS, rhsInterface, schurSize);

    // Compute the solution for all the points
    for (int j = 0; j < schurSize; j++)
        solution[(j + 1) * dimSubmatrix + j] = xInt[j];
    for (int j = 0; j < dimSubmatrix; j++)
    {
        solution[j] = yStorage[j];
        solution[j] -= yi[j] * xInt[0];
    }
    for (int i = 1; i < schurSize; i++)
    {
        for (int j = 0; j < dimSubmatrix; j++)
        {
            solution[i * (dimSubmatrix + 1) + j] = yStorage[i * dimSubmatrix + j];
            solution[i * (dimSubmatrix + 1) + j] -= xi[j] * xInt[i - 1];
            solution[i * (dimSubmatrix + 1) + j] -= yi[j] * xInt[i];
        }
    }
    for (int j = 0; j < dimLatestSubmatrix; j++)
    {
        solution[(schurSize) * (dimSubmatrix + 1) + j] = yStorage[(schurSize) * dimSubmatrix + j];
        solution[(schurSize) * (dimSubmatrix + 1) + j] -= xm[j] * xInt[numDecomp - 2];
    }

    // Free memory
    delete[] xInt;
}

void SchurSolver::computeError()
{
    for (int i = 0; i < n; i++)
        error += (solution[i] - exactSolution[i]) * (solution[i] - exactSolution[i]);
    error = std::sqrt(error);
    error = error / n;
}

[[deprecated]] void SchurSolver::thomasAlgorithm(_OUT double *solution, _INOUT double *diag,
                                                 double *upperDiag, double *lowerDiag,
                                                 _INOUT double *rhs, int dim)
{
    double w;

    for (int i = 1; i < dim; i++)
    {
        w = lowerDiag[i - 1] / diag[i - 1];
        diag[i] -= w * upperDiag[i - 1];
        rhs[i] -= w * rhs[i - 1];
    }
    solution[dim - 1] = rhs[dim - 1] / diag[dim - 1];
    for (int i = dim - 2; i >= 0; i--)
    {
        solution[i] = (rhs[i] - upperDiag[i] * solution[i + 1]) / diag[i];
    }
}

void SchurSolver::thomasAlgorithm(_OUT double *solution, const double *const diag,
                                  const double *const upperDiag, const double *const lowerDiag,
                                  const double *const rhs, const int &dim) const
{
    double w;
    double *diagCopy = new double[dim];
    double *rhsCopy = new double[dim];

    // Copy diag and rhs
    std::memcpy(diagCopy, diag, dim);
    std::memcpy(rhsCopy, rhs, dim);

    for (int i = 1; i < dim; i++)
    {
        w = lowerDiag[i - 1] / diagCopy[i - 1];
        diagCopy[i] -= w * upperDiag[i - 1];
        rhsCopy[i] -= w * rhsCopy[i - 1];
    }
    solution[dim - 1] = rhsCopy[dim - 1] / diagCopy[dim - 1];
    for (int i = dim - 2; i >= 0; i--)
    {
        solution[i] = (rhsCopy[i] - upperDiag[i] * solution[i + 1]) / diagCopy[i];
    }

    // Free memory
    delete[] diagCopy;
    delete[] rhsCopy;
}

void SchurSolver::schurSubsystemsSolver(_OUT double x[], _OUT double y[],
                                        double diag, double upperDiag,
                                        double lowerDiag, int dim)
{
    // Wikipedia notation
    double w;
    double *b = new double[dim];
    double *d = new double[dim];

    std::fill_n(b, dim, diag);
    d[0] = lowerDiag;

    for (int i = 1; i < dim; i++)
    {
        w = lowerDiag / b[i - 1];
        b[i] = b[i] - w * upperDiag;
        d[i] = d[i] - w * d[i - 1];
    }
    x[dim - 1] = d[dim - 1] / b[dim - 1];
    y[dim - 1] = upperDiag / b[dim - 1];
    for (int i = dim - 2; i >= 0; i--)
    {
        x[i] = (d[i] - upperDiag * x[i + 1]) / b[i];
        y[i] = (-upperDiag * y[i + 1]) / b[i];
    }
}

void SchurSolver::extractRhsInterface()
{
    int stride = (n - numDecomp + 1) / numDecomp;

    for (int i = 1; i < numDecomp; i++)
    {
        rhsInterface[i - 1] = rhs[i * (stride + 1) - 1];
    }
}

void SchurSolver::extractRhsAj(double rhsAi[], int j)
{
    int stride = (n - numDecomp + 1) / numDecomp, lenghtRhsAi = stride;
    if (j == numDecomp - 1)
        lenghtRhsAi = n - numDecomp + 1 - ((numDecomp - 1) * stride);

    for (int i = 0; i < lenghtRhsAi; i++)
    {
        rhsAi[i] = rhs[i + j * (stride + 1)];
    }
}

void SchurSolver::buildRhs()
{
    // Global part of the rhs: forcing term
    // If we are not the first processor, we also have a left shared boundary
    if (rank > 0)
    {
        const unsigned int global_index = rank - 1;
        localRhs[0] = rhs[global_index];
    }
    // Non-shared part of the rhs
    for (unsigned int i = 1; i < local_n + 1; ++i)
    {
        const unsigned int global_index = rank + i - 1;
        localRhs[i] = rhs[global_index];
    }
    // If we are not the last processor, we also have a right shared boundary
    if (rank < size - 1)
    {
        const unsigned int global_index = rank + local_n;
        localRhs[n + 1] = rhs[global_index];
    }
}

void SchurSolver::buildLateralElements()
{
    if (rank == 0)
    {
        lateralElements_D[A_PLUS] = 0.0;
        lateralElements_D[C_MINUS] = upperElem;
        bottomElements_E[C_PLUS] = 0.0;
        bottomElements_E[A_MINUS] = lowerElem;
    }
    else if (rank == size - 1)
    {
        lateralElements_D[A_PLUS] = lowerElem;
        lateralElements_D[C_MINUS] = 0.0;
        bottomElements_E[C_PLUS] = upperElem;
        bottomElements_E[A_MINUS] = 0.0;
    }
    else
    {
        lateralElements_D[A_PLUS] = lowerElem;
        lateralElements_D[C_MINUS] = upperElem;
        bottomElements_E[C_PLUS] = upperElem;
        bottomElements_E[A_MINUS] = lowerElem;
    }
}

void SchurSolver::buildMatrix()
{
    for (unsigned int i = 0; i < local_n; ++i)
        matrixDiag[i] = diagElem;

    for (unsigned int i = 0; i < local_n; ++i)
        matrixUpperDiag[i] = upperElem;

    for (unsigned int i = 0; i < local_n; ++i)
        matrixLowerDiag[i] = lowerElem;
}
