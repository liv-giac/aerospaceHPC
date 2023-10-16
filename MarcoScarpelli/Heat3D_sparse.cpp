#include <array>
#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <iomanip>
#include <memory>
#include <chrono>
#include <map>

enum exit_codes
{
    SUCCESS = 0,
    FAILURE = 1
};

using matrix = std::map<size_t, double>;
using vector = std::unique_ptr<double[]>;

inline constexpr size_t get_idx_vec(const unsigned int &i, const unsigned int &j, const unsigned int &k, const unsigned int &size)
{
    return i * size * size + j * size + k;
}

inline constexpr size_t get_idx_matr(const unsigned int &row_idx, const unsigned int &col_idx, const unsigned int &size)
{
    return row_idx * size + col_idx;
}

/// @brief Matrix-vector multiplication
/// @param A Matrix
/// @param b Vector
/// @param result Result vector
/// @param size Size of matrix and vector
void matr_vec_mult(const matrix &A, const vector &b, vector &result, const unsigned int &size)
{
#pragma omp parallel for
    for (unsigned int i = 0; i < size; ++i)
    {
        result[i] = 0.0;
    }

    for (const auto &it : A)
    {
        const unsigned int row_idx = it.first / size;
        const unsigned int col_idx = it.first % size;

        result[row_idx] += it.second * b[col_idx];
    }
}

inline double get_scalar_prod(const vector &v1, const vector &v2, const unsigned int &size)
{
    double prod = 0.0;

    for (unsigned int i = 0; i < size; ++i)
    {
        prod += v1[i] * v2[i];
    }
    return prod;
}

inline double get_norm(const vector &v, const unsigned int &size)
{
    return std::sqrt(get_scalar_prod(v, v, size));
}

inline double get_l2_norm(const vector &v, const unsigned int &size)
{
    return get_norm(v, size) / std::sqrt(size);
}

inline void vec_copy(vector &src, vector &dst, const unsigned int &size)
{
    for (unsigned int i = 0; i < size; ++i)
    {
        dst[i] = src[i];
    }
}

unsigned int solve_system(const matrix &A, const vector &b, vector &x, const unsigned int &size, const double &tol = 1e-6, const unsigned int &max_iter = 1000)
{
    vector resid_vector(new double[size]{0});
    vector p_vector(new double[size]{0});
    vector buffer(new double[size]{0});

    unsigned int iter = 0;

    // Init residual
    matr_vec_mult(A, x, buffer, size);

#pragma omp parallel for
    for (unsigned int i = 0; i < size; ++i)
    {
        resid_vector[i] = b[i] - buffer[i];
    }

    // Init p
    vec_copy(resid_vector, p_vector, size);

    double prev_resid_scalar_prod = get_scalar_prod(resid_vector, resid_vector, size);

    while (std::sqrt(prev_resid_scalar_prod) > tol && iter <= max_iter)
    {
        matr_vec_mult(A, p_vector, buffer, size); // buffer = A * p

        const double alpha = prev_resid_scalar_prod / get_scalar_prod(p_vector, buffer, size); // alpha = <r, r> / <p, A * p>

#pragma omp parallel for
        for (unsigned int i = 0; i < size; ++i)
        {
            x[i] += alpha * p_vector[i];
            resid_vector[i] -= alpha * buffer[i];
        }

        /*
            In theory we should exit here if the residual is small enough, but we will continue
            to avoid using an if statement; the while condition will check it.
        */

        // Scalar product here is on the new residual
        const double beta = get_scalar_prod(resid_vector, resid_vector, size) / prev_resid_scalar_prod; // beta = <r_new, r_new> / <r_old, r_old>

#pragma omp parallel for
        for (unsigned int i = 0; i < size; ++i)
        {
            p_vector[i] = resid_vector[i] + beta * p_vector[i];
        }

        prev_resid_scalar_prod = get_scalar_prod(resid_vector, resid_vector, size);

        ++iter;
    }
    std::cout << "Done in " << iter << " iterations, final residual: " << get_norm(resid_vector, size) << std::endl;

    if (get_norm(resid_vector, size) > tol || iter == max_iter)
        return exit_codes::FAILURE;

    return exit_codes::SUCCESS;
}

int main()
{
    const unsigned int spaces[] = {30};
    const unsigned int times[] = {10};

    const double length = 1.0; // [m]

    constexpr size_t iterations = std::end(spaces) - std::begin(spaces);

    for (unsigned int idx = 0; idx < iterations; ++idx)
    {
        const unsigned int space_steps = spaces[idx];
        const unsigned int time_steps = times[idx];

        const double deltaX = length / ((double)(space_steps + 1)); // [m]
        const double deltaT = 1.0 / ((double)(time_steps));         // [s]

        const size_t num_elements = space_steps * space_steps * space_steps;
        matrix matrix;

        vector solution_vec(new double[num_elements]{0});
        vector rhs_vec(new double[num_elements]{0});

        // Assign coefficients

        const double minus_one_over_2_deltaX_sqr = -1.0 / (2.0 * deltaX * deltaX);
        const double three_over_deltaX_sqr = 3.0 / (deltaX * deltaX);
        const double diag_coeff = 1.0 / deltaT + three_over_deltaX_sqr;

        // Init rhs vector
        for (unsigned int i = 0; i < space_steps; ++i)
        {
            for (unsigned int j = 0; j < space_steps; ++j)
            {
                for (unsigned int k = 0; k < space_steps; ++k)
                {
                    rhs_vec[get_idx_vec(i, j, k, space_steps)] = sin(deltaX * (i + 1)) * sin(deltaX * (j + 1)) * sin(deltaX * (k + 1));
                }
            }
        }

        for (unsigned int time_step = 0; time_step < time_steps; ++time_step)
        {
            std::cout << "Time step: " << time_step << std::endl;
// Init rhs vector
#pragma omp parallel for
            for (unsigned int i = 0; i < space_steps; ++i)
            {
                for (unsigned int j = 0; j < space_steps; ++j)
                {
                    for (unsigned int k = 0; k < space_steps; ++k)
                    {
                        if (i > 0)
                            rhs_vec[get_idx_vec(i, j, k, space_steps)] -= minus_one_over_2_deltaX_sqr * solution_vec[get_idx_vec(i - 1, j, k, space_steps)];
                        if (j > 0)
                            rhs_vec[get_idx_vec(i, j, k, space_steps)] -= minus_one_over_2_deltaX_sqr * solution_vec[get_idx_vec(i, j - 1, k, space_steps)];
                        if (k > 0)
                            rhs_vec[get_idx_vec(i, j, k, space_steps)] -= minus_one_over_2_deltaX_sqr * solution_vec[get_idx_vec(i, j, k - 1, space_steps)];
                        if (i < space_steps - 1)
                            rhs_vec[get_idx_vec(i, j, k, space_steps)] -= minus_one_over_2_deltaX_sqr * solution_vec[get_idx_vec(i + 1, j, k, space_steps)];
                        if (j < space_steps - 1)
                            rhs_vec[get_idx_vec(i, j, k, space_steps)] -= minus_one_over_2_deltaX_sqr * solution_vec[get_idx_vec(i, j + 1, k, space_steps)];
                        if (k < space_steps - 1)
                            rhs_vec[get_idx_vec(i, j, k, space_steps)] -= minus_one_over_2_deltaX_sqr * solution_vec[get_idx_vec(i, j, k + 1, space_steps)];

                        rhs_vec[get_idx_vec(i, j, k, space_steps)] += (1.0 / deltaT - three_over_deltaX_sqr) * solution_vec[get_idx_vec(i, j, k, space_steps)];
                    }
                }
            }
            // Init matrix

            // Middle section
            //  TODO
            for (unsigned int i = 0; i < space_steps; ++i)
            {
                for (unsigned int j = 0; j < space_steps; ++j)
                {
                    for (unsigned int k = 0; k < space_steps; ++k)
                    {
                        if (i > 0)
                            matrix[get_idx_matr(get_idx_vec(i - 1, j, k, space_steps), get_idx_vec(i, j, k, space_steps), num_elements)] += minus_one_over_2_deltaX_sqr;
                        if (j > 0)
                            matrix[get_idx_matr(get_idx_vec(i, j - 1, k, space_steps), get_idx_vec(i, j, k, space_steps), num_elements)] += minus_one_over_2_deltaX_sqr;
                        if (k > 0)
                            matrix[get_idx_matr(get_idx_vec(i, j, k - 1, space_steps), get_idx_vec(i, j, k, space_steps), num_elements)] += minus_one_over_2_deltaX_sqr;
                        if (i < space_steps - 1)
                            matrix[get_idx_matr(get_idx_vec(i + 1, j, k, space_steps), get_idx_vec(i, j, k, space_steps), num_elements)] += minus_one_over_2_deltaX_sqr;
                        if (j < space_steps - 1)
                            matrix[get_idx_matr(get_idx_vec(i, j + 1, k, space_steps), get_idx_vec(i, j, k, space_steps), num_elements)] += minus_one_over_2_deltaX_sqr;
                        if (k < space_steps - 1)
                            matrix[get_idx_matr(get_idx_vec(i, j, k + 1, space_steps), get_idx_vec(i, j, k, space_steps), num_elements)] += minus_one_over_2_deltaX_sqr;

                        matrix[get_idx_matr(get_idx_vec(i, j, k, space_steps), get_idx_vec(i, j, k, space_steps), num_elements)] += diag_coeff;
                    }
                }
            }

            if (solve_system(matrix, rhs_vec, solution_vec, num_elements) == exit_codes::FAILURE)
            {
                std::cout << "ERROR: Solution not found!" << std::endl;
                return exit_codes::FAILURE;
            }

            std::cout << std::endl;
            std::cout << "--------------------------------------------------" << std::endl;

            // Print solution
            /*for (unsigned int i = 0; i < num_elements; ++i)
            {
                std::cout << std::setw(5) << solution_vec[i];
            }
            std::cout << std::endl;*/
        }
    }
}
