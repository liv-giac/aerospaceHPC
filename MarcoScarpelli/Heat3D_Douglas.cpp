
#include <array>
#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>
#include <iomanip>
#include <memory>
#include <chrono>
#include <map>

/*
///////////////////////////////////////////////////////////////////////////////////////
    VARIOUS DEFINITIONS
*/

double _cur_space_steps;
double _cur_time_steps;
double _cur_deltaX;
double _cur_deltaT;
unsigned int _cur_time_step;

enum exit_codes
{
    SUCCESS = 0,
    FAILURE = 1
};

enum solve_direction
{
    X,
    Y,
    Z
};

using matrix = std::map<size_t, double>;
using vector = std::unique_ptr<double[]>;

/*
///////////////////////////////////////////////////////////////////////////////////////
    VALUES
*/

constexpr double get_exact_solution(const double &x, const double &y, const double &z, const double &t)
{
    return std::sin(x) * std::sin(y) * std::sin(z) * std::sin(t);
}

constexpr double get_forcing_term(const double &x, const double &y, const double &z, const double &t)
{
    return (std::cos(t) + 3 * std::sin(t)) * std::sin(x) * std::sin(y) * std::sin(z);
}

inline double get_diag_coeff()
{
    return 1.0 + _cur_deltaT / (_cur_deltaX * _cur_deltaX);
}

inline double get_off_diag_coeff()
{
    return -_cur_deltaT / (2.0 * _cur_deltaX * _cur_deltaX);
}

/*
///////////////////////////////////////////////////////////////////////////////////////
    INDEXING
*/

constexpr size_t get_idx_vec(const unsigned int &i, const unsigned int &j, const unsigned int &k, const unsigned int &size)
{
    return i * size * size + j * size + k;
}

size_t get_idx_vec(const unsigned int &i, const unsigned int &j, const unsigned int &k, const unsigned int &size, const solve_direction &dir)
{
    switch (dir)
    {
    case Z:
        return i * size * size + j * size + k;
    case X:
        return j * size * size + k * size + i;
    case Y:
        return k * size * size + i * size + j;
    }
}

constexpr size_t get_idx_matr(const unsigned int &row_idx, const unsigned int &col_idx, const unsigned int &size)
{
    return row_idx * size + col_idx;
}

/*
///////////////////////////////////////////////////////////////////////////////////////
    VECTOR-MATRIX OPERATIONS
*/
double get_scalar_prod(const vector &v1, const vector &v2, const unsigned int &size)
{
    double prod = 0.0;
    for (unsigned int i = 0; i < size; ++i)
    {
        prod += v1[i] * v2[i];
    }
    return prod;
}

void vec_copy(vector &src, vector &dst, const unsigned int &size)
{
#pragma omp parallel for
    for (unsigned int i = 0; i < size; ++i)
    {
        dst[i] = src[i];
    }
}

/*
///////////////////////////////////////////////////////////////////////////////////////
    NORM
*/

double get_norm(const vector &v, const unsigned int &size)
{
    return std::sqrt(get_scalar_prod(v, v, size));
}

double get_l2_norm(const vector &v, const unsigned int &size)
{
    return get_norm(v, size) / std::sqrt(size);
}

/*
///////////////////////////////////////////////////////////////////////////////////////
    SOLVER
*/

/**
 * @brief Solve the linear system using Thomas' algorithm.
 *
 * @param mat Input matrix (tridiagonal)
 * @param x Solution vector
 * @param b Right-hand side vector
 *
 * @details We need to solve a tridiagonal system, but this is not the actual main system we need to solve. Depending on the direction parameter (x, y or z) we need to compute the directional derivative along that direction, which means accessing neighbouring elements along that direction only. This implies that we need to access the solution and right-hand side vectors in a way that depends on the direction, while streamlining the accessor logic to avoid code duplication. This is offloaded to the get_idx_vec function, which takes the direction as a parameter and returns the correct index.
 */
void solve_system(vector &x, const vector &d, const solve_direction &dir)
{
    const unsigned int size = _cur_space_steps * _cur_space_steps * _cur_space_steps;

    vector c_new(new double[size - 1]);
    vector d_new(new double[size]);

    const unsigned int __cur_space_steps = _cur_space_steps;

/* We access the rhs and solution vectors using a triple nested loop, where the only direction accessing neighbouring elements is the last one; the indexing function will convert it to the correct index depending on the direction. This means we only need to check for boundary breaking on the last index.
 */
#pragma omp parallel for
    for (unsigned int k = 0; k < __cur_space_steps; ++k)
    {
        for (unsigned int j = 0; j < _cur_space_steps; ++j)
        {
            c_new[get_idx_vec(0, j, k, _cur_space_steps, dir)] = get_off_diag_coeff() / get_diag_coeff();
            d_new[get_idx_vec(0, j, k, _cur_space_steps, dir)] = d[get_idx_vec(0, j, k, _cur_space_steps, dir)] / get_diag_coeff();

            for (unsigned int i = 1; i < _cur_space_steps - 1; ++i)
            {
                c_new[get_idx_vec(i, j, k, _cur_space_steps, dir)] = get_off_diag_coeff() / (get_diag_coeff() - get_off_diag_coeff() * c_new[get_idx_vec(i - 1, j, k, _cur_space_steps, dir)]);
                d_new[get_idx_vec(i, j, k, _cur_space_steps, dir)] = (d[get_idx_vec(i, j, k, _cur_space_steps, dir)] - get_off_diag_coeff() * d_new[get_idx_vec(i - 1, j, k, _cur_space_steps, dir)]) / (get_diag_coeff() - get_off_diag_coeff() * c_new[get_idx_vec(i - 1, j, k, _cur_space_steps, dir)]);
            }
            d_new[get_idx_vec(_cur_space_steps - 1, j, k, _cur_space_steps, dir)] = (d[get_idx_vec(_cur_space_steps - 1, j, k, _cur_space_steps, dir)] - get_off_diag_coeff() * d_new[get_idx_vec(_cur_space_steps - 2, j, k, _cur_space_steps, dir)]) / (get_diag_coeff() - get_off_diag_coeff() * c_new[get_idx_vec(_cur_space_steps - 2, j, k, _cur_space_steps, dir)]);
        }
    }

#pragma omp parallel for
    for (unsigned int k = __cur_space_steps; k > 0; --k)
    {
        for (unsigned int j = _cur_space_steps; j > 0; --j)
        {
            // Indices are like so to avoid underflow
            x[get_idx_vec(_cur_space_steps - 1, j - 1, k - 1, _cur_space_steps, dir)] = d_new[get_idx_vec(_cur_space_steps - 1, j - 1, k - 1, _cur_space_steps, dir)];

            for (unsigned int i = _cur_space_steps - 1; i > 0; --i)
            {
                // x[i] = rhs'[i] - c'[i]x[i+1]
                x[get_idx_vec(i - 1, j - 1, k - 1, _cur_space_steps, dir)] = d_new[get_idx_vec(i - 1, j - 1, k - 1, _cur_space_steps, dir)] - c_new[get_idx_vec(i - 1, j - 1, k - 1, _cur_space_steps, dir)] * x[get_idx_vec(i, j - 1, k - 1, _cur_space_steps, dir)];
            }
        }
    }
}

/*
///////////////////////////////////////////////////////////////////////////////////////
    INITIALISATIONS
*/

void build_rhs(vector &rhs, const vector &cur_solution)
{
    const unsigned int __cur_space_steps = _cur_space_steps;

#pragma omp parallel for
    for (unsigned int i = 0; i < __cur_space_steps; ++i)
    {
        for (unsigned int j = 0; j < _cur_space_steps; ++j)
        {
            for (unsigned int k = 0; k < _cur_space_steps; ++k)
            {
                // dT * (f(n+.5) + (dxx+dyy+dzz)u)
                rhs[get_idx_vec(i, j, k, _cur_space_steps)] =
                    (_cur_deltaT / (_cur_deltaX * _cur_deltaX)) * ((i < _cur_space_steps - 1 ? cur_solution[get_idx_vec(i + 1, j, k, _cur_space_steps)] : 0.0) +
                                                                   (i > 0 ? cur_solution[get_idx_vec(i - 1, j, k, _cur_space_steps)] : 0.0) +
                                                                   (j < _cur_space_steps - 1 ? cur_solution[get_idx_vec(i, j + 1, k, _cur_space_steps)] : 0.0) +
                                                                   (j > 0 ? cur_solution[get_idx_vec(i, j - 1, k, _cur_space_steps)] : 0.0) +
                                                                   (k < _cur_space_steps - 1 ? cur_solution[get_idx_vec(i, j, k + 1, _cur_space_steps)] : 0.0) +
                                                                   (k > 0 ? cur_solution[get_idx_vec(i, j, k - 1, _cur_space_steps)] : 0.0) - 6.0 * cur_solution[get_idx_vec(i, j, k, _cur_space_steps)]);

                // Add forcing term to boundary elements only
                if (i == 0 || i == _cur_space_steps - 1 || j == 0 || j == _cur_space_steps - 1 || k == 0 || k == _cur_space_steps - 1)
                {
                    rhs[get_idx_vec(i, j, k, _cur_space_steps)] += _cur_deltaT * get_forcing_term(i * _cur_deltaX, j * _cur_deltaX, k * _cur_deltaX, _cur_deltaT * (double)(2 * _cur_time_step + 1) / 2.0);
                }
            }
        }
    }
}

inline void apply_bcs(vector &rhs, const solve_direction &dir)
{
    const unsigned int __cur_space_steps = _cur_space_steps;

    double x = 1.0, y = 1.0, z = 1.0;

    switch (dir)
    {
    case X:
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < __cur_space_steps; ++k)
        {
            for (unsigned int j = 0; j < _cur_space_steps; ++j)
            {
                y = j * _cur_deltaX;
                z = k * _cur_deltaX;
                rhs[get_idx_vec(_cur_space_steps - 1, j, k, _cur_space_steps, dir)] += (get_exact_solution(x, y, z, _cur_deltaT * (double)(_cur_time_step + 1)) - get_exact_solution(x, y, z, _cur_deltaT * (double)(_cur_time_step))) * _cur_deltaT / (_cur_deltaX * _cur_deltaX);
            }
        }
        break;
    }
    case Y:
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < __cur_space_steps; ++k)
        {
            for (unsigned int j = 0; j < _cur_space_steps; ++j)
            {
                x = j * _cur_deltaX;
                z = k * _cur_deltaX;
                rhs[get_idx_vec(_cur_space_steps - 1, j, k, _cur_space_steps, dir)] += (get_exact_solution(x, y, z, _cur_deltaT * (double)(_cur_time_step + 1)) - get_exact_solution(x, y, z, _cur_deltaT * (double)(_cur_time_step))) * _cur_deltaT / (_cur_deltaX * _cur_deltaX);
            }
        }
        break;
    }
    case Z:
    {
#pragma omp parallel for
        for (unsigned int k = 0; k < __cur_space_steps; ++k)
        {
            for (unsigned int j = 0; j < _cur_space_steps; ++j)
            {
                x = j * _cur_deltaX;
                y = k * _cur_deltaX;
                rhs[get_idx_vec(_cur_space_steps - 1, j, k, _cur_space_steps, dir)] += (get_exact_solution(x, y, z, _cur_deltaT * (double)(_cur_time_step + 1)) - get_exact_solution(x, y, z, _cur_deltaT * (double)(_cur_time_step))) * _cur_deltaT / (_cur_deltaX * _cur_deltaX);
            }
        }
        break;
    }
    }
}

int main()
{
    const unsigned int spaces[] = {100};
    const unsigned int times[] = {10};

    const double length = 1.0; // [m]

    constexpr size_t iterations = std::end(spaces) - std::begin(spaces);

    std::vector<double> errors;
    errors.reserve(iterations);
    std::vector<double> deltaXs;
    deltaXs.reserve(iterations);
    std::vector<double> deltaTs;
    deltaTs.reserve(iterations);
    std::vector<unsigned long long> solve_times;
    solve_times.reserve(iterations);
    std::vector<unsigned long long> setup_times;
    setup_times.reserve(iterations);

    for (unsigned int idx = 0; idx < iterations; ++idx)
    {
        auto start_setup = std::chrono::high_resolution_clock::now();
        const unsigned int space_steps = spaces[idx];
        const unsigned int time_steps = times[idx];

        _cur_space_steps = space_steps;
        _cur_time_steps = time_steps;

        const double deltaX = length / ((double)(space_steps + 1)); // [m]
        const double deltaT = 1.0 / ((double)(time_steps));         // [s]

        _cur_deltaX = deltaX;
        _cur_deltaT = deltaT;

        deltaTs.push_back(deltaT);
        deltaXs.push_back(deltaX);

        std::cout << "Space steps: " << space_steps << ", time steps: " << time_steps << std::endl;

        const size_t num_elements = space_steps * space_steps * space_steps;
        matrix matrix;

        vector solution_vec(new double[num_elements]{0});
        vector solution_vec_old(new double[num_elements]{0});
        vector delta_solution(new double[num_elements]{0});
        vector rhs_vec(new double[num_elements]);

        auto end_setup = std::chrono::high_resolution_clock::now();
        auto elapsed_setup = std::chrono::duration_cast<std::chrono::nanoseconds>(end_setup - start_setup).count();
        setup_times.push_back(elapsed_setup);

        // std::cout << "Program initialised" << std::endl;

        // Assign coefficients

        auto start = std::chrono::high_resolution_clock::now();

        for (_cur_time_step = 0; _cur_time_step <= time_steps; ++_cur_time_step)
        {
            // Assign current solution vector to old
            vec_copy(solution_vec, solution_vec_old, num_elements);

            // Build initial rhs
            build_rhs(rhs_vec, solution_vec);

            // Solve 1
            /*
                Equation is:
                    (1- deltaT/2 * laplace_x)(solution_star - solution) = rhs
            */
            apply_bcs(delta_solution, X);
            solve_system(delta_solution, rhs_vec, X);

            // Solve 2
            /*
                Equation is:
                    (1- deltaT/2 * laplace_y)(solution_star - solution) = rhs
            */
            apply_bcs(delta_solution, Y);
            solve_system(delta_solution, delta_solution, Y);

            // Solve 3
            /*
                Equation is:
                    (1- deltaT/2 * laplace_z)(solution_star - solution) = rhs
            */

            apply_bcs(delta_solution, Z);
            solve_system(delta_solution, delta_solution, Z);

// Get solution
#pragma omp parallel for
            for (unsigned int i = 0; i < num_elements; ++i)
            {
                solution_vec[i] = solution_vec_old[i] + delta_solution[i];
            }

            std::cout << "Solved time step " << _cur_time_step << std::endl;
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        std::cout << "-------------------------------------------" << std::endl;

        // Compute error norm
        double error = 0.0;
#pragma omp parallel for
        for (unsigned int i = 0; i < space_steps; ++i)
        {
            for (unsigned int j = 0; j < space_steps; ++j)
            {
                for (unsigned int k = 0; k < space_steps; ++k)
                {
                    const double exact = get_exact_solution(i * deltaX, j * deltaX, k * deltaX, 1.0);
                    const double approx = solution_vec[get_idx_vec(i, j, k, space_steps)];
                    error += std::pow(exact - approx, 2);
                }
            }
        }
        error = std::sqrt(error) / num_elements;
        errors.push_back(error);

        std::cout << "Error: " << error << std::endl;

        std::cout << "===========================================" << std::endl;

        solve_times.push_back(elapsed);
    }

    for (unsigned int i = 0; i < iterations; ++i)
    {
        std::cout << "Space steps: " << spaces[i] << ", time steps: " << times[i] << std::endl
                  << "Setup time: " << setup_times[i] << "ns , solve time " << solve_times[i] << "ns" << std::endl;
        // Also print time in seconds
        std::cout << "Setup time: " << setup_times[i] / 1e9 << "s , solve time " << solve_times[i] / 1e9 << "s" << std::endl;
        std::cout << "Solve time per element, per iteration:" << solve_times[i]/(spaces[i]*spaces[i]*spaces[i]*times[i]) << "ns" << std::endl;
        std::cout << "Solve time per element, per iteration:" << solve_times[i]/(spaces[i]*spaces[i]*spaces[i]*times[i])/1e9 << "s" << std::endl;
    }
}
