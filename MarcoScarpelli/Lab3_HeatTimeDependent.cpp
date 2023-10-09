#include <array>
#include <iostream>
#include <cmath>
#include <chrono>
#include <vector>

// inline double exact_solution(const double &x)
//{
//     return std::sin(x);
// }

// For cache friendliness
struct compute_element
{
    double solution_elem_prev;
    double solution_elem;
};

int main()
{
    // Heat equation solver
    std::vector<double> errors;
    std::vector<double> deltaXs;
    std::vector<double> deltaTs;
    std::vector<unsigned long long> time_per_iteration;

    auto setup_start = std::chrono::high_resolution_clock::now();

    // unsigned int spaces[] = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 100000000};
    unsigned int spaces[] = {10, 20, 40, 100};
    unsigned int times[] = {400, 1600, 6400, 40000};

    for (unsigned int idx = 0; idx < 4; ++idx)
    {
        // Space settings
        const unsigned int space_steps = spaces[idx]; // 100000000;
        const double length = 1.0;                    // [m]

        const unsigned int time_steps = times[idx];
        const double max_time = 1.0; // [m]

        // Constants calculation
        const double deltaX = length / ((double)(space_steps - 1)); // [m]
        const double deltaT = max_time / ((double)(time_steps));    // [s]

        deltaXs.push_back(deltaX);
        deltaTs.push_back(deltaT);

        compute_element *elements = new compute_element[space_steps];

        // Initialise solution at t = 0
        for (unsigned int cur_space_elem = 0; cur_space_elem < space_steps; ++cur_space_elem)
        {
            elements[cur_space_elem].solution_elem_prev = 0.0;
        }

        for (unsigned int time_step = 0; time_step <= time_steps; ++time_step)
        {
            auto setup_start = std::chrono::high_resolution_clock::now();

            const double cur_time = (double)time_step * deltaT;
            const double deltaT_over_deltaX_sq = deltaT / (deltaX * deltaX);

            // Boundary conditions
            elements[0].solution_elem = 0.0;
            elements[space_steps - 1].solution_elem = std::sin(1.0) * std::sin((double)(time_step)*deltaT);

            auto setup_end = std::chrono::high_resolution_clock::now();

            /*
            ---------------
                SOLVING
            ---------------
            */
            auto solve_start = std::chrono::high_resolution_clock::now();

            // elements[0].solution_elem = elements[0].solution_elem_prev + deltaT_over_deltaX_sq * (-2.0 * elements[0].solution_elem_prev + elements[1].solution_elem_prev) + deltaT * std::sin((double)(space_steps - 1) * deltaX) * (std::sin(cur_time) + std::cos(cur_time));
            // elements[0].solution_elem = 0.0;
            for (unsigned int cur_space_elem = 1; cur_space_elem < space_steps - 1; ++cur_space_elem)
            {
                const double cur_space = (double)cur_space_elem * deltaX;

                elements[cur_space_elem].solution_elem = elements[cur_space_elem].solution_elem_prev + deltaT_over_deltaX_sq * (elements[cur_space_elem - 1].solution_elem_prev - 2.0 * elements[cur_space_elem].solution_elem_prev + elements[cur_space_elem + 1].solution_elem_prev) + deltaT * std::sin(cur_space) * (std::sin(cur_time) + std::cos(cur_time));
            }
            // elements[space_steps - 1].solution_elem = std::sin(1.0) * std::sin(cur_time);

            auto solve_end = std::chrono::high_resolution_clock::now();
            unsigned long long solve_time = std::chrono::duration_cast<std::chrono::nanoseconds>(solve_end - solve_start).count() / space_steps;

            for (unsigned int cur_space_elem = 0; cur_space_elem < space_steps; ++cur_space_elem)
            {
                elements[cur_space_elem].solution_elem_prev = elements[cur_space_elem].solution_elem;
            }

            time_per_iteration.push_back(solve_time);
        }
        /*
        ---------------
            ERROR CALCULATION
        ---------------
        */

        std::cout << "Calculating error..." << std::endl;

        auto error_start = std::chrono::high_resolution_clock::now();

        double error = 0.0;

        // L2 error norm
        for (unsigned int cur_space_elem = 0; cur_space_elem < space_steps; ++cur_space_elem)
        {
            double cur_t = (double)(time_steps)*deltaT;

            error += std::pow(elements[cur_space_elem].solution_elem - (std::sin((double)(cur_space_elem)*deltaX) * (std::sin(cur_t))), 2);
            // std::cout << std::pow(elements[cur_space_elem].solution_elem - std::sin((double)(cur_space_elem + 1) * deltaX), 2) << std::endl;
        }
        error = std::sqrt(error);
        error /= (double)space_steps;

        auto error_end = std::chrono::high_resolution_clock::now();

        std::cout << "Spaces: " << space_steps << std::endl;
        std::cout << "Time steps: " << time_steps << std::endl;

        std::cout << "L2 error norm: " << error << std::endl;

        // unsigned long long setup_time = std::chrono::duration_cast<std::chrono::nanoseconds>(setup_end - setup_start).count() / space_steps;
        unsigned long long error_time = std::chrono::duration_cast<std::chrono::nanoseconds>(error_end - error_start).count() / space_steps;
        unsigned long long total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(error_end - setup_start).count() / space_steps;

        unsigned long long solve_time = 0;
        for (const unsigned long long &t : time_per_iteration)
        {
            solve_time += t;
        }

        std::cout << "--------------------------------" << std::endl;

        std::cout << "Total solve time: " << solve_time << " ns = " << (double)solve_time / 1.e9 << "s" << std::endl;

        solve_time /= time_steps;

        // std::cout << "Setup time per cell: " << setup_time << " ns = " << (double)setup_time / 1.e9 << "s" << std::endl;
        std::cout << "Average solve time per iteration: " << solve_time << " ns = " << (double)solve_time / 1.e9 << "s" << std::endl;
        std::cout << "Average solve time per iteration, per cell: " << solve_time << " ns = " << (double)solve_time / ((double)space_steps * 1.e9) << "s" << std::endl;
        std::cout << "Error time per cell: " << error_time << " ns = " << (double)error_time / 1.e9 << "s" << std::endl;
        std::cout << "Total time per cell: " << total_time << " ns = " << (double)total_time / 1.e9 << "s" << std::endl;

        errors.push_back(error);

        delete[] elements;

        std::cout << std::endl;
    }

    std::cout << "deltaXs: " << std::endl;
    for (const double &d : deltaXs)
    {
        std::cout << d << ",";
    }
    std::cout << std::endl;

    std::cout << "deltaTs: " << std::endl;
    for (const double &d : deltaTs)
    {
        std::cout << d << ",";
    }
    std::cout << std::endl;

    std::cout << "Error norms: " << std::endl;
    for (const double &d : errors)
    {
        std::cout << d << ",";
    }
    std::cout << std::endl;

    std::cout << "Error order: " << std::endl;
    for (size_t i = 0; i < errors.size() - 1; i++)
    {
        std::cout << (std::log(errors[i + 1]) - std::log(errors[i])) / (std::log(deltaXs[i + 1]) - std::log(deltaXs[i])) << std::endl;
    }
}
