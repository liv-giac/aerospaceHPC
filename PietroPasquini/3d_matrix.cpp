#include <iostream>
#include <fstream>
#include <string>
#include <ios>
#include <math.h>
#include <Eigen/Dense>
#include <unistd.h>

/**
 * process_mem_usage(double &, double &) - takes two doubles by reference,
 * attempts to read the system-dependent data for a process' virtual memory
 * size and resident set size, and return the results in KB.
 *
 * On failure, returns 0.0, 0.0
 */
void process_mem_usage(double &vm_usage, double &resident_set)
{
    using std::ifstream;
    using std::ios_base;
    using std::string;

    vm_usage = 0.0;
    resident_set = 0.0;

    // 'file' stat seems to give the most reliable results
    //
    ifstream stat_stream("/proc/self/stat", ios_base::in);

    // dummy vars for leading entries in stat that we don't care about
    //
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // the two fields we want
    //
    unsigned long vsize;
    long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >> stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

    stat_stream.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}



int main(int argc, char** argv)
{
    std::ofstream out_file;
    out_file.open("memory_used.csv", std::ios_base::app);
    unsigned int n = std::stoi(argv[1]);
    unsigned int N = std::pow(n, 3);
    double dt = 0.1;
    double dx = 0.1;
    Eigen::MatrixXd mat(N, N);

    // Main diagonal
    mat.setIdentity();
    mat *= 1.0 / dt + 3.0 / std::pow(dx, 2.0);

    // 1 off subdiagonals (x)
    mat.topRightCorner(N - 1, N - 1) += Eigen::MatrixXd::Identity(N - 1, N - 1) * 1.0 / (2.0 * std::pow(dx, 2.0));
    mat.bottomLeftCorner(N - 1, N - 1) += Eigen::MatrixXd::Identity(N - 1, N - 1) * 1.0 / (2.0 * std::pow(dx, 2.0));
    for (unsigned int i = n; i < N; i += n)
    {
        mat(i, i - 1) = 0.0;
        mat(i - 1, i) = 0.0;
    }

    // 3 off diagonals (y)
    mat.topRightCorner(N - n, N - n) += Eigen::MatrixXd::Identity(N - n, N - n) * 1.0 / (2.0 * std::pow(dx, 2.0));
    mat.bottomLeftCorner(N - n, N - n) += Eigen::MatrixXd::Identity(N - n, N - n) * 1.0 / (2.0 * std::pow(dx, 2.0));

    // 9 off diagonals (z)
    mat.topRightCorner(N - std::pow(n, 2), N - std::pow(n, 2)) += Eigen::MatrixXd::Identity(N - std::pow(n, 2), N - std::pow(n, 2)) * 1.0 / (2.0 * std::pow(dx, 2.0));
    mat.bottomLeftCorner(N - std::pow(n, 2), N - std::pow(n, 2)) += Eigen::MatrixXd::Identity(N - std::pow(n, 2), N - std::pow(n, 2)) * 1.0 / (2.0 * std::pow(dx, 2.0));

    // std::cout << mat << std::endl;

    // Memory usage
    double vm, rss;
    process_mem_usage(vm, rss);
    vm /= 1024.0;
    rss /= 1024.0;
    std::cout << "Memory used for matrix of a mesh of size " << n << "x" << n << " is: " << std::endl;
    if (vm <= 1024.0)
    {
        std::cout << "VM: " << vm << " MB ; RSS: " << rss << " MB" << std::endl;
    }
    else
    {
        std::cout << "VM: " << vm / 1024.0 << " GB ; RSS: " << rss / 1024.0 << " GB" << std::endl;
    }

    out_file << n << "," << vm << std::endl;
    
    out_file.close();

    return 0;
}