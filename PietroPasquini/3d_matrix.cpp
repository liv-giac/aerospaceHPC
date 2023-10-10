#include <iostream>
#include <math.h>
#include <Eigen/Dense>

int main()
{
    unsigned int n = 3;
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
    for (unsigned int i = n; i < N; i+=n)
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

    std::cout << mat << std::endl;

    return 0;
}