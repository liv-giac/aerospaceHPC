#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <chrono>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

int main(int argc, char** argv)
{
    double X = 1.0, dx, T = 1.0,dt;
    vector<int> N{11 , 21 , 41};

    vector<double> errors;
    errors.reserve(N.size());
    
    vector<double> dxs;
    dxs.reserve(N.size());

    //se ho n nodi ho n-1 intervalli
    for  (int  n : N ){
    
    dx = 0.0;
    dt = 0.0;

    cout << "========== N = "<< n << " ============="<<endl<<endl;

    dx = X / (double)(n - 1);
    dt = dx * dx * 0.25;

    dxs.push_back(dx);

    cout << "Calculated dx = " << dx << endl;
    cout << "Calculated dt = " << dt << endl<<endl;
    

    vector<double> T_next;
    T_next.reserve(n);
    vector<double> T_now(n ,0.0); // Temperatura al tempo 0

    int iter = 0;
    double t = 0.0;
    double dt_dx2 = dt / (dx*dx);
    double sincos = 0.0;
    auto start = high_resolution_clock::now();
    for ( t = dt ; t <= T ; t = t + dt ){
        
        //j=0
        T_next[0] = 0.0; 

        sincos = ( sin(t) + cos(t));

        for( int j = 1; j < n-1; j++){
            T_next[j] = T_now[j] + dt_dx2*(T_now[j-1] -2*T_now[j] + T_now[j+1]) + dt * sin(j*dx) * sincos;
        }

        T_next[n-1] = sin(1.0) * sin(t); 

        for ( int i = 0 ; i < n ; i++){
            T_now[i] = T_next[i];
        }
        
        iter++;

    }

    auto stop =high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    cout << "Number of iteration = "<< iter << endl;
    cout << "Absolute time elapsed = " << duration.count() << " microseconds" << endl;
    cout << "Time elapsed = " << duration.count() / ((1.0/dx) * (1.0/dt))  << " microseconds (divided by "<< ((1.0/dx) * (1.0/dt)) <<")"<< endl;

    double err = 0.0;
    double sum = 0.0;
    double exact_sol = 0.0;

    for(int k = 0; k <n ; k++){

        exact_sol= (sin(k*dx)*sin(1));

        sum += (T_now[k] - exact_sol) *  (T_now[k] - exact_sol) ;
    }

    err = sqrt(sum) / (double) n;

    cout << "Error = "<< err << endl<<endl;

    errors.push_back(err);

    }

    double n = 0.0;
    double log1 = 0.0;
    double log2 = 0.0;

    for (int i = 0 ; i < errors.size() - 1 ; i++){


        log1 = log(errors[i]) - log(errors[i+1]);
        log2 = log(dxs[i]) - log (dxs[i+1]);

        n = log1 / log2;

        cout << "Convergence order "<< i <<" "<< i+1 << " = "<< n<<endl;

    }

    return 0;

}
