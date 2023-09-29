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

// implementazione del Thomas algorithm
void solve(double * a, double* b, double * c, double * d, int n) {
    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}

int main(int argc, char** argv)
{
    double X = 1.0, dx;
    int n;
    cout << "Matrix dimension: "<< endl;
    cin >> n;
    
    dx = X / (double)(n + 1);
    cout << "Calculated dx = " << dx << endl << endl;
    SparseMatrix<double> A(n,n);
    A.reserve(VectorXi::Constant(n,3));// per ottimizzare perfomance assembly della matrice (preso da documentazione eigen)
    VectorXd b(n);
    VectorXd x(n);
    
    cout << "======================== EIGEN ==============================="<< endl;
    auto start_A = high_resolution_clock::now();

    cout<< "Assembling matrix"<<endl;
    for (int i=0; i<n; i++) {
    A.coeffRef(i, i) += 2.0 / (dx * dx);
	if(i>0) A.coeffRef(i, i-1) += -1.0 / (dx * dx);
        if(i<n-1) A.coeffRef(i, i+1) += -1.0 / (dx * dx) ;	
    }
    auto stop_A = high_resolution_clock::now();
    
    
    // bcs :
    cout << "Assembling rhs"<<endl;
    cout << "Imposing BCS"<<endl;
    auto start_b = high_resolution_clock::now();

    b[0]=std::sin(1*dx) + 0.0; // primo elemento
    b[n-1]=std::sin((n)*dx) + std::sin(1)/(dx*dx); // ultimo elemento
    cout << "Calculating RHS"<<endl;
    for(int i=1;i<n-1;++i){
        b[i] = std::sin( (i+1) * dx);
    }
    auto stop_b = high_resolution_clock::now();


    cout <<"Solving..."<<endl<<endl;
    auto start_s = high_resolution_clock::now();

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    x = solver.solve(b);

    auto stop_s = high_resolution_clock::now();

    double err = 0.0;
    double sum = 0.0;
    for(int i=0;i<n;i++){
        sum += (x[i] - std::sin((i+1)*dx)) * (x[i] - std::sin((i+1)*dx));
    }
    err = std::sqrt(sum);
    cout << "Error = "<< err << endl<<endl;
    
    auto duration_A = duration_cast<microseconds>(stop_A - start_A);
    cout << "Time A assembly = "<< duration_A.count() << " microseconds"<<endl;
    auto duration_b = duration_cast<microseconds>(stop_b - start_b);
    cout << "Time b assembly = "<< duration_b.count() << " microseconds" << endl;
    auto duration_s = duration_cast<microseconds>(stop_s - start_s);
    cout << "Time elapsed solving = " << duration_s.count() << " microseconds" << endl<<endl;
    auto duration_eigen = duration_A.count() + duration_b.count() + duration_s.count();
    cout << "Total time = "<< duration_eigen  << "microseconds"<<endl<<endl;

    // vettori per algoritmo di Thomas, ciao thomas :)
    cout<< "=======================THOMAS METHOD=========================="<<endl;
    vector<double> a;
    a.reserve(n);
    vector<double> bT;
    bT.reserve(n);
    vector<double> c;
    c.reserve(n);
    vector<double> d;
    d.reserve(n);

    //double a[n],bT[n],c[n],d[n];
    cout << "Assembling support vector for matrix and rhs..."<<endl;
    start_A = high_resolution_clock::now();
    a[0]= 0.0;
    bT[0]= 2.0/(dx*dx);
    c[0]=-1.0/(dx*dx);
    for(int i=1; i< n-1 ; i++){
        a[i]=-1.0/(dx*dx);
        bT[i]=2.0/(dx*dx);
        c[i]=-1.0/(dx*dx);
    }
    a[n-1]= -1.0/(dx*dx);
    bT[n-1]= 2.0/(dx*dx);
    c[n-1]= 0.0;

    d[0]=std::sin(1*dx) + 0.0; 
    d[n-1]=std::sin((n)*dx) + std::sin(1)/(dx*dx);
    for(int i=1;i<n-1;++i){
        d[i] = std::sin( (i+1) * dx);
    }
    stop_A = high_resolution_clock::now();
    
    cout << "Solving the system..."<<endl;
    start_s = high_resolution_clock::now();
    solve(a.data(),bT.data(),c.data(),d.data(),n); // b è il rhs
    stop_s = high_resolution_clock::now();

    err = 0.0;
    sum = 0.0;
    for(int i=0;i<n;i++){
        sum += (d[i] - std::sin((i+1)*dx)) * (d[i] - std::sin((i+1)*dx));
    }
    err = std::sqrt(sum);
    cout << "Error Thomas = "<< err << endl<<endl;

    duration_s = duration_cast<microseconds>(stop_s - start_s);
    duration_A = duration_cast<microseconds>(stop_A - start_A);
    cout<< "Time elapsed assembling = "<< duration_A.count() << " microseconds"<<endl;
    cout<< "Time elapsed solving = "<< duration_s.count() << " microseconds "<<endl;
    auto duration_thomas = duration_A.count() + duration_s.count();
    cout << "Total time = "<< duration_thomas <<" microseconds"<<endl<<endl;

    double percentage = ((duration_eigen / duration_thomas) - 1.0) * 100.0;
    cout << "Thomas method is "<< percentage <<" % faster than eigen"<<endl;




    return 0;

}
