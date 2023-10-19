#include <iostream>
#include <vector>
#include <cstring>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;


void solve(double* a, double* b, double* c, double* d, int n) {

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


int main(int argc, char** argv) {
	double X = 1.0, dx;
    double T = 1.0, dt;
    int n, N;
    int t = 10; // numero
    //double s;
    cout << "Inserire grandezza"<< endl;
    cin >> n;
    N = pow(n, 3);
    vector<double> a; 
    a.reserve(N); 
    vector<double> b; 
    b.reserve(N);
    vector<double> c; 
    c.reserve(N);
    vector<double> dold; 
    dold.reserve(N);
    vector<double> d; 
    d.reserve(N); 
    vector<double> sol;
    sol.reserve(N);
    vector<double> exactsol;
    exactsol.reserve(N);
    dx = X / (double)(n - 1.0);
    dt = T / t;  
    cout << "Calculated dx = " << dx << endl << endl;
    auto start_s = high_resolution_clock::now();
    
    //Condizione iniziale sul tempo
    for (unsigned int i = 0; i < n; i++)
        for (unsigned int j = 0; j < n; j++)
            for (unsigned int k = 0; k < n; k++)
                {
                     sol[i + j*n + k*n*n] = 0.0;
                     //cold[i + j*n + k*n*n] = 0.0;
                }     




    //ciclo sul tempo
    for (double h = dt; h <= T; h += dt)
    {

            for (unsigned int i = 0; i < n; i++)
              for (unsigned int j = 0; j < n; j++)
                for (unsigned int k = 0; k < n; k++)
                {
                  a[i + j*n + k*n*n] = -dt/(2.0*dx*dx);
                  b[i + j*n + k*n*n] = 1.0 + dt/(dx*dx);
                  c[i + j*n + k*n*n] = -dt/(2.0*dx*dx);
                  a[0] = 0.0;
                  c[N-1] = 0.0;
                  d[i + j*n + k*n*n] = dt*(cos(dt/2 + h) + 3*sin(dt/2 + h))*sin(i*dx)*sin(j*dx)*sin(k*dx) + dt*(- 6*sol[i + j*n + k*n*n]) + sol[i + j*n + k*n*n] 
                  - dt/2.0*(- 2*sol[i + j*n + k*n*n]);

                  if (i != 0) 
                  {
                     d[i + j*n + k*n*n] += dt*(sol[i-1 + j*n + k*n*n]) 
                  - dt/2.0*(sol[i-1 + j*n + k*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale

                  if (j != 0) 
                  {
                     d[i + j*n + k*n*n] += dt*( 
                  sol[i + (j-1)*n + k*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale

                  if (k != 0) 
                  {
                     d[i + j*n + k*n*n] += dt*(sol[i + j*n + (k-1)*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale

                  if (i != n-1) 
                  {
                     d[i + j*n + k*n*n] += dt*(sol[i+1 + j*n + k*n*n])
                  - dt/2.0*(sol[i+1 + j*n + k*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale
                  else 
                     d[i + j*n + k*n*n] += sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(h);

                  if (j != n-1) 
                  {
                     d[i + j*n + k*n*n] += dt*(
                  sol[i + (j+1)*n + k*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale
                  else 
                     d[i + j*n + k*n*n] += sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(h);

                  if (k != n-1) 
                  {
                     d[i + j*n + k*n*n] += dt*(sol[i + j*n + (k+1)*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale
                  else 
                     d[i + j*n + k*n*n] += sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(h);
                }  
            // a[n-1] = -1.0/(dx*dx);
            // b[n-1] = 2.0/(dx*dx);
            // c[n-1] = 0.0;
            // d[n-1] = std::sin(n*dx) + std::sin(1)/(dx*dx);
            // auto stop_A = high_resolution_clock::now();
            // auto duration_A = duration_cast<microseconds>(stop_A - start_A);
            // cout << "Time A assembly = "<< duration_A.count()/(double)n << " microseconds"<<endl;
            // auto start_s = high_resolution_clock::now();
            // for (unsigned int i = 0; i < N; i++)
            //                             cout << c[i] << endl;  
            //primo sistema
            solve(a.data(),b.data(),c.data(),d.data(),N);   // viene modificato c, d ha la soluzione

            //c.swap(cold);
            dold.swap(d);
            //copio la risultante nella soluzione
            //memcpy(d, sol);
            //sol.swap(d);
         
            //secondo sistema
             for (unsigned int i = 0; i < n; i++)
              for (unsigned int j = 0; j < n; j++)
                for (unsigned int k = 0; k < n; k++)
                {
                //   a[i + j*n + k*n*n] = -dt/(2.0*dx*dx);
                //   b[i + j*n + k*n*n] = 1.0 + dt/(dx*dx);
                  c[i + j*n + k*n*n] = -dt/(2.0*dx*dx);
                  //a[0] = 0.0;
                  c[N-1] = 0.0;
                  d[i + j*n + k*n*n] = dold[i + j*n + k*n*n]
                  - dt/2.0*(- 2*sol[i + j*n + k*n*n]);

                  if (j != 0) 
                  {
                     d[i + j*n + k*n*n] -= dt/2.0*(sol[i + (j-1)*n + k*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale
                  if (i == n-1) 
                     d[i + j*n + k*n*n] += sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(h);

                  if (j != n-1) 
                  {
                     d[i + j*n + k*n*n] -= dt/2.0*(sol[i + (j+1)*n + k*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale
                  else 
                     d[i + j*n + k*n*n] += sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(h);

                  if (k == n-1) 
                     d[i + j*n + k*n*n] += sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(h);
                }     

            solve(a.data(),b.data(),c.data(),d.data(),N);   // viene modificato c, d ha la soluzione

            //c.swap(cold);
            dold.swap(d);
            //terzo sistema
            for (unsigned int i = 0; i < n; i++)
              for (unsigned int j = 0; j < n; j++)
                for (unsigned int k = 0; k < n; k++)
                {
                  c[i + j*n + k*n*n] = -dt/(2.0*dx*dx);
                  //a[0] = 0.0;
                  c[N-1] = 0.0;
                  d[i + j*n + k*n*n] = dold[i + j*n + k*n*n]
                  - dt/2.0*(- 2*sol[i + j*n + k*n*n]);

                  if (k != 0) 
                  {
                     d[i + j*n + k*n*n] -= dt/2.0*(sol[i + j*n + (k-1)*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale
                  if (i == n-1) 
                     d[i + j*n + k*n*n] += sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(h);

                  if (k != n-1) 
                  {
                     d[i + j*n + k*n*n] -= dt/2.0*(sol[i + j*n + (k+1)*n*n]);
                  }   //condizioni al contorno spaziali si devono aggiornare ad ognin iterazione temporale
                  else 
                     d[i + j*n + k*n*n] += sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(h);

                  if (j == n-1) 
                     d[i + j*n + k*n*n] += sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(h);
                }   

                 solve(a.data(),b.data(),c.data(),d.data(),N); 

                    for (unsigned int i = 0; i < sol.size(); i++) {
                         sol[i] = sol[i] + d[i];
                        }

    }
    auto stop_s = high_resolution_clock::now();
    //auto stop_total = high_resolution_clock::now();
    auto duration_s = duration_cast<microseconds>(stop_s - start_s);
    cout << "Time solving = "<< duration_s.count() << " microseconds"<<endl;

    for (unsigned int i = 0; i < n; i++)
              for (unsigned int j = 0; j < n; j++)
                for (unsigned int k = 0; k < n; k++)
                {
                   exactsol[i + j*n + k*n*n] = sin(i*dx)*sin(j*dx)*sin(k*dx)*sin(T);         //
                }    
	// for (int i = 0; i < n; i++) {
	// 	cout << d[i] << endl;
	// }
	cout << endl << "n= " << n << endl << "n is not changed hooray !!";
    double err = 0.0;
    double sum = 0.0;
    for(int i=0;i<N;i++){
        sum += (sol[i] - exactsol[i]) * (sol[i] - exactsol[i]);
    }
    err = std::sqrt(sum)/N;
    cout << "Error = "<< err << endl;
    // auto duration_total = duration_cast<microseconds>(stop_total - start_total);
    // cout << "Total time = "<< duration_total.count()/(double)n << "microseconds"<<endl;
    
    // double m;
    // m = (std::log(err) - std::log(3.5e-6))/(std::log(dx) - std::log(1.0/101));
    // cout << "Order = "<< m << endl;
	return 0;
}