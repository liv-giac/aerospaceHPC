#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(){
    const int Nx = 3;
    const int Ny = 3;
    const int Nz = 3;
    // X is k, J is y, Z is i
    std::vector< int > N ={Nx,Ny,Nz};

    std::vector<std::string> source;
    std::vector<std::string> destination;

    destination.resize(Nx*Ny*Nz);
    source.resize(Nx*Ny*Nz);
    
    // this goes along X so i'm using just the Nx and the Ny
    for(int i = 0 ; i < Nz ; i++)
        for(int j = 0; j < Ny; j++)
            for ( int k = 0 ; k < Nx ; k++){
                source[k + j * Nx + i * Nx * Ny] = ("T_" + std::to_string(k) + "_" + std::to_string(j) + "_" + std::to_string(i)).c_str();            
                cout << source[k + j * Nx + i * Nx * Ny] << endl;
            }

    int h = 0;
    cout << "H = 0 from X to Y" << endl;
    for(int i = 0; i < N[(h+2)%3]; ++i){
        for(int j = 0 ; j < N[(h+1)%3]; ++j){
            for(int k = 0; k < N[h%3]; ++k){
                destination[j + i*N[(h+1)%3] + k * N[(h+1)%3]*N[(h+2)%3]]=source[ k + j*N[h%3] + i*N[h%3]*N[(h+1)%3]];
                cout << destination[j + i*N[(h+1)%3] + k * N[(h+1)%3]*N[(h+2)%3]]<<endl;                
            }
        }
    }

    h = 1;
    cout << "H = 1 from Y to Z" << endl;
    for(int i = 0; i < N[(h+2)%3]; ++i){
        for(int j = 0 ; j < N[(h+1)%3]; ++j){
            for(int k = 0; k < N[h%3]; ++k){
                destination[j + i*N[(h+1)%3] + k * N[(h+1)%3]*N[(h+2)%3]]=source[ k + j*N[h%3] + i*N[h%3]*N[(h+1)%3]];
                cout << destination[j + i*N[(h+1)%3] + k * N[(h+1)%3]*N[(h+2)%3]]<<endl;                
            }
        }
    }

   h = 2;
    cout << "H = 2 from Z to X" << endl;
    for(int i = 0; i < N[(h+2)%3]; ++i){
        for(int j = 0 ; j < N[(h+1)%3]; ++j){
            for(int k = 0; k < N[h%3]; ++k){
                destination[j + i*N[(h+1)%3] + k * N[(h+1)%3]*N[(h+2)%3]]=source[ k + j*N[h%3] + i*N[h%3]*N[(h+1)%3]];
                cout << destination[j + i*N[(h+1)%3] + k * N[(h+1)%3]*N[(h+2)%3]]<<endl;                
            }
        }
    }




    return 0;
}