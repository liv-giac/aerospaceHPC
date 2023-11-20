#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(){
    const int Nx = 2;
    const int Ny = 3;
    const int Nz = 3;

    std::vector< int > N ={Nx,Ny,Nz};

    std::vector<std::string> source;
    std::vector<std::string> destination;

    destination.resize(Nx*Ny*Nz);
    source.resize(Nx*Ny*Nz);
    
    // this goes along X
    for(int i = 0 ; i < Nz ; i++)
        for(int j = 0; j < Ny; j++)
            for ( int k = 0 ; k < Nx ; k++){
                source[k + j * Ny + i * Nz * Nz] = ("T_" + std::to_string(k) + "_" + std::to_string(j) + "_" + std::to_string(i)).c_str();            
                cout << source[k + j * Ny + i * Nz * Nz] << endl;
            }
    // Print the starting point
    for (int i = 0; i < source.size(); ++i) {
        std::cout << source[i] << std::endl;
    }

    int h = 1;
    cout << "N[ h %3 ] = "<< N[h%3]<<endl;
    cout << "N[(h+1)%3] = "<<N[(h+1)%3] << endl;
    cout << "N[(h+2)%3] = "<<N[(h+2)%3] << endl;

    for(int i = 0; i < N[(h+2)%3]; ++i){
        for(int j = 0 ; j < N[(h+1)%3]; ++j){
            for(int k = 0; k < N[h%3]; ++k){
                // cout <<"Index destination :"<< j + i*N[(h+1)%3] + k * N[(h+1)%3]*N[(h+2)%3] << endl;
                // cout <<"Index source :" << k + j*N[h%3] +i*N[h%3]*N[h%3] << endl;
                destination[j + i*N[(h+1)%3] + k * N[(h+1)%3]*N[(h+2)%3]]=source[ k + j*N[h%3] + i*N[h%3]*N[(h+1)%3]];
            }
        }
    }
    for(int i = 0 ; i < destination.size(); ++i){
        cout << destination[i] << endl;
    }



    return 0;
}