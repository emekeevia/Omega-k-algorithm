#include <iostream>
#include <vector>
#include "extra_tools.h"



void RVP(vector<vector<std::complex<double>>> &Raw_data, double f_r, double K_r){

}

void RMA(){ //Range migration algorithm or omega-K algorithm
    vector<vector<std::complex<double>>> Raw_data = read_file(false);
    // Sensor parameters (ERS satellite)
    double fs = 18.962468e+6;  //Range Sampling Frequency [Hz]
    double K_r = 4.18989015e+11;     // FM Rate Range Chirp [1/s^2] --> up-chirp
    double tau_p = 37.12e-6;       // Chirp duration [s]
    double V = 7098.0194;                // Effective satellite velocity [m/s]
    double Lambda = 0.05656;            // Length of carrier wave [m]
    double R_0 = 852358.15;              // Range to center of antenna footprint [m]
    double ta = 0.6;                     // Aperture time [s]
    double prf = 1679.902;               // Pulse Repitition Frequency [Hz]

    RVP(Raw_data, fs, K_r);
}

int main() {
    RMA();
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
