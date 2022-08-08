#include <iostream>
#include <vector>
#include <complex.h>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include "extra_tools.h"



void RVP_correct(vector<vector<std::complex<double>>> &Raw_data, double fs, double K_r, size_t size_range, size_t size_azimuth){

    double f_r = 2 * K_r/fs;//range frequency

    vector<std::complex<double>> RD(size_range);
    std::complex<double> phs_compensation = exp(-M_PI*(f_r*f_r/K_r)*static_cast<complex<double>>(I));
    fftw_plan plan_f, plan_b;
    for(size_t i = 0; i < size_azimuth;i++){
        plan_f = fftw_plan_dft_1d(size_range, (fftw_complex*) &Raw_data[i],
                                            (fftw_complex*) &RD, FFTW_FORWARD, FFTW_ESTIMATE); //making draft plan

        fftw_execute(plan_f); // Fourier Transform
        RD[i] = RD[i] * phs_compensation;
        plan_f = fftw_plan_dft_1d(size_range, (fftw_complex*) &RD,
                                  (fftw_complex*) &Raw_data[i], FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan_b);
    }
    fftw_destroy_plan(plan_f);
    fftw_destroy_plan(plan_b);
}

void Azimuth_FFT(vector<vector<std::complex<double>>> &Raw_data, size_t size_range, size_t size_azimuth){
    fftw_plan plan_f;
    vector<complex<double>> temp(size_azimuth);
    for(size_t j = 0; j < size_range;j++){
        for(size_t i = 0; i < size_azimuth;i++){
            temp[i] = Raw_data[i][j];
        }
        plan_f = fftw_plan_dft_1d(size_azimuth, (fftw_complex*) &temp,
                                  (fftw_complex*) &temp, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan_f);
        for(size_t i = 0; i < size_azimuth;i++){
            Raw_data[i][j] = temp[i];
        }
    }
    fftw_destroy_plan(plan_f);
}
void RMA(){ //Range migration algorithm or omega-K algorithm
    vector<vector<std::complex<double>>> Raw_data = read_file(false);
    size_t size_azimuth = Raw_data.size();
    size_t size_range = Raw_data[0].size();
    // Sensor parameters (ERS satellite)
    double fs = 18.962468e+6;  //Range Sampling Frequency [Hz]
    double K_r = 4.18989015e+11;     // FM Rate Range Chirp [1/s^2] --> up-chirp
    double tau_p = 37.12e-6;       // Chirp duration [s]
    double V = 7098.0194;                // Effective satellite velocity [m/s]
    double Lambda = 0.05656;            // Length of carrier wave [m]
    double R_0 = 852358.15;              // Range to center of antenna footprint [m]
    double ta = 0.6;                     // Aperture time [s]
    double prf = 1679.902;// Pulse Repitition Frequency [Hz]
    double c = 3.0e8;
    double f_c = 0.0; //Central frequency
    double alpha = 0.0; //squint angle of sight

    //(1) RVP correction
    RVP_correct(Raw_data, fs, K_r, size_range, size_azimuth);

    //(2) Azimuth FFT
    Azimuth_FFT(Raw_data, size_range, size_azimuth);

    //(3) Matching filter
    vector<vector<std::complex<double>>> filter(size_azimuth, vector<std::complex<double>>(size_range));
    vector<double> KR(size_range);
    double  r_c, KX, K_xi_min, K_xi_max, sqrt_KX_KR;

    K_xi_min = 0.0;
    K_xi_max = 0.0;

    vector<double> t = fill_up(-ta/2, ta/2, 1/prf);// time axis in azimuth
    vector<double> tau = fill_up(-tau_p/2, tau_p/2, 1/fs);
    KX = 4*M_PI*f_c * cos(alpha)/c;
    for(size_t i = 0; i < size_azimuth;i++){
        r_c = 0.0;//must depend on azimuth
        for(size_t j = 0;j < size_range;j++){
            KR[j] = 4*M_PI*K_r/c * ((f_c/K_r) + tau[j] - 2*r_c/c);
            if(KR[j] * KR[j] - KX * KX > 0.0){
                sqrt_KX_KR = sqrt(KR[j] * KR[j] - KX * KX);
            }else{
                sqrt_KX_KR = 0.0;
            }
            if(sqrt_KX_KR > K_xi_max){
                K_xi_max = sqrt_KX_KR;
            }

            if(sqrt_KX_KR < K_xi_min){
                K_xi_min = sqrt_KX_KR;
            }
            filter[i][j] = exp(KR[j] * r_c + r_c * sqrt_KX_KR *static_cast<complex<double>>(I));
        }
    }

    //(4) Stolt interpolation
    vector<double> K_xi = fill_up(K_xi_min, K_xi_max, size_range);
    //(5) 2D-IFFT

}

int main() {
    RMA();
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
