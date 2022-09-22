#include "tester.h"

void equality(vector<vector<complex<double>>>& in_mas, string file_with_comp_data, string step_name){
    vector<vector<complex<double>>> sample_mas = read_file(false, file_with_comp_data);
    size_t azimuth_size = in_mas.size();
    size_t range_size = in_mas[0].size();
    double epsilon = 0.01;
    double in, sample;
    bool flag;

    for(size_t i = 0; i < azimuth_size;i++){
        for(size_t j = 0; j < range_size;j++){
            in = abs(in_mas[i][j]);
            sample = abs(sample_mas[i][j]);
            if(in != 0.0 && sample != 0.0){
                if(abs(in - sample)/max(in, sample) >= epsilon){
                    flag = true;
                }
            }else {
                if(abs(in - sample)>= epsilon){
                    flag = true;
                }
            }
        }
    }
    if(flag){
        cerr << step_name << ": error";
    }else{
        cerr << step_name << ": OK";
    }
}