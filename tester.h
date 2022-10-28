#pragma once

#include <iostream>
#include <cmath>
#include <array>
#include <string>
#include <sstream>
#include <vector>
#include <fftw3.h>
#include <complex>
#include <complex.h>
#include <algorithm>
#include <fstream>
#include "extra_tools.h"

using namespace std;

double metrik_inf(complex<double> my,complex<double> example);
double metrik_2(complex<double> my,complex<double> example);

void equality(vector<vector<complex<double>>>& in_mas, string file_with_comp_data, string step_name,string metric);
void simple_equality(vector<vector<complex<double>>>& in_mas, string file_with_comp_data, string step_name,string metric);
double relative_er(double one, double two);
