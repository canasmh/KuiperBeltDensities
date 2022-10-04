#ifndef READER_H
#define READER_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class StreamingInstabilityData
{
public:
    float rho_ice;
    float rho_sil;
    double unit_mass;
    float initial_porosity;
    vector<double> mass;
    vector<double> density;
    vector<double> ice_fraction;
    vector<double> porosity;

    double npar = 1.536e7;
    double eps = 0.03;
    double total_density = 16749076.820152447;
    double dx = 0.00078125;
    double dy = 0.00078125;
    double dz = 0.00078125;

    double mgas_code = total_density * dx * dy * dz;
    double mpar_code = mgas_code * eps / npar;
    double mpar;

    StreamingInstabilityData(float rho_ice, float rho_sil, float unit_mass, float initial_porosity, string filepath);
};

#endif