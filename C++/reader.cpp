#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "reader.h"
using namespace std;

StreamingInstabilityData::StreamingInstabilityData(float rho_ice, float rho_sil, float unit_mass, float initial_porosity, string filepath) {
    int i = 1;
    int n_ice = 0;
    int n_sil = 0;
    StreamingInstabilityData::rho_ice = rho_ice;
    StreamingInstabilityData::rho_sil = rho_sil;
    StreamingInstabilityData::unit_mass = unit_mass;
    StreamingInstabilityData::initial_porosity = initial_porosity;
    StreamingInstabilityData::mpar = mpar_code * unit_mass;

    ifstream SIData;
    SIData.open(filepath);

    while (SIData.good()) {
        string line;
        getline(SIData, line, ',');

        if (i <= 3) {
            i++;
            if (i == 3) {
                i++;
            }
            continue;
        } else if (i % 3 == 1) {
            i++;
            continue;
        } else if (i % 3 == 2) {
            n_ice = stoi(line);
            i++;
        } else if (i % 3 == 0) {
            n_sil = stoi(line);
            i++;
            i++;
            double tmp_ice_frac = n_ice / ((n_ice + n_sil) * 1.0);
            double tmp_mass = ((n_ice + n_sil) * 1.0) * mpar;
            double tmp_density = initial_porosity * ((rho_ice * tmp_ice_frac) + (rho_sil * (1 - tmp_ice_frac)));

            ice_fraction.push_back(tmp_ice_frac);
            mass.push_back(tmp_mass);
            density.push_back(tmp_density);
            porosity.push_back(initial_porosity);
        }   
    }
    SIData.close();
};

void StreamingInstabilityData::add_masses(unsigned int n_bins, unsigned int m_per_bins, double min_dens, double max_dens, double min_mass, double max_mass) {
    double dmass = (max_mass - min_mass) / (n_bins * 1.0);
    for (unsigned int i=0; i < n_bins; i++) {
        double current_mass = min_mass;
        for (unsigned int j=0; j < m_per_bins; j++) {
            double new_density = min_dens + ((double) rand() / (RAND_MAX)) * (max_dens - min_dens);
            double new_ice_fraction = (new_density - StreamingInstabilityData::rho_sil * StreamingInstabilityData::initial_porosity) / ((StreamingInstabilityData::rho_ice - StreamingInstabilityData::rho_sil) * StreamingInstabilityData::initial_porosity);
            StreamingInstabilityData::mass.push_back(current_mass);
            StreamingInstabilityData::density.push_back(new_density);
            StreamingInstabilityData::ice_fraction.push_back(new_ice_fraction);
            StreamingInstabilityData::porosity.push_back(StreamingInstabilityData::initial_porosity);

        }
    }

}
