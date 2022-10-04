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

    double npar = 1.536e7;
    double eps = 0.03;
    double total_density = 16749076.820152447;
    double dx = 0.00078125;
    double dy = 0.00078125;
    double dz = 0.00078125;

    double mgas_code = total_density * dx * dy * dz;
    double mpar_code = mgas_code * eps / npar;
    double mpar;

    StreamingInstabilityData(float rho_ice, float rho_sil, float unit_mass, float initial_porosity, string filepath) {
        int i = 1;
        int n_ice = 0;
        int n_sil = 0;
        rho_ice = rho_ice;
        rho_sil = rho_sil;
        unit_mass = unit_mass;
        initial_porosity = initial_porosity;
        mpar = mpar_code * unit_mass;

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
                cout << "tmp_ice_frac= " << tmp_ice_frac << endl;
                cout << "tmp_mass= " << tmp_mass << endl;

                ice_fraction.push_back(tmp_ice_frac);
                mass.push_back(tmp_mass);
                density.push_back(tmp_density);
            }   
        }
        SIData.close();
    }
};

int main() {

    StreamingInstabilityData kbos(1.0, 3.0, 2.823973078884959e+28, 0.5, "../data/si-data.csv");
    cout << "Mass     " << "Density     " << "Ice fraction    %   " << endl;
    for (int i=0; i < kbos.mass.size(); i++) {

        cout << kbos.mass[i] << " " << kbos.density[i] << " " << kbos.ice_fraction[i] * 100 << endl;
    }

    return 0;
};