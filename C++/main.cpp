#include <vector>
#include <iostream>
#include <string>
#include "reader.h"
#include "constants.h"

using namespace std;

int main() {
    cout << "G= " << G << endl;
    cout << "M_PLUTO= " << M_PLUTO << endl;
    cout << "M_SUN= " << M_SUN << endl;
    cout << "M_EARTH= " << M_EARTH << endl;
    cout << "M_CERES " << M_CERES << endl;

    double density_max = 0.7;
    double density_min = 0.5;

    for (unsigned int i=0; i < 100; i++) {
        double random_num = density_min + ((double) rand() / (RAND_MAX)) * (density_max - density_min);
        cout << random_num << endl;
    }
    // string filePath = "../data/si-data.csv";
    // StreamingInstabilityData kbos(1.0, 3.0, 2.823973078884959e+28, 0.5, filePath);
    // cout << "mass " << "density " << "ice fraction (%) " << "porosity " << endl;
    // for (unsigned int i=0; i < kbos.mass.size(); i++) {

    //     cout << kbos.mass[i] << " " << kbos.density[i] << " " << kbos.ice_fraction[i] << " " << kbos.porosity[i] << endl;
    // };

    return 0;
}