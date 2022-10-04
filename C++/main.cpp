#include <vector>
#include <iostream>
#include <string>
#include "reader.h"

using namespace std;

int main() {
    string filePath = "../data/si-data.csv";
    StreamingInstabilityData kbos(1.0, 3.0, 2.823973078884959e+28, 0.5, filePath);
    cout << "mass " << "density " << "ice fraction (%) " << "porosity " << endl;
    for (unsigned int i=0; i < kbos.mass.size(); i++) {

        cout << kbos.mass[i] << " " << kbos.density[i] << " " << kbos.ice_fraction[i] << " " << kbos.porosity[i] << endl;
    };
    
    return 0;
}