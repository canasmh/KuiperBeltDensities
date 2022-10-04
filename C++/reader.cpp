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

    StreamingInstabilityData(string filepath) {
        int i = 1;
        int n_ice = 0;
        int n_sil = 0;

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

            }   
        }
        SIData.close();
    }
};

int main() {

    StreamingInstabilityData kbos("../data/si-data.csv");


    

    // ifstream myFile;
    // myFile.open("../data/si-data.csv");

    // while (myFile.good()) {
    //     string line;
    //     getline(myFile, line, ',');
    //     cout << line << endl;
    // }

    // myFile.close();

    return 0;
};