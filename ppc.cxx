
#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>  // ADD THIS - needed for std::vector
#include <cstdlib> // ADD THIS - needed for getenv, atof, etc.

using namespace std;
#include "ini.cxx"

void initialize(){
    m.set();
}

float photon_yield(string loss_type, int energy_gev, float track_length) {
    float rho = 0.9216f; // density of ice in icecube
    float logE = logf(energy_gev);
    float em_cascade_value=5.321*0.910f/rho;  // important for em cascade photons
    float eff_tracl_length = 0;
    float num_photons = 0.0f;
    float ppm = 2000.0f;
    if (loss_type == "amu-") {
        // calcualte photons for muon track
        float additional_track = 1+ max(0.0f, 0.1880f+0.0206f*logE)*0.910f/rho;
        num_photons = track_length>0?track_length*additional_track:0;
    } 
    if (loss_type == "em") {
        num_photons=energy_gev*em_cascade_value;
    }
    return num_photons*ppm;
}

int main(int argc, char* argv[]) {
    // Arguements: loss type str, energy gev, track length m
    cout << "initalizing: " << endl;
    initialize();
    string loss_type= argv[1];
    int energy_gev = stoi(argv[2]);
    float track_length = stof(argv[3]);

    cout << "loss type: " << loss_type << ", has energy: " << energy_gev << ", and track length: " << track_length << endl;
    float photons = photon_yield(loss_type, energy_gev, track_length);
    cout << "photons: " << photons << endl;

    return 0;
}

