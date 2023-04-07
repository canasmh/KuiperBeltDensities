#ifndef CONSTANTS_H
#define CONSTANTS_H

// List of constants (cgs unless otherwise stated)
double G = 6.67384e-8;
double M_SUN = 2e33;
double M_EARTH = 6e27;
double M_PLUTO = 1.309e25;
double M_CERES = 9.39e23;
double AU_TO_CM = 1.49597871e13;
double SIGMA_MOL = 2e-15;
double AMU = 1.6605402e-24;  // Atomic Mass unit
double R_GAS = 8.314e7;  // The gas constant
double GAMMA = 1.4;  // ratio of specific heats
double GAMMA1 = 1. / GAMMA;  // Inverse of ratio
double M_MOL = 2.3;  // Mean molecular weight
double R_GAS_MU = R_GAS / M_MOL;
double C_P = GAMMA * R_GAS_MU / (GAMMA - 1);
double C_V = C_P / GAMMA;
double K_B = 1.380649e-16;  // Boltzmann Constant
double ALPHA_VISCOSITY = 2./3;
double YRS_TO_SEC = 365.25 * 24 * 60 * 60;

#endif