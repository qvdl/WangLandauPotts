//
// Created by quinten on 12/13/25.
//

#ifndef SI2530PROJECT_FUNCTIONS_H
#define SI2530PROJECT_FUNCTIONS_H

#include <vector>
#include "qpottsmodel.h"

std::vector<double> computeExactDOS(int dim, size_t L, size_t Q, double J);
std::vector<double> computeExactDOS(QPottsModel& lattice);
std::vector<double> compute_lnVec(const std::vector<double>& vec);

double computeMSE(const std::vector<double> &v1, const std::vector<double> &v2);
double computeTc(double J, int Q);

void print_npArrays(const QPottsModel& model, const std::vector<double>& g_exact, const std::vector<double>& ln_gAbs_wl);
void runWLandPrint(    QPottsModel& model,
    int Q,
    int L,
    int dimension,
    int N );

void writeWLOutputCSV(const std::string& filename, const std::vector<int>& energies, const std::vector<double>& log_densities);
void runWLonceWriteCSV(    QPottsModel& model,
    int Q,
    int L,
    int dimension,
    const std::string& output_file);

#endif //SI2530PROJECT_FUNCTIONS_H