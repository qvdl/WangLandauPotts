//
// Created by quinten on 12/13/25.
//

#ifndef SI2530PROJECT_FUNCTIONS_H
#define SI2530PROJECT_FUNCTIONS_H

#include <vector>
#include "qpottsmodel.h"

std::vector<double> computeExactDOS(int dim, size_t L, size_t Q, double J);
std::vector<double> computeExactDOS(QPottsModel& lattice);

double compute_MSE(const std::vector<double> &v1, const std::vector<double> &v2);

#endif //SI2530PROJECT_FUNCTIONS_H