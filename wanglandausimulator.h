//
// Created by quinten on 11/12/25.
//

#ifndef SI2530PROJECT_WANGLANDAUSIMULATOR_H
#define SI2530PROJECT_WANGLANDAUSIMULATOR_H


#include <fstream>
#include <vector>

#include "qpottsmodel.h"
#include "rng.h"

class WangLandauSimulator {
private:
    QPottsModel& lattice;
    int min_energy;
    int max_energy;
    double ln_f;
    std::vector<size_t> histogram;
    std::vector<double> log_densities;
    const int check_interval;

public:
    explicit WangLandauSimulator(QPottsModel &lattice);

    // GET FUNCTIONS
    std::vector<size_t> getHistogram() { return histogram; }
    std::vector<double> getLogDensities() { return log_densities; }

    // HELP FUNCTIONS
    void resetHistogram();
    [[nodiscard]] bool isHistogramFlat(double criterion) const;
    void plotHistogram(const std::string& filename) const;
    void plotLogDensities(const std::string& filename) const;
    [[nodiscard]] int energy_to_index(double energy) const;

    // ALGORITHM IMPLEMENTATION
    void runSimulation(double tolerance, double flatness_criterion);
};

#endif // SI2530PROJECT_WANGLANDAUSIMULATOR_H
