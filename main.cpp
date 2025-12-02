#include <iostream>
#include "qpottsmodel.h"
#include "wanglandausimulator.h"


int main() {
    auto potts = QPottsModel(2, 3, 3, 1.0);
    auto wl_simulator = WangLandauSimulator<QPottsModel>(potts);
    wl_simulator.runSimulation(1'000'000, 1e-8, 0.8);
    std::cout << "Histogram is flat? --- " << wl_simulator.isHistogramFlat(0.8);

    const std::string hist_filename = "/home/quinten/PycharmProjects/Computational_Physics/SI2530PROJECT/histogram.csv";
    const std::string density_filename = "/home/quinten/PycharmProjects/Computational_Physics/SI2530PROJECT/densities.csv";
    wl_simulator.plotHistogram(hist_filename);
    wl_simulator.plotDensities(density_filename);
}