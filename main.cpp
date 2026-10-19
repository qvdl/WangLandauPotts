#include <iostream>
#include "qpottsmodel.h"
#include "wanglandausimulator.h"
#include "ThermoState.h"


int main() {
    auto model = QPottsModel(2, 5, 3, 1.0);
    auto wl_simulator = WangLandauSimulator(model);
    wl_simulator.runSimulation(1e-8, 0.8);
    const auto log_d = wl_simulator.getLogDensities();

    const std::string density_filename = "densities.csv";
    wl_simulator.plotLogDensities(density_filename);

    constexpr double temp=2.71;
    ThermoState thermo(wl_simulator.getMinEnergy(), log_d, temp, 1);
}
