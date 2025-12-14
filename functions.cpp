//
// Created by quinten on 12/13/25.
//

#include <cmath>
#include "functions.h"
#include "qpottsmodel.h"
#include <algorithm>

std::vector<double> computeExactDOS(const int dim, const size_t L, const size_t Q, const double J) {
    QPottsModel lattice(dim, L, Q, J);
    return computeExactDOS(lattice);
}


std::vector<double> computeExactDOS(QPottsModel& lattice) {
    const auto Q = lattice.getQ();
    const int max_energy = QPottsModel::getMaxEnergy();
    const int min_energy = lattice.getMinEnergy();
    const auto nBins = max_energy-min_energy+1;

    // DOS
    std::vector<double> g(nBins, 0.0);

    // Q ^ L ^ d  different configurations
    const size_t N = lattice.getTotalSites();
    const auto total = static_cast<size_t>(std::pow(Q, N));

    // configuration
    std::vector<int> digits(N, 0);

    for (size_t conf=0; conf<total; ++conf) {
        for (size_t i=0; i<N; ++i) { lattice.setSpin(i, digits[i]); }

        const double current_energy = lattice.computeEnergy();
        const int bin = static_cast<int>(std::llround(current_energy)) - min_energy;
        g[static_cast<size_t>(bin)] += 1.0;

        size_t k=0;
        while (k<N) {
            if (++digits[k] < Q) break;
            digits[k] = 0;
            ++k;
        }
    }

    return g;
}


double compute_MSE(const std::vector<double> &v1, const std::vector<double> &v2) {
    double mse = 0.0;
    int nonempty_bins = 0;
    for (int E=0; E<v1.size(); ++E) {
        const double p_exact = v1[E];
        const double p_wl = v2[E];
        if (p_exact == 0) { continue; }

        mse += (p_exact-p_wl) * (p_exact-p_wl);
        nonempty_bins++;
    }

    return mse / static_cast<double>(nonempty_bins);
}