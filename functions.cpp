//
// Created by quinten on 12/13/25.
//

#include <cmath>
#include "functions.h"
#include "qpottsmodel.h"
#include <iostream>
#include <chrono>
#include <fstream>
#include "wanglandausimulator.h"

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

    auto t0 = std::chrono::high_resolution_clock::now();
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
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
    std::cout << "Checking " << std::pow(Q, N) << " configurations took " << duration << " ms\n\n";

    return g;
}

std::vector<double> compute_lnVec(const std::vector<double>& vec) {
    std::vector<double> lnVec(vec.size(), 0.0);

    for (size_t i=0; i<vec.size(); ++i) {
        if (vec[i]>0) { lnVec[i] = std::log(vec[i]); }
    }

    return lnVec;
}

double computeMSE(const std::vector<double> &v1, const std::vector<double> &v2) {
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

double computeTc(const double J, const int Q) {
    return J / std::log( 1 + std::sqrt(Q) );
}

void print_npArrays(const QPottsModel& model, const std::vector<double>& g_exact, const std::vector<double>& ln_gAbs_wl) {
    const auto Q = model.getQ();
    const auto L = model.getSize();
    const auto dim = model.getDim();
    // energy range
    const int min_energy = model.getMinEnergy();
    constexpr int max_energy = 0;
    const std::size_t nBins = max_energy - min_energy + 1;

    // name suffix: QLd
    const std::string suffix =
        std::to_string(Q) +
        std::to_string(L) +
        std::to_string(dim);

    // ---- Energies ----
    std::cout << "E_" << suffix << " = np.array([";
    for (std::size_t i = 0; i < nBins; ++i) {
        std::cout << (min_energy + static_cast<int>(i));
        if (i + 1 < nBins) std::cout << ", ";
    }
    std::cout << "])\n\n";

    // ---- Exact g(E) ----
    std::cout << "lng_exact" << suffix << " = np.array([";
    for (std::size_t i = 0; i < nBins; ++i) {
        std::cout << g_exact[i];
        if (i + 1 < nBins) std::cout << ", ";
    }
    std::cout << "])\n\n";

    // ---- WL g(E) ----
    std::cout << "lng_wl" << suffix << " = np.array([";
    for (std::size_t i = 0; i < nBins; ++i) {
        if (ln_gAbs_wl[i] > 0.0) { std::cout << std::exp(ln_gAbs_wl[i]); }
        else { std::cout << "0"; }
        if (i + 1 < nBins) std::cout << ", ";
    }
    std::cout << "])\n";
}

void writeWLOutputCSV(const std::string &filename, const std::vector<int> &energies, const std::vector<double> &log_densities) {
    std::ofstream file(filename);
    file << "E,ln_g\n";  // header for Python

    for (size_t i = 0; i < energies.size(); ++i) {
        file << energies[i] << "," << log_densities[i] << "\n";
    }

    file.close();
}


void runWLandPrint(QPottsModel &model, int Q, int L, int dimension, int N) {
    {
        std::vector<std::vector<double>> all_g;
        all_g.reserve(N);

        // ---- run N independent WL simulations ----
        for (int k = 0; k < N; ++k) {
            WangLandauSimulator wl(model);

            const auto t0 = std::chrono::high_resolution_clock::now();
            wl.runSimulation_parallel();
            const auto duration =
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - t0
                ).count();

            std::cout << "Running WL in parallel with (Q, L, d) = ("
                      << Q << ", " << L << ", " << dimension
                      << ") took " << duration / 1000.0 << " seconds\n";

            all_g.emplace_back(wl.compute_lngAbs());
        }

        // ---- energy range ----
        const int min_energy = model.getMinEnergy();
        const int max_energy = QPottsModel::getMaxEnergy();
        const std::size_t nBins = max_energy - min_energy + 1;

        // ---- NumPy name suffix ----
        const std::string suffix =
            std::to_string(Q) +
            std::to_string(L) +
            std::to_string(dimension);

        // ---- print as 2D NumPy array ----
        std::cout << "lng_wl" << suffix << " = np.array([\n";
        for (std::size_t r = 0; r < all_g.size(); ++r) {
            std::cout << "  [";
            for (std::size_t i = 0; i < nBins; ++i) {
                std::cout << all_g[r][i];
                if (i + 1 < nBins) std::cout << ", ";
            }
            std::cout << "]";
            if (r + 1 < all_g.size()) std::cout << ",";
            std::cout << "\n";
        }
        std::cout << "])\n";
    }
}

void runWLonceWriteCSV(QPottsModel &model, int Q, int L, int dimension, const std::string &output_file) {
    {
        // Construct the WL-simulator
        WangLandauSimulator wl(model);

        // Run the WL algorithm
        const auto t0 = std::chrono::high_resolution_clock::now();
        wl.runSimulation_parallel();
        const auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - t0
            ).count();

        std::cout << "Running WL in parallel with (Q, L, d) = ("
                  << Q << ", " << L << ", " << dimension << ") took "
                  << duration / 1000.0 << " seconds\n\n";

        // Get ln(g_abs) from the WL output
        const auto ln_gAbs_wl = wl.compute_lngAbs();

        // Energies vector (inclusive range is usually what you want)
        std::vector<int> energies;
        energies.reserve(static_cast<std::size_t>(QPottsModel::getMaxEnergy() - wl.getMinEnergy() + 1));

        for (int e = wl.getMinEnergy(); e <= QPottsModel::getMaxEnergy(); ++e) {
            energies.emplace_back(e);
        }

        writeWLOutputCSV(output_file, energies, ln_gAbs_wl);
    }
}
