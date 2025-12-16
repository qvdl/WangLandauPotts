#include "wanglandausimulator.h"

#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <omp.h>


WangLandauSimulator::WangLandauSimulator(QPottsModel &lattice)
    : lattice(lattice),
      min_energy(lattice.getMinEnergy()),
      max_energy(QPottsModel::getMaxEnergy()),
      ln_f(1.0),
      check_interval(10'000)
{
    const size_t n_bins = max_energy - min_energy + 1;
    histogram.resize(n_bins, 0);
    log_densities.resize(n_bins, 0.0);
}



// HELP FUNCTIONS
void WangLandauSimulator::resetHistogram() { std::ranges::fill(histogram, 0); }

bool WangLandauSimulator::isHistogramFlat(const double criterion) const {
    int visited_count = 0; // number of bins we have visited
    int total_hist = 0;    // sum of all counts

    for (auto& count : histogram) {
        if (count > 0) {
            visited_count++;
            total_hist += static_cast<int>(count);
        }
    }

    // if we have not visited enough energies yet, histogram is not flat
    if (static_cast<double>(visited_count) / static_cast<double>(histogram.size()) < 0.25) { return false; }

    // compute average
    const double avg_hist = static_cast<double>(total_hist) / static_cast<double>(visited_count);

    // check flatness for visited energies
    for (size_t E=0; E<histogram.size(); ++E) {
        if (histogram[E] > 0 && histogram[E] < criterion * avg_hist) { return false; }
    }

    return true;
}

void WangLandauSimulator::plotHistogram(const std::string& filename) const {
    std::ofstream file(filename);
    file << "energy,count\n";
    for (int e = min_energy; e <= max_energy; ++e) {
        const int idx = e - min_energy;
        file << e << "," << histogram[static_cast<size_t>(idx)] << "\n";
    }
    file.close();
}

void WangLandauSimulator::plotLogDensities(const std::string& filename) const {
    std::ofstream file(filename);
    file << "energy,density\n";
    for (int e = min_energy; e <= max_energy; ++e) {
        const int idx = e - min_energy;
        file << e << "," << log_densities[static_cast<size_t>(idx)] << "\n";
    }
    file.close();
}

int WangLandauSimulator::energy_to_index(const double energy) const {
    if (energy < min_energy-1e-12|| energy > max_energy+1e-12) {
        throw std::invalid_argument("Energy is out of bound!");
    }
    return static_cast<int>(std::round(energy) - min_energy);
}

// ALGORITHM IMPLEMENTATION
void WangLandauSimulator::runSimulation(const double tolerance, const double flatness_criterion) {
    double current_energy = lattice.computeEnergy();
    int counter = 0;
    int current_energy_bin = energy_to_index(current_energy);

    // std::cout << "Starting Wang-Landau. Range: [" << min_energy << ", " << max_energy << "]\n";

    while (ln_f > tolerance) {
        const size_t random_site_idx = RNG::uniformIndex(0, lattice.getTotalSites() - 1);
        const int current_spin = lattice.getSpin(random_site_idx);
        int random_spin;

        do {
            random_spin = RNG::uniformSpin(static_cast<int>(lattice.getQ()));
        } while (random_spin == current_spin);

        const double dE = lattice.computeEnergyDifference(random_site_idx, random_spin);
        const double new_energy = current_energy + dE;
        const int new_energy_bin = energy_to_index(new_energy);

        // ln(g_old) - ln(g_new)
        double const lnP = log_densities[current_energy_bin] - log_densities[new_energy_bin];

        // min{ 1 , e^(lnP) }
        double const P_acc = std::min(1.0, std::exp(lnP));

        if (lnP >= 0.0 || P_acc > RNG::uniformProb()) {
            // ACCEPT
            lattice.setSpin(random_site_idx, random_spin);
            current_energy = new_energy;
            current_energy_bin = new_energy_bin;
        }

        // UPDATE DENSITY EN HISTOGRAM
        histogram[current_energy_bin] += 1;
        log_densities[current_energy_bin] += ln_f;

        // CHECK FLATNESS
        if (counter % check_interval == 0 && counter > 0 && isHistogramFlat(flatness_criterion)) {
            // std::cout << "Histogram flat. Reducing f. (ln_f = " << ln_f * 0.5 << ")\n";
            resetHistogram();
            ln_f *= 0.5;
        }

        ++counter;
    }
}


void WangLandauSimulator::runSimulation_parallel(double tolerance, double flatness_criterion)
{
#pragma omp parallel
{
    // Thread-local persistent walker
    QPottsModel local_lattice = lattice;             // copy once
    double current_energy = local_lattice.computeEnergy();
    int current_energy_bin = energy_to_index(current_energy);

    const std::size_t n_sites = local_lattice.getTotalSites();
    const int q = static_cast<int>(local_lattice.getQ());

    bool done = false;
    while (!done) {
        // --- chunk of MC attempts per thread ---
        for (int step = 0; step < check_interval; ++step) {

            const std::size_t i = RNG::uniformIndex(0, n_sites - 1);
            const int s_old = local_lattice.getSpin(i);
            int s_new; do { s_new = RNG::uniformSpin(q); } while (s_new == s_old);

            const double dE = local_lattice.computeEnergyDifference(i, s_new);
            const double newE = current_energy + dE;
            const int newBin = energy_to_index(newE);

            // hard reject if out of range
            if (newBin >= 0 && newBin < (int)log_densities.size()) {
                const double lnP = log_densities[current_energy_bin] - log_densities[newBin];
                if (lnP >= 0.0 || std::exp(lnP) > RNG::uniformProb()) {
                    local_lattice.setSpin(i, s_new);
                    current_energy = newE;
                    current_energy_bin = newBin;
                }
            }
            // WL updates at the visited (current) bin
#pragma omp atomic update
            histogram[current_energy_bin] += 1;
#pragma omp atomic update
            log_densities[current_energy_bin] += ln_f;
        }

        // --- sync + flatness/√f ---
#pragma omp barrier
#pragma omp single
        {
            if (isHistogramFlat(flatness_criterion)) {
                ln_f *= 0.5;         // f -> sqrt(f)
                resetHistogram();
            }
            if (ln_f <= tolerance) {
                // tell everyone to exit
                // (captured by all threads after the barrier)
            }
        }
#pragma omp barrier
        if (ln_f <= tolerance) done = true;
    } // while
} // parallel
}


std::vector<double> WangLandauSimulator::compute_lngAbs() const {
    const double ln_g_E0 = log_densities[0];
    const double ln_Q = std::log(lattice.getQ());

    std::vector<double> ln_g_abs(log_densities.size(), 0.0);

    for (std::size_t i=0; i<log_densities.size(); ++i) {
        if (log_densities[i] > 0) { ln_g_abs[i] = log_densities[i] - ln_g_E0 + ln_Q; }
    }

    return ln_g_abs;
}
