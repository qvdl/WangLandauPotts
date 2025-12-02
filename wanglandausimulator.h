//
// Created by quinten on 11/12/25.
//

#ifndef SI2530PROJECT_WANGLANDAUSIMULATOR_H
#define SI2530PROJECT_WANGLANDAUSIMULATOR_H


#include <algorithm>
#include <random>
#include <numeric>
#include <fstream>


template<typename ModelType>
class WangLandauSimulator {
private:
    ModelType& lattice;
    int min_energy;
    int max_energy;
    double ln_f;
    std::vector<size_t> histogram;
    std::vector<double> log_densities;
    std::mt19937 gen{std::random_device{}()};
    std::uniform_int_distribution<size_t> spin_index_dist;
    std::uniform_int_distribution<int> spin_value_dist;
    std::uniform_real_distribution<double> prob_dist{0.0, 1.0};

public:
    explicit WangLandauSimulator(ModelType &lattice) : lattice(lattice), min_energy(lattice.getMinEnergy()),
                                                       max_energy(lattice.getMaxEnergy()), ln_f(1.0),
                                                       spin_index_dist(0, lattice.getTotalSites() - 1),
                                                       spin_value_dist(0, lattice.getQ() - 1) {
        histogram.resize(max_energy - min_energy + 1, 0);
        log_densities.resize(max_energy - min_energy + 1, 0);
    }

    // GET FUNCTIONS
    std::vector<size_t> getHistogram() { return histogram; }
    std::vector<double> getLogDensities() { return log_densities; }

    // HELP FUNCTIONS
    void resetHistogram() { std::ranges::fill(histogram, 0); }

    [[nodiscard]] bool isHistogramFlat(const double criterion) const {
        std::vector<size_t> reachable_counts;
        for (const auto& count : histogram) {
            if (count > 0) {
                reachable_counts.emplace_back(count);
            }
        }

        if (reachable_counts.empty()) { return false; }

        const double mean = std::accumulate(reachable_counts.begin(), reachable_counts.end(), 0.0) / static_cast<double>(reachable_counts.size());
        const double min = *std::ranges::min_element(reachable_counts.begin(), reachable_counts.end());
        return (min > mean*criterion);
    }

    void plotHistogram(const std::string& filename) const {
        std::ofstream file(filename);
        file << "energy,count\n";
        for (int e=min_energy; e<=max_energy; ++e) {
            int idx = e - min_energy;
            file << e << "," << histogram[idx] << "\n";
        }
        file.close();
    }

    void plotDensities(const std::string& filename) const {
        std::ofstream file(filename);
        file << "energy,density\n";
        for (int e=min_energy; e<=max_energy; ++e) {
            int idx = e - min_energy;
            file << e << "," << std::exp(log_densities[idx]) << "\n";
        }
        file.close();
    }

    [[nodiscard]] int energy_to_index(double energy) const {
        const int idx = static_cast<int>(std::round(energy) - min_energy);

        if (idx < 0 || idx >= histogram.size()) {
            // DEBUG: print wat er misgaat
            std::cerr << "ENERGY BUG: energy=" << energy
                      << " min=" << min_energy << " max=" << max_energy
                      << " idx=" << idx << " hist_size=" << histogram.size() << std::endl;
            throw std::out_of_range("Energy index out of bounds");
        }
        return idx;
    }

    // ALGORITHM IMPLEMENTATION
    void runSimulation(const size_t mc_sweeps, const double tolerance, double flatness_criterion) {
        double current_energy = lattice.computeEnergy();
        int counter=0;

        try {
            if (energy_to_index(current_energy) < 0) {
                throw std::out_of_range("Initial energy is outside the defined range.");
            }
        } catch (const std::out_of_range& e) {
            std::cerr << "Simulation aborted: " << e.what() << std::endl;
            return;
        }

        while (ln_f > tolerance && counter < mc_sweeps) {
            size_t random_idx = spin_index_dist(gen);
            int random_spin = spin_value_dist(gen);

            const double dE = lattice.computeEnergyDifference(random_idx, random_spin);
            const double new_energy = current_energy + dE;

            // ACCEPT or REJECT
            int current_idx;
            int new_idx;

            try {
                current_idx = energy_to_index(current_energy);
                new_idx = energy_to_index(new_energy);
            } catch (const std::out_of_range& e) {
                // if the new energy is out of bounds, do nothing
                continue;
            }

            // ln(g_old) - ln(g_new)
            double const lnP = log_densities[current_idx] - log_densities[new_idx];

            // min( 1 , e^(lnP) )
            double const P_acc = std::min(1.0, std::exp(lnP));

            if (P_acc > prob_dist(gen)) {
                // ACCEPT
                lattice.setSpin(random_idx, random_spin);
                current_energy=new_energy;
                current_idx=new_idx;
            }

            // UPDATE DENSITY EN HISTOGRAM
            histogram[current_idx] += 1;
            log_densities[current_idx] += ln_f;


            // CHECK FLATNESS
            if (counter % 10'000 == 0 && counter > 0 && isHistogramFlat(flatness_criterion)) {
                resetHistogram();
                ln_f *= 0.5;
            }

            ++counter;
        }
    }
};


#endif //SI2530PROJECT_WANGLANDAUSIMULATOR_H