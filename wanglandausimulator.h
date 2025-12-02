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
        log_densities.resize(max_energy - min_energy + 1, 0.0);
    }

    // GET FUNCTIONS
    std::vector<size_t> getHistogram() { return histogram; }
    std::vector<double> getLogDensities() { return log_densities; }

    // HELP FUNCTIONS
    void resetHistogram() { std::ranges::fill(histogram, 0); }

    [[nodiscard]] bool isHistogramFlat(const double criterion) const {
        int visited_count = 0; // number of bins we have visited
        int total_hist = 0; // sum of all counts

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

    void plotHistogram(const std::string& filename) const {
        std::ofstream file(filename);
        file << "energy,count\n";
        for (int e=min_energy; e<=max_energy; ++e) {
            int idx = e - min_energy;
            file << e << "," << histogram[idx] << "\n";
        }
        file.close();
    }

    void plotLogDensities(const std::string& filename) const {
        std::ofstream file(filename);
        file << "energy,density\n";
        for (int e=min_energy; e<=max_energy; ++e) {
            int idx = e - min_energy;
            file << e << "," << log_densities[idx] << "\n";
        }
        file.close();
    }

    [[nodiscard]] int energy_to_index(const double energy) const {
        const int idx = static_cast<int>(std::round(energy) - min_energy);

        if (idx < 0 || idx >= histogram.size()) {
            throw std::out_of_range("Energy index out of bounds!");
        }

        return idx;
    }

    // ALGORITHM IMPLEMENTATION
    void runSimulation(const double tolerance, const double flatness_criterion) {
        double current_energy = lattice.computeEnergy();
        int counter=0;
        int check_interval = 10'000;

        int current_energy_bin = energy_to_index(current_energy);
        if (current_energy_bin < 0 || current_energy_bin >= log_densities.size()) {
            std::cerr << "Start energy out of bounds!\n";
            return;
        }

        std::cout << "Starting Wang-Landau. Range: [" << min_energy << ", " << max_energy << "]\n";

        while (ln_f > tolerance) {
            size_t random_site_idx = spin_index_dist(gen);
            int current_spin = lattice.getSpin(random_site_idx);
            int random_spin;

            do {
                random_spin = spin_value_dist(gen);
            } while (random_spin == current_spin);

            const double dE = lattice.computeEnergyDifference(random_site_idx, random_spin);
            const double new_energy = current_energy + dE;
            int new_energy_bin = energy_to_index(new_energy);

            if (new_energy_bin < 0 || new_energy_bin >= log_densities.size()) {
                std::cout << "Energy out of bounds...";
                continue;
            }

            // ln(g_old) - ln(g_new)
            double const lnP = log_densities[current_energy_bin] - log_densities[new_energy_bin];

            // min( 1 , e^(lnP) )
            double const P_acc = std::min(1.0, std::exp(lnP));

            if (lnP >= 0.0 || P_acc > prob_dist(gen)) {
                // ACCEPT
                lattice.setSpin(random_site_idx, random_spin);
                current_energy=new_energy;
                current_energy_bin=new_energy_bin;
            }

            // UPDATE DENSITY EN HISTOGRAM
            histogram[current_energy_bin] += 1;
            log_densities[current_energy_bin] += ln_f;


            // CHECK FLATNESS
            if (counter % check_interval == 0 && counter > 0 && isHistogramFlat(flatness_criterion)) {
                std::cout << "Histogram flat. Reducing f. (ln_f = " << ln_f * 0.5 << ")\n";
                resetHistogram();
                ln_f *= 0.5;
            }

            ++counter;
        }
    }
};


#endif //SI2530PROJECT_WANGLANDAUSIMULATOR_H