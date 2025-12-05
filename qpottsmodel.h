//
// Created by quinten on 11/12/25.
//

#ifndef SI2530PROJECT_LATTICE_H
#define SI2530PROJECT_LATTICE_H

#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>


class QPottsModel {
private:
    size_t dim;
    size_t L;
    size_t Q;
    double J;
    size_t total_sites;
    std::vector<int> spins;

    [[nodiscard]] size_t getIndex(const std::vector<size_t>& coordinates) const;
    [[nodiscard]] std::vector<size_t> getCoordinates(size_t index) const;

public:
    QPottsModel(const size_t dimension, const size_t size) : QPottsModel(dimension, size, 2, 1.0) {}
    QPottsModel(const size_t dimension, const size_t size, const size_t q, const int j) : dim(dimension), L(size), Q(q), J(j) {
        if (dimension == 0) { throw std::invalid_argument("Dimension should be at least 1!"); }
        if (size == 0) { throw std::invalid_argument("Size should be at least 1!"); }
        if (Q<2) { throw std::invalid_argument("Q should be larger than 2!"); }
        if (J<=0) { throw std::invalid_argument("J must be > 0!"); }

        size_t total = 1;
        for (size_t d=0; d<dim; ++d) {
            total *= L;
        }

        total_sites=total;
        spins.resize(total_sites, 0);
    }

    [[nodiscard]] size_t getDim() const { return dim; }
    [[nodiscard]] size_t getSize() const { return L; }
    [[nodiscard]] std::vector<int> getSpins() const { return spins; }
    [[nodiscard]] size_t getTotalSites() const { return total_sites; }
    [[nodiscard]] size_t getQ() const { return Q; }
    [[nodiscard]] int getSpin(const size_t index) const { return spins[index]; }
    [[nodiscard]] int getSpin(const std::vector<size_t>& coordinates) const { return spins[getIndex(coordinates)]; }
    [[nodiscard]] int getMinEnergy() const { return static_cast<int>(-J * dim * getTotalSites()); }
    static int getMaxEnergy() { return 0; }

    void setSpin(size_t index, int new_spin);
    void setSpin(const std::vector<size_t>& coordinates, int new_spin);

    [[nodiscard]] double computeEnergy() const;
    [[nodiscard]] double computeEnergyDifference(size_t idx, int new_spin) const;
    [[nodiscard]] double computeEnergyDifference(const std::vector<size_t>& coordinates, int new_spin) const;
};


#endif //SI2530PROJECT_LATTICE_H