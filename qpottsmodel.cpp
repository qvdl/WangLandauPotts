//
// Created by quinten on 11/12/25.
//

#include "qpottsmodel.h"


size_t QPottsModel::getIndex(const std::vector<size_t>& coordinates) const {
    if (coordinates.size() != dim) { throw std::invalid_argument("Coordinates should have dimension " + std::to_string(dim)); }

    size_t index = 0;
    size_t multiplier = 1;

    for(int d = 0; d < dim; d++) {
        index += coordinates[d] * multiplier;
        multiplier *= L;
    }
    return index;
}

std::vector<size_t> QPottsModel::getCoordinates(const size_t index) const {
    if (index >= spins.size()) { throw std::invalid_argument("Index out of bounds!"); }

    std::vector<size_t> coords(dim);
    size_t temp = index;

    for(int d = 0; d < dim; d++) {
        coords[d] = temp % L;
        temp /= L;
    }

    return coords;
}

void QPottsModel:: setSpin(const size_t index, const int new_spin) {
    if (index >= total_sites) { throw std::invalid_argument("Index out of range!"); }
    if (new_spin >= Q) { throw std::invalid_argument("Spin cannot be larger than " + std::to_string(Q)); }
    spins[index] = new_spin;
}

void QPottsModel::setSpin(const std::vector<size_t> &coordinates, const int new_spin) {
    if (new_spin >= Q) { throw std::invalid_argument("Spin cannot be larger than " + std::to_string(Q)); }
    spins[getIndex(coordinates)] = new_spin;
}

double QPottsModel::computeEnergy() const {
    double energy=0.0;
    for (int idx=0; idx<getTotalSites(); ++idx) {
        const std::vector<size_t> coordinates = getCoordinates(idx);

        for (int d=0; d<dim; ++d) {
            std::vector<size_t> neighbor_coordinates = coordinates;
            neighbor_coordinates[d] = (coordinates[d]+1) % L;

            const size_t nb_idx = getIndex(neighbor_coordinates);
            if (spins[idx] == spins[nb_idx]) { energy -= J; }
        }
    }

    return energy;
}

double QPottsModel::computeEnergyDifference(const size_t idx, const int new_spin) const {
    const auto coordinates = getCoordinates(idx);
    return computeEnergyDifference(coordinates, new_spin);
}

double QPottsModel::computeEnergyDifference(const std::vector<size_t> &coordinates, const int new_spin) const {
    if (new_spin >= Q) { throw std::invalid_argument("Spin proposal not possible: new spin should be smaller than Q=" + std::to_string(Q)); }

    double delta=0.0;
    const size_t idx = getIndex(coordinates);
    const int old_spin = spins[idx];

    if (old_spin == new_spin) { return delta; }

    for (int d=0; d<dim; ++d) {
        // check energy difference in positive direction
        std::vector<size_t> pos_neighbor_coords = coordinates;
        pos_neighbor_coords[d] = (coordinates[d] + 1) % L;

        const size_t pos_neighbor_idx = getIndex(pos_neighbor_coords);
        const int pos_neighbor_spin = spins[pos_neighbor_idx];

        if (old_spin == pos_neighbor_spin) { delta += J; }
        if (new_spin == pos_neighbor_spin) { delta -= J; }

        // check energy difference in negative direction
        std::vector<size_t> neg_neighbor_coords = coordinates;
        neg_neighbor_coords[d] = (coordinates[d] + L - 1) % L;

        const size_t neg_neighbor_idx = getIndex(neg_neighbor_coords);
        const int neg_neighbor_spin = spins[neg_neighbor_idx];

        if (old_spin == neg_neighbor_spin) { delta += J; }
        if (new_spin == neg_neighbor_spin) { delta -= J; }
    }

    return delta;
}