//
// Created by quinten on 12/5/25.
//

#include "ThermoState.h"
#include <iostream>

ThermoState::ThermoState(const int min_energy,
                const std::vector<double>& log_g, // ln g(E) for E = Emin + i
                const double T,
                const double kB)
        : Emin_(min_energy), T_(T)
{
    const std::size_t n = log_g.size();
    P_.assign(n, 0.0);

    // Precompute log-weights ln w_i = ln g(E_i) - beta * E_i
    const double beta = 1.0 / (kB * T_);
    std::vector<double> lnw(n);
    // Also keep E_i = Emin + i (as double) for expectations
    // (no need to store E permanently)
    for (std::size_t i = 0; i < n; ++i) {
        const auto Ei = static_cast<double>(Emin_ + static_cast<int>(i));
        lnw[i] = log_g[i] - beta * Ei;
    }

    // log-sum-exp for normalization: logZ = m + log( sum exp(lnw - m) )
    const double m = *std::max_element(lnw.begin(), lnw.end());
    double sum_exp = 0.0;
    for (std::size_t i = 0; i < n; ++i) sum_exp += std::exp(lnw[i] - m);
    const double logZ = m + std::log(sum_exp);

    // Build P(E|T), and compute <E>, <E^2>
    double E1 = 0.0, E2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto Ei = static_cast<double>(Emin_ + static_cast<int>(i));
        const double Pi = std::exp(lnw[i] - logZ);
        P_[i] = Pi;
        E1 += Ei * Pi;
        E2 += Ei * Ei * Pi;
    }

    // Store thermodynamic quantities
    E_  = E1;
    Cv_ = (E2 - E1 * E1) / (kB * T_ * T_);
    F_  = -kB * T_ * logZ;
    S_  = (E1 - F_) / T_;            // also equal to kB*(logZ + beta*E1)
    Z_  = std::exp(logZ);
}
