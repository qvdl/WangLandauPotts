//
// Created by quinten on 12/5/25.
//

#ifndef SI2530PROJECT_THERMOSTATE_H
#define SI2530PROJECT_THERMOSTATE_H

#include <vector>
#include <cmath>
#include <algorithm> // std::max_element
#include <cstddef>   // size_t

class ThermoState {
private:
    int Emin_{0};
    double T_{0.0};

    double E_{0.0};
    double Cv_{0.0};
    double F_{0.0};
    double S_{0.0};
    double Z_{0.0};
    std::vector<double> P_;

public:
    ThermoState(int min_energy, const std::vector<double>& log_g, double T, double kB = 1.0);

    [[nodiscard]] double T()  const { return T_; }
    [[nodiscard]] int    Emin() const { return Emin_; }

    [[nodiscard]] double E()  const { return E_; }
    [[nodiscard]] double Cv() const { return Cv_; }
    [[nodiscard]] double F()  const { return F_; }
    [[nodiscard]] double S()  const { return S_; }
    [[nodiscard]] double Z()  const { return Z_; }

    [[nodiscard]] const std::vector<double>& P() const { return P_; } // P(E|T) over E = Emin + i
};

#endif // SI2530PROJECT_THERMOSTATE_H
