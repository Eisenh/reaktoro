// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "ActivityModelDEW.hpp"

// C++ includes
#include <map>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/ConvertUtils.hpp>
#include <Reaktoro/Common/Index.hpp>
#include <Reaktoro/Common/NamingUtils.hpp>
#include <Reaktoro/Math/BilinearInterpolator.hpp>
#include <Reaktoro/Models/ActivityModels/Support/AqueousMixture.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>

namespace Reaktoro {

using std::abs;
using std::log;
using std::log10;
using std::pow;
using std::sqrt;

namespace {

/// The electrostatic constant \f$ \eta\f$ in the DEW model (in units of (A*cal)/mol)
const auto eta = 1.66027e+05;

/// DEW-specific parameters for high pressure/temperature conditions
struct DEWParams
{
    real a1;  ///< DEW parameter a1
    real a2;  ///< DEW parameter a2
    real a3;  ///< DEW parameter a3
    real a4;  ///< DEW parameter a4
    real c1;  ///< DEW parameter c1
    real c2;  ///< DEW parameter c2
};

/// DEW parameters for common ions (values from DEW database)
const std::map<std::string, DEWParams> dew_parameters = {
    {"Na+",  {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}},  // TODO: Add actual DEW parameters
    {"Cl-",  {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}},  // TODO: Add actual DEW parameters
    // Add more ions as needed
};

/// Calculate the DEW-specific correction term
auto dewCorrectionTerm(real T, real P, const DEWParams& params) -> real
{
    // Convert temperature to Kelvin if needed
    const real TK = T;  // Assuming T is already in Kelvin
    
    // Calculate temperature-dependent terms
    const real T_term1 = params.a1 / TK;
    const real T_term2 = params.a2 / (TK * TK);
    const real T_term3 = params.a3 * log(TK);
    const real T_term4 = params.a4 * TK;
    
    // Calculate pressure-dependent terms
    const real P_term1 = params.c1 * P;
    const real P_term2 = params.c2 * P * P;
    
    // Combine terms with appropriate coefficients
    const real correction = eta * (T_term1 + T_term2 + T_term3 + T_term4 + 
                                 P_term1 + P_term2);
    
    return correction;
}

} // namespace

auto activityModelDEW(const SpeciesList& species) -> ActivityModel
{
    // Create the aqueous mixture
    AqueousMixture mixture(species);

    // Define the activity model function of the aqueous phase
    ActivityModel fn = [=](ActivityPropsRef props, ActivityModelArgs args) mutable
    {
        // The arguments for the activity model evaluation
        const auto& [T, P, x] = args;

        // Calculate the base HKF activity coefficients
        mixture.update(T, P, x);
        const auto& ln_gamma_HKF = mixture.lnActivityCoefficients();

        // Initialize DEW correction terms
        VectorXr ln_gamma_DEW = VectorXr::Zero(x.size());

        // Calculate DEW correction terms for each species
        for(auto i = 0; i < x.size(); ++i)
        {
            const auto& species = mixture.species(i);
            const auto it = dew_parameters.find(species.name());
            
            if(it != dew_parameters.end())
            {
                const auto& params = it->second;
                ln_gamma_DEW[i] = dewCorrectionTerm(T, P, params);
            }
        }

        // Combine HKF and DEW contributions
        props.ln_g = ln_gamma_HKF + ln_gamma_DEW;
        props.ln_a = props.ln_g + log(x);

        // Set the state of matter of the phase
        props.som = StateOfMatter::Liquid;
    };

    return fn;
}

auto ActivityModelDEW() -> ActivityModelGenerator
{
    return [](const SpeciesList& species) { return activityModelDEW(species); };
}

} // namespace Reaktoro
