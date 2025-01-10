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

// Catch includes
#include <catch2/catch.hpp>

// Reaktoro includes
#include <Reaktoro/Models/ActivityModels/ActivityModelDEW.hpp>
#include <Reaktoro/Water/WaterConstants.hpp>
using namespace Reaktoro;

namespace {

/// Return mole fractions for the species.
inline auto moleFractions(const SpeciesList& species) -> ArrayXr
{
    auto idx = [&](auto formula) { return species.indexWithFormula(formula); };

    ArrayXr n(species.size());
    n = 0.1;
    n[idx("H2O")] = 55.508;
    n[idx("H+" )] = 1e-7;
    n[idx("OH-")] = 1e-7;
    n[idx("Na+")] = 0.3;
    n[idx("Cl-")] = 0.3;

    return n / n.sum();
}

// Check if the activities of the aqueous species are correct assuming activity coefficients are.
inline auto checkActivities(ArrayXrConstRef x, ActivityPropsConstRef props)
{
    const auto iH2O = 0;

    // The concentrations of the species (molalities for solutes, mole fraction for solvent water)
    ArrayXr c = x/(x[iH2O] * waterMolarMass);
    c[iH2O] = x[iH2O];

    for(auto i = 0; i < x.size(); ++i)
    {
        INFO("i = " << i);
        CHECK( exp(props.ln_a[i] - props.ln_g[i]) == Approx(c[i]) );
    }
}

} // anonymous namespace

TEST_CASE("Testing ActivityModelDEW", "[ActivityModelDEW]")
{
    const auto species = SpeciesList("H2O H+ OH- Na+ Cl- Ca++ HCO3- CO3-- CO2 NaCl HCl NaOH");

    // Test at standard conditions
    SECTION("Testing at standard conditions")
    {
        const auto T = 300.0;
        const auto P = 12.3e5;
        const auto x = moleFractions(species);

        ActivityModel fn = ActivityModelDEW()(species);
        ActivityProps props = ActivityProps::create(species.size());
        fn(props, {T, P, x});

        // Verify activity coefficients are calculated
        for(auto i = 0; i < species.size(); ++i)
        {
            INFO("i = " << i);
            CHECK( props.ln_g[i] != 0.0 );
        }

        checkActivities(x, props);
    }

    // Test at high pressure/temperature conditions
    SECTION("Testing at high pressure/temperature conditions")
    {
        const auto T = 500.0;
        const auto P = 5000.0e5; // 5000 bar
        const auto x = moleFractions(species);

        ActivityModel fn = ActivityModelDEW()(species);
        ActivityProps props = ActivityProps::create(species.size());
        fn(props, {T, P, x});

        // Verify activity coefficients are calculated
        for(auto i = 0; i < species.size(); ++i)
        {
            INFO("i = " << i);
            CHECK( props.ln_g[i] != 0.0 );
        }

        checkActivities(x, props);
    }
}
