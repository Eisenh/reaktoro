// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Thermodynamics/EOS/GaseousStateEquation.hpp>

namespace Reaktoro {

/// Return an equation of state for a gaseous phase based on the ideal model.
/// @param mixture The gaseous mixture
/// @return The equation of state function for the gaseous phase
/// @see GaseousMixture, GaseousStateEquation
auto gaseousStateEquationIdeal(const GaseousMixture& mixture) -> GaseousStateEquation;

} // namespace Reaktoro
