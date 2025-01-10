#include <Reaktoro/pybind11.hxx>
#include <Reaktoro/Models/ActivityModels/ActivityModelDEW.hpp>

using namespace Reaktoro;

void exportActivityModelDEW(py::module& m)
{
    using namespace Reaktoro;

    m.def("ActivityModelDEW", &ActivityModelDEW,
        "Return an activity model for aqueous species based on the DEW (Deep Earth Water) model.\n"
        "The DEW model extends the HKF model with additional terms for high pressure/temperature conditions.\n"
        "\n"
        "Returns\n"
        "-------\n"
        "ActivityModel\n"
        "    The activity model function for aqueous species based on the DEW model.");
}
