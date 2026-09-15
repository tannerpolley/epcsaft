#include "properties/helmholtz_model.hpp"
#include "surface/parameter_builder.hpp"

#include <cctype>
#include <cmath>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>

#include "../../build/subprojects/packagefiles/asl/funcadd.h"

namespace {

using epcsaft::DerivativeOrder;
using epcsaft::NeutralProperty;
using epcsaft::NeutralPropertyDerivatives;

struct PreparedModel final {
    std::mutex mutex;
    std::tuple<std::string, std::string, std::string, std::string> inputs;
    std::unique_ptr<epcsaft::HelmholtzModel> model;
};

PreparedModel& prepared_model() {
    static PreparedModel result;
    return result;
}

thread_local std::string callback_error;

double invalid_value() {
    return std::numeric_limits<double>::quiet_NaN();
}

void invalidate(arglist* arguments) {
    if (arguments == nullptr) return;
    if (arguments->derivs != nullptr) {
        for (int index = 0; index < arguments->nr; ++index) {
            arguments->derivs[index] = invalid_value();
        }
    }
    if (arguments->hes != nullptr) {
        const int count = arguments->nr * (arguments->nr + 1) / 2;
        for (int index = 0; index < count; ++index) {
            arguments->hes[index] = invalid_value();
        }
    }
}

[[noreturn]] void reject(const std::string& message) {
    throw std::invalid_argument(message);
}

double parse_decimal(const char* value, const char* name) {
    if (value == nullptr || *value == '\0') reject(std::string("missing ") + name);
    const std::string text(value);
    for (const char character : text) {
        if (std::isspace(static_cast<unsigned char>(character)) != 0) {
            reject(std::string("invalid ") + name);
        }
    }
    std::size_t end = 0;
    double result = 0.0;
    try {
        result = std::stod(text, &end);
    } catch (const std::exception&) {
        reject(std::string("invalid ") + name);
    }
    if (end != text.size() || !std::isfinite(result)) {
        reject(std::string("invalid ") + name);
    }
    return result;
}

void validate_arguments(arglist* arguments) {
    if (arguments == nullptr || arguments->n != 6 || arguments->nr != 2
        || arguments->at == nullptr || arguments->ra == nullptr
        || arguments->sa == nullptr) {
        reject("epcsaft neutral function requires six arguments (two real, four string)");
    }
    for (int position = 0; position < arguments->n; ++position) {
        const int expected = position < 2 ? position : -(position - 1);
        if (arguments->at[position] != expected) reject("invalid neutral argument map");
    }
}

double real_argument(arglist* arguments, int position) {
    const int index = arguments->at[position];
    if (index < 0 || index >= arguments->nr) reject("invalid real argument map");
    return arguments->ra[index];
}

const char* string_argument(arglist* arguments, int position) {
    const int index = arguments->at[position];
    if (index >= 0) reject("invalid string argument map");
    const char* value = arguments->sa[-index - 1];
    if (value == nullptr) reject("missing string argument");
    return value;
}

epcsaft::HelmholtzModel& model_for(
    const std::tuple<std::string, std::string, std::string, std::string>& inputs
) {
    auto& prepared = prepared_model();
    if (prepared.model != nullptr) {
        if (prepared.inputs != inputs) reject("immutable neutral parameter tuple changed");
        return *prepared.model;
    }

    const auto [m, sigma, epsilon, molar_mass] = inputs;
    epcsaft::ParameterDocument document;
    document.component_ids = {"component-0"};
    document.electrolyte_choice = "none";
    document.ionic_dispersion_choice = "";
    document.permittivity_choice = "none";
    auto add = [&](const char* field, double magnitude, const char* unit) {
        epcsaft::ParameterEntry entry;
        entry.identity = std::string("component/0/") + field;
        entry.field = field;
        entry.index = {0};
        entry.magnitude = magnitude;
        entry.unit = unit;
        document.entries.push_back(std::move(entry));
    };
    add("segment_count", parse_decimal(m.c_str(), "segment_count"), "dimensionless");
    add("sigma_constant", parse_decimal(sigma.c_str(), "sigma_constant"), "angstrom");
    add(
        "dispersion_energy_over_k",
        parse_decimal(epsilon.c_str(), "dispersion_energy_over_k"),
        "kelvin"
    );
    add("molar_mass", parse_decimal(molar_mass.c_str(), "molar_mass"), "kilogram / mole");
    add("charge_number", 0.0, "elementary-charge");
    prepared.model = std::make_unique<epcsaft::HelmholtzModel>(
        epcsaft::build_parameters(document)
    );
    prepared.inputs = inputs;
    return *prepared.model;
}

double evaluate(
    arglist* arguments,
    NeutralProperty property
) {
    validate_arguments(arguments);
    const auto inputs = std::make_tuple(
        std::string(string_argument(arguments, 2)),
        std::string(string_argument(arguments, 3)),
        std::string(string_argument(arguments, 4)),
        std::string(string_argument(arguments, 5))
    );
    const double temperature = real_argument(arguments, 0);
    const double density = real_argument(arguments, 1);
    auto& prepared = prepared_model();
    std::lock_guard<std::mutex> lock(prepared.mutex);
    const auto order = arguments->hes != nullptr
        ? DerivativeOrder::Hessian
        : (arguments->derivs != nullptr ? DerivativeOrder::Gradient : DerivativeOrder::Value);
    const NeutralPropertyDerivatives result = model_for(inputs).neutral_property(
        temperature, density, property, order
    );
    if (arguments->derivs != nullptr) {
        for (int position = 0; position < arguments->n; ++position) {
            const int index = arguments->at[position];
            if (index >= 0) {
                if (index >= arguments->nr) reject("invalid derivative map");
                arguments->derivs[index] = result.gradient.at(position);
            }
        }
    }
    if (arguments->hes != nullptr) {
        arguments->hes[0] = result.hessian.at(0);
        arguments->hes[1] = result.hessian.at(1);
        arguments->hes[2] = result.hessian.at(3);
    }
    return result.value;
}

double callback(arglist* arguments, NeutralProperty property) {
    if (arguments != nullptr) arguments->Errmsg = nullptr;
    invalidate(arguments);
    try {
        return evaluate(arguments, property);
    } catch (const std::exception& error) {
        callback_error = arguments != nullptr && arguments->hes != nullptr
            ? std::string("\"") + error.what()
            : (arguments != nullptr && arguments->derivs != nullptr
                ? std::string("'") + error.what()
                : error.what());
        if (arguments != nullptr) arguments->Errmsg = callback_error.data();
        return invalid_value();
    } catch (...) {
        callback_error = "unknown epcsaft neutral callback error";
        if (arguments != nullptr) arguments->Errmsg = callback_error.data();
        return invalid_value();
    }
}

double residual_helmholtz(arglist* arguments) {
    return callback(arguments, NeutralProperty::ResidualHelmholtz);
}

double pressure(arglist* arguments) {
    return callback(arguments, NeutralProperty::Pressure);
}

double recordings(arglist* arguments) {
    if (arguments == nullptr || arguments->n != 0 || arguments->nr != 0) {
        callback_error = "epcsaft neutral recordings requires zero arguments";
        if (arguments != nullptr) arguments->Errmsg = callback_error.data();
        return invalid_value();
    }
    auto& prepared = prepared_model();
    std::lock_guard<std::mutex> lock(prepared.mutex);
    return prepared.model == nullptr
        ? 0.0
        : static_cast<double>(prepared.model->counters().recordings);
}

}  // namespace

extern "C" void funcadd_ASL(AmplExports* exports) {
    if (exports == nullptr || exports->Addfunc == nullptr) return;
    exports->Addfunc(
        "epcsaft_neutral_ar", residual_helmholtz,
        FUNCADD_REAL_VALUED | FUNCADD_STRING_ARGS, 6, nullptr, exports
    );
    exports->Addfunc(
        "epcsaft_neutral_pressure", pressure,
        FUNCADD_REAL_VALUED | FUNCADD_STRING_ARGS, 6, nullptr, exports
    );
    exports->Addfunc(
        "epcsaft_neutral_recordings", recordings,
        FUNCADD_REAL_VALUED, 0, nullptr, exports
    );
}
