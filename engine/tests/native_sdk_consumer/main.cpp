#include <epcsaft/native_model_v1.h>
#include <epcsaft/native_sdk_v1.h>

#include <cstdint>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
    if (argc != 2) {
        return 2;
    }
    char error[EPCSAFT_NATIVE_SDK_V1_ERROR_SIZE] = {};
    epcsaft_native_model_handle_v1* model =
        epcsaft_native_model_load_v1(argv[1], error, sizeof(error));
    if (model == nullptr) {
        std::cerr << error << '\n';
        return 3;
    }
    const epcsaft_native_sdk_v1* sdk = epcsaft_native_model_sdk_v1(model);
    const bool valid = sdk != nullptr
        && sdk->abi_version == EPCSAFT_NATIVE_SDK_V1_ABI_VERSION
        && sdk->component_count == 2
        && sdk->evaluate_mixture_phase != nullptr;
    std::vector<double> gradient(3);
    std::vector<double> hessian(9);
    epcsaft_mixture_phase_block_result_v1 result{
        sizeof(result),
        EPCSAFT_NATIVE_STATUS_OK_V1,
        3,
        gradient.size(),
        hessian.size(),
        gradient.data(),
        hessian.data(),
    };
    const double amounts[] = {0.4, 0.6};
    const int status = valid
        ? sdk->evaluate_mixture_phase(
            sdk->model_context,
            250.0,
            amounts,
            2,
            1.0e-3,
            &result
        )
        : EPCSAFT_NATIVE_STATUS_INVALID_ARGUMENT_V1;
    std::cout << "component_count=" << (sdk == nullptr ? 0 : sdk->component_count)
              << '\n'
              << "pressure_pa=" << result.pressure_pa << '\n';
    epcsaft_native_model_destroy_v1(model);
    return valid && status == EPCSAFT_NATIVE_STATUS_OK_V1 ? 0 : 4;
}
