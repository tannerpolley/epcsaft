#pragma once

#include <string>

#include "flash.hpp"

namespace epcsaft_engine::equilibrium {

[[nodiscard]] std::string flash_result_to_json(const FlashResult& result);

}  // namespace epcsaft_engine::equilibrium
