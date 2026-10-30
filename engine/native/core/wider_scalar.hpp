#pragma once
#include <cmath>
#include <limits>
#include <vector>
#ifdef _WIN32
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/eigen.hpp>
#endif

namespace epcsaft {
#ifdef _WIN32
static_assert(BOOST_VERSION == 108900);
using WiderScalar = boost::multiprecision::cpp_bin_float_double_extended;
static_assert(std::numeric_limits<WiderScalar>::radix == 2 && std::numeric_limits<WiderScalar>::digits == 64);
static_assert(std::numeric_limits<WiderScalar>::min_exponent == -16382 && std::numeric_limits<WiderScalar>::max_exponent == 16383);
namespace wide_math {
using std::abs;
using std::exp;
using std::expm1;
using std::isfinite;
using std::log;
using std::log1p;
using std::sqrt;
using std::pow;
using boost::multiprecision::abs;
using boost::multiprecision::exp;
using boost::multiprecision::expm1;
using boost::multiprecision::isfinite;
using boost::multiprecision::log;
using boost::multiprecision::log1p;
using boost::multiprecision::sqrt;
using boost::multiprecision::pow;
}
#else
using WiderScalar = long double;
namespace wide_math = std;
#endif
template <typename Iterator>
std::vector<double> narrow_values(Iterator first, Iterator last) {
#ifdef _WIN32
    std::vector<double> result;
    for (; first != last; ++first) result.push_back(static_cast<double>(*first));
    return result;
#else
    return {first, last};
#endif
}
}
