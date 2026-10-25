#include "regression_checks_shared.hpp"
using namespace regression_checks;
#include <fstream>
#include <sstream>
int main() {try {
    const auto target=manufactured().observations.front().observed;
    std::ifstream in("build/current-calorics/main-2026-10-06/native-rows/native-rows.tsv");
    std::string line;double value=0;bool found=false;
    while(std::getline(in,line))if(line.find("regression_checks R2-J2\tr4/recovery/prediction/0\t")==0){value=std::stod(line.substr(line.rfind('\t')+1));found=true;}
    CHECK(found);const double contribution=0.5*std::pow(std::log(value/target),2);
    const double bound=target*std::expm1(std::sqrt(2e-14));
    CHECK(recording_support::compare_custom("r4/recovery/prediction/0","bubble pressure","Pa",value,target,
        recording_support::fields({recording_support::field("target",bound),recording_support::field("squared_log_residual_half",contribution),recording_support::field("original_total_cost_tolerance",1e-14)}),
        "Numerical","native manufactured pressure at truth parameters","A moved prediction against the unchanged manufactured target and original fit cost",
        "0.5*log(prediction/reference)^2<=1e-14 (necessary consequence of original total-cost criterion)",contribution<=1e-14));
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
