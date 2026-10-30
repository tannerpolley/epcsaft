#include "core/wider_scalar.hpp"
#include <iostream>
#include <iomanip>
int main(){using W=epcsaft::WiderScalar;using L=std::numeric_limits<W>;std::cout<<std::setprecision(25)<<"radix="<<L::radix<<" digits="<<L::digits<<" min_exponent="<<L::min_exponent<<" max_exponent="<<L::max_exponent<<" epsilon="<<L::epsilon()<<" denorm_min="<<L::denorm_min()<<'\n';}
