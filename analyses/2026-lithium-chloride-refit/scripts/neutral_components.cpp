#include "packet-inputs.hpp"
#include "properties/helmholtz_model.hpp"
#include "contributions/helmholtz.hpp"
#include <cppad/cppad.hpp>
#include <cstdio>
#include <cmath>
#include <numeric>
#include <vector>
#include <stdexcept>
int main() {
    auto p=aqueous(), q=extended();
    std::vector<double> n={1/p.molar_mass_kg_per_mol[0],.2,.1,.4,.1,1.35,.05};
    const double total=std::accumulate(n.begin(),n.end(),0.);
    for(auto& v:n) v/=total;
    const auto a=epcsaft::HelmholtzModel(p).evaluate(298.15,1./30000,n,epcsaft::DerivativeOrder::Hessian);
    n.push_back(0);
    const auto b=epcsaft::HelmholtzModel(q).evaluate(298.15,1./30000,n,epcsaft::DerivativeOrder::Hessian);
    std::printf("quantity,row,column,A,B,absolute_difference,relative_difference\n");
    auto compare=[](const char* name,int i,int j,double a,double b) {
        double error=std::abs(a-b), scale=std::max(std::abs(a),std::abs(b));
        double relative=scale==0?0:error/scale;
        std::printf("%s,%d,%d,%.17g,%.17g,%.17g,%.17g\n",name,i,j,a,b,error,relative);
        if(!(relative<=1e-15)) throw std::runtime_error("#781 derivative mismatch");
    };
    compare("residual",0,0,a.value,b.value);
    for(int i=0;i<7;++i) {
        compare("gradient",i,0,a.gradient[i+2],b.gradient[i+2]);
        compare("zero-sum-first",i,0,a.gradient[i+2]-a.gradient[2],b.gradient[i+2]-b.gradient[2]);
        for(int j=0;j<7;++j) {
            compare("hessian",i,j,a.hessian[(i+2)*9+j+2],b.hessian[(i+2)*10+j+2]);
            auto direction=[&](const auto& h,int width){return h[(i+2)*width+j+2]-h[(i+2)*width+2]-h[2*width+j+2]+h[2*width+2];};
            compare("zero-sum-second",i,j,direction(a.hessian,9),direction(b.hessian,10));
        }
    }
    std::FILE* f=std::fopen("build/issue-781/occupied-neutral.csv","w");
    std::fprintf(f,"neutral_mol,neutral_kg_per_mol,F_mol,F_reference_mol,relative_error,molar_value,molar_reference,mu_water,mu_Li,mu_Mg,mu_Cl\n");
    std::vector<double> reference_mu;
    for(double amount:{1.,2.}) for(double mass:{.170338,.340676}) {
        n={1/p.molar_mass_kg_per_mol[0],.2,.1,.4,.1,1.35,.05,amount};
        q.molar_mass_kg_per_mol[7]=mass;
        const double F=epcsaft::cation_mixing_helmholtz(n,q);
        std::vector<CppAD::AD<double>> independent(n.begin(),n.end());
        CppAD::Independent(independent);
        CppAD::ADFun<double> tape(independent,{epcsaft::cation_mixing_helmholtz(independent,q)});
        const auto mu=tape.Jacobian(n);
        if(reference_mu.empty()) reference_mu=mu;
        if(mu!=reference_mu || mu[7]!=0) throw std::runtime_error("neutral changed chemical potentials");
        const auto h=tape.Hessian(n,0);
        for(int i=0;i<8;++i) if(h[i*8+7]!=0 || h[7*8+i]!=0) throw std::runtime_error("neutral Hessian");
        const double N=std::accumulate(n.begin(),n.end(),0.);
        for(auto& v:n) v/=N;
        const auto state=epcsaft::HelmholtzModel(q).state(298.15,30000,n);
        const double relative=std::abs((F+.032688)/F);
        if(relative>1e-15 || std::abs((state.cation_mixing-F/N)/(F/N))>1e-15) throw std::runtime_error("water mass mismatch");
        std::fprintf(f,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",amount,mass,F,-.032688,relative,state.cation_mixing,F/N,mu[0],mu[1],mu[3],mu[5]);
    }
    std::fclose(f);

}
