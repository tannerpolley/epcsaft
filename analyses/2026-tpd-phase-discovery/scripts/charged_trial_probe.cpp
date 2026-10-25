#include "properties_support.hpp"
#include "equilibrium_support.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/ipopt.hpp"
#include "core/numerical_policy.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <numeric>
using namespace epcsaft::equilibrium;
void vec(const std::string& n,const std::vector<double>& v){std::cout<<n; for(auto x:v)std::cout<<" "<<x;std::cout<<"\n";}
void dump(const std::string& tag,const CompiledSystem& c,const EquilibriumResult& r){
 std::cout<<"RESULT "<<tag<<" success="<<r.success<<" status="<<r.solver_status<<" message="<<r.message<<" objective="<<r.objective<<" residual="<<r.convergence.residual<<" criterion="<<c.problem.solver_options.value(epcsaft::policy::row)<<" inertia="<<r.reduced_hessian_minimum<<" iterations="<<r.iterations<<" evaluations="<<r.evaluations<<" hessians="<<r.hessian_calls<<"\n";
 vec("coordinates",r.coordinates);vec("raw_rows",r.residuals);vec("classified_rows",r.classified_residuals);vec("multipliers",r.multipliers);vec("composition",r.mole_fractions);vec("density",r.molar_densities);vec("balance",c.balance);for(auto b:c.B)vec("B",b);
 if(!r.multipliers.empty()){double d=0;for(size_t i=0;i<c.rank;++i)d+=r.multipliers[i]*c.balance[i];std::cout<<"dual="<<d<<" identity_defect="<<r.objective-d<<"\n";}
}
int main(int argc,char**){std::cout<<std::setprecision(17)<<std::unitbuf;try{
 auto p=equilibrium_support::permute_parameters(properties_support::aqueous_electrolyte(),{1,2,0});epcsaft::HelmholtzModel m(p);
 auto ph=equilibrium_support::phase("L",PhaseKind::Liquid);ph.support={0,1,2};std::vector<std::vector<double>> t;
 auto c=compile_stability_trial(m,298.15,1e5,ph,{.01,.01,.98},{0,0,0},t);
 if(argc>1){
 epcsaft::core::Capture capture;epcsaft::core::Scope scope("raw_diagnostic");
 auto start=boundary_start(c);vec("boundary_start",start);
 Ipopt::SmartPtr<GibbsMinimizer> nlp=new GibbsMinimizer(c,start,&scope);
 const auto stops=run_ipopt(nlp,c.problem.solver_options.value(epcsaft::policy::row)*(c.feasibility_cell()?1.0:epcsaft::policy::EQUILIBRIUM_INNER_TOLERANCE_FACTOR.value),(int)epcsaft::policy::EQUILIBRIUM_MAX_ITERATIONS.value,0,true,1/c.feed_total);
 std::cout<<"RAW initialize="<<application_return_name(stops.first)<<" returned="<<application_return_name(stops.second)<<" solver="<<solver_status_name(nlp->status())<<" iterations="<<nlp->iterations()<<" evaluations="<<nlp->evaluations()<<" hessians="<<nlp->hessian_calls()<<" prep_error="<<nlp->preparation_error().value_or("none")<<"\n";
 auto raw=nlp->solution();vec("raw_coordinates",raw);vec("optimizer_multipliers",nlp->constraint_multipliers());
 std::vector<double> lo(raw.size()),hi(raw.size()),gl(c.rows.size()),gu(c.rows.size());nlp->get_bounds_info(raw.size(),lo.data(),hi.data(),c.rows.size(),gl.data(),gu.data());vec("lower_bounds",lo);
 auto at=c.prepare(raw);auto ev=c.evaluate_at(at);vec("raw_composition",at.composition[0]);vec("raw_volumes",at.volume);vec("raw_eta",at.eta);vec("chemical_without_plane",ev.values);
 std::cout<<"feasible_solvent_derivative="<<ev.values[2]-.5*(ev.values[0]+ev.values[1])<<" objective="<<ev.objective<<" floor="<<ev.floor<<"\n";
 auto full=raw;full.push_back(ev.values[0]);full.push_back(ev.values[1]);
 auto rows=c.evaluate(full).values;auto j=c.jacobian(full);using Dense=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
 Eigen::VectorXd delta=Eigen::Map<const Dense>(j.data(),full.size(),full.size()).fullPivLu().solve(-Eigen::Map<const Eigen::VectorXd>(rows.data(),rows.size()));
 vec("full_coordinates",full);vec("full_raw_rows",rows);vec("full_classified_rows",c.classified_residuals(c.prepare(full),rows));std::vector<double> d(delta.data(),delta.data()+delta.size());vec("Newton_delta",d);
 for(double alpha:{1.,.5,.25}){auto x=full;for(size_t i=0;i<x.size();++i)x[i]+=alpha*d[i];auto a=c.prepare(x);std::cout<<"POLISH alpha="<<alpha<<" solvent="<<a.composition[0][2]<<" cutoff="<<epcsaft::policy::EQUILIBRIUM_BOX_FLOOR.value<<"\n";vec("polish_composition",a.composition[0]);}
 for(double water:{1e-8,1e-12,1e-20}){std::vector<double> w{(1-water)/2,(1-water)/2,water};auto x=raw;for(size_t i=0;i<3;++i)x[i]=std::log(w[i]);const auto root=epcsaft::solve_density(m,298.15,1e5,w,epcsaft::DensityBranch::Liquid,std::nullopt);x[4]=std::log(m.packing_fraction_coefficient(298.15,w)*root.molar_density);auto a=c.prepare(x);auto e=c.evaluate_at(a);std::cout<<"DIRECT water="<<water<<" density="<<root.molar_density<<" objective="<<e.objective<<" feasible_solvent_derivative="<<e.values[2]-.5*(e.values[0]+e.values[1])<<" pressure_row="<<e.values[3]<<" floor="<<e.floor<<"\n";}
 return 0;
 }
 {epcsaft::core::Capture capture;auto r=solve_compiled(c);dump("artificial-zero-plane",c,r);equilibrium_support::print_trace(capture.trace());}
 auto q=equilibrium_support::flash(3);q.phases={ph};q.feed.values={.01,.01,.98};q.temperature=FreeValue::fixed(298.15);q.pressure=FreeValue::fixed(1e5);q.solver_options.values[epcsaft::policy::row]=1e-10;
 auto declared=compile_problem(m,q);auto candidate=solve_compiled(declared);dump("declared-liquid",declared,candidate);
 if(candidate.success){auto a=declared.prepare(candidate.coordinates);auto shift=a.standard.derivative[0];for(size_t i=0;i<shift.size();++i)for(size_t k=0;k<declared.rank;++k)shift[i]-=declared.B[k][i]*candidate.multipliers[k];vec("frozen_shift",shift);
 auto actual=compile_stability_trial(m,298.15,1e5,ph,{.01,.01,.98},shift,t);
 {epcsaft::core::Capture capture;auto r=solve_compiled(actual);dump("actual-candidate-plane",actual,r);equilibrium_support::print_trace(capture.trace());}
 }
}catch(const std::exception& e){std::cout<<"EXCEPTION "<<e.what()<<"\n";return 2;}}
