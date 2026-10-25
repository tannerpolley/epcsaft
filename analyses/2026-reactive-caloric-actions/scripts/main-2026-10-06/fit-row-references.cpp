#include "regression_j2_inventory_checks.cpp"
using namespace regression_checks;
#include <map>
std::map<std::string,double> recorded_rows() {
    std::ifstream f("build/current-calorics/main-2026-10-06/native-rows/native-rows.tsv");
    std::map<std::string,double> rows;std::string line;
    while(std::getline(f,line)) if(line.rfind("regression_checks R2-J2\t",0)==0) {
        const auto a=line.find('\t'),b=line.find('\t',a+1);rows[line.substr(a+1,b-a-1)]=std::stod(line.substr(line.rfind('\t')+1));
    }return rows;
}
int main(){try {
    const auto rows=recorded_rows();const auto value=[&](const std::string& id){return rows.at("r4/recovery/"+id);};
    const auto known=manufactured();auto request=recovery_request(known,value("epsilon"),value("kij"));
    const auto u=request.coordinates.start(request.parameters);const HelmholtzModel model(request.parameters,request.coordinates.active_parameters());
    double initial=0;
    const auto start=recovery_request(known,148,0);const HelmholtzModel start_model(start.parameters);
    for(const auto& o:start.observations) {const auto p=predict(start_model,o);CHECK(p.residual);initial+=0.5 * *p.residual * *p.residual;}
    CHECK(agree("r4/recovery/initial_cost","cost","1",value("initial_cost"),initial,0,1e-12,"original R3 cost criterion","native direct weighted residual sum","r4"));
    for(std::size_t i=0;i<known.observations.size();++i) {
        const auto& o=known.observations[i];const double pred=value("prediction/"+std::to_string(i));
        const double residual=o.form==ResidualForm::LogRatio?std::log(pred/o.observed):(pred-o.observed)/o.scale;
        const double bound=o.form==ResidualForm::LogRatio?o.observed*std::expm1(std::sqrt(2e-14)):o.scale*std::sqrt(2e-14);
        const auto id="r4/recovery/prediction/"+std::to_string(i);
        CHECK(recording_support::compare_custom(id.c_str(),"manufactured prediction","native",pred,o.observed,
            fields({field("target",bound),field("half_squared_residual",0.5*residual*residual)}),"Numerical","native manufactured truth",
            "necessary consequence of original total fit cost criterion","0.5 residual^2<=1e-14",0.5*residual*residual<=1e-14));
        for(std::size_t c=0;c<2;++c) {
            const auto shifted=[&](double delta){auto q=u;q[c]+=delta;const HelmholtzModel shifted_model(request.coordinates.apply(request.parameters,q));return *predict(shifted_model,o).value;};
            auto rungs=ladder(shifted,0.005);const double scale=request.coordinates.coordinates[c].scale;
            rungs.coarse/=scale;rungs.fine/=scale;const double slope=residual_slope(o,pred);
            const double J=value("optimizer_jacobian/"+std::to_string(2*i+c));
            const auto id="r4/recovery/optimizer_jacobian/"+std::to_string(2*i+c);
            const bool ok=ladder_agrees(id+"/physical-gradient","dy_dtheta",J/(slope*scale),rungs,"original R2-J1 0.005 optimizer step and FIRST_ORDER criterion","rebuilt-model gradient","r2-j1");
            const double ref=rungs.fine*slope*scale,allowance=FIRST_ORDER*std::max(1.,std::abs(rungs.fine))*std::abs(slope*scale);
            CHECK(recording_support::compare_custom(id.c_str(),"optimizer residual Jacobian","1",J,ref,
                fields({field("target",allowance),field("physical_scale",scale),field("residual_slope",slope)}),"Numerical","R2-J1 native rebuilt-model ladder",
                "existing physical-gradient criterion with exact residual and coordinate conversion","original R2-J1 ladder criterion",ok));
        }
    }
    double gram[4]={},cov[4]={};for(std::size_t a=0;a<2;++a)for(std::size_t b=0;b<2;++b)for(std::size_t i=0;i<5;++i)
        gram[2*a+b]+=value("optimizer_jacobian/"+std::to_string(2*i+a))*value("optimizer_jacobian/"+std::to_string(2*i+b));
    for(std::size_t k=0;k<4;++k)cov[k]=value("covariance/"+std::to_string(k));
    double worst=0;for(std::size_t a=0;a<2;++a)for(std::size_t b=0;b<2;++b)
        worst=std::max(worst,std::abs(cov[2*a]*gram[b]+cov[2*a+1]*gram[2+b]-(a==b?1.:0.)));
    const double det=gram[0]*gram[3]-gram[1]*gram[2],inverse[4]={gram[3]/det,-gram[1]/det,-gram[2]/det,gram[0]/det};
    for(std::size_t k=0;k<4;++k) {
        const auto id="r4/recovery/covariance/"+std::to_string(k);
        CHECK(recording_support::compare_custom(id.c_str(),"optimizer covariance","1",cov[k],inverse[k],
            fields({field("target",1e-10),field("Gram_identity_error",worst)}),"Identity","R3 inverse Gram matrix",
            "covariance times independently qualified Jacobian Gram is identity","max_abs(Cov JtJ-I)<=1e-10",worst<=1e-10));
        const double scales[2]={10,.1};const auto physical="r4/recovery/physical_covariance/"+std::to_string(k);
        CHECK(agree(physical,"physical covariance","native",value("physical_covariance/"+std::to_string(k)),scales[k/2]*cov[k]*scales[k%2],0,1e-12,"R3 covariance coordinate identity","Cov_theta=D Cov_u Dt","r3"));
    }
    CHECK(agree("r4/recovery/variance_factor","variance factor","1",value("variance_factor"),2*value("final_cost")/3,0,1e-12,"five rows, rank two; original R3 criterion","sum(w r^2)/(n-rank)","r3"));
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
