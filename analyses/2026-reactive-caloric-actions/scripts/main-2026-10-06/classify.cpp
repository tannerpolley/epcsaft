#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
using Row=std::map<std::string,std::string>;
std::vector<std::string> split(const std::string& s,char delimiter) {
    std::vector<std::string> out;std::string field;std::istringstream in(s);
    while(std::getline(in,field,delimiter)) out.push_back(field);
    if(!s.empty()&&s.back()==delimiter) out.push_back("");return out;
}
std::vector<Row> read(const std::string& path,char sep='\t') {
    std::ifstream in(path);if(!in) throw std::runtime_error("missing table: "+path);
    std::string line;std::getline(in,line);const auto headers=split(line,sep);std::vector<Row> out;
    while(std::getline(in,line)) {auto fields=split(line,sep);Row row;
        for(std::size_t i=0;i<headers.size()&&i<fields.size();++i)row[headers[i]]=fields[i];out.push_back(row);}
    return out;
}
double number(Row r,const std::string& k) {
    if(r[k].empty()||r[k]=="None")return std::numeric_limits<double>::quiet_NaN();return std::stod(r[k]);
}
double limit(Row r) {
    if(!r["target"].empty())return number(r,"target");
    if(!r["allowance"].empty())return number(r,"allowance");
    const double a=number(r,"atol"),b=number(r,"rtol"),m=std::abs(number(r,"model")),v=std::abs(number(r,"reference"));
    if(r["rule"]=="symmetric_close"||r["rule"]=="strict_symmetric_close")return std::max(a,b*std::max(m,v));
    if(r["rule"]=="reference_scale"||r["rule"]=="strict_reference_scale")return a+b*v;
    if(r["rule"]=="absolute"||r["rule"]=="strict_absolute")return a;
    if(r["rule"]=="reference_or_one"||r["rule"]=="strict_reference_or_one")return a+b*std::max(1.,v);
    return std::numeric_limits<double>::quiet_NaN();
}
int main(int argc,char** argv) {try {
    if(argc!=2)return 2;const std::string out=argv[1];std::cout<<std::setprecision(17);
    struct Count {int pass=0,fail=0,observations=0;double maximum=0;std::string maximum_id;};std::map<std::string,Count> counts;
    const auto rows=read(out+"/comparisons.tsv");
    std::ofstream critical(out+"/original-order.tsv");critical<<std::setprecision(17);
    critical<<"order\trow\texact\treference\trefinement\ttarget\terror\terror_over_target\trefinement_over_quarter_target\n";
    std::vector<Row> selected;double exact=0;
    for(auto r:rows) {
        auto& c=counts[r["test_id"]];if(r["verdict"]=="True")++c.pass;else if(r["verdict"]=="False")++c.fail;else ++c.observations;
        if(!r["verdict"].empty()) {double ratio=number(r,"error_over_target");
            const auto id=r["comparison_id"];
            const bool agreement=r["rule"]!="abs(V_row-V_anchored)>10*target" && (r["criterion"]!="custom" || r["test_id"]=="native.properties.reference-limits" || id.find("direct/action/")==0 || id.find("pressure/action/")==0 || id.find("electrolyte-coordinate/")==0);
            if(!std::isfinite(ratio)&&agreement)ratio=std::abs(number(r,"model")-number(r,"reference"))/limit(r);
            if(std::isfinite(ratio)&&ratio>c.maximum){c.maximum=ratio;c.maximum_id=r["comparison_id"]; }}
        const auto id=r["comparison_id"];
        if(id.find("F6/N-bubble-P-fix/A2/T+water-feed/total-cv_G-V/")==0) {if(!r["verdict"].empty())exact=number(r,"model");else if(r["resolved"]=="1")selected.push_back(r);}
    }
    std::ofstream summary(out+"/qualification.tsv");
    summary<<std::setprecision(17)<<"test_id\tpass\tfail\tobservations\tmax_error_over_target\tmaximum_id\n";
    for(auto r:rows) {
        const auto id=r["comparison_id"];
        if(id.rfind("F1/",0)!=0 && id.rfind("F3/",0)!=0 && id.rfind("F4/",0)!=0 && id.rfind("F5/",0)!=0 && id.rfind("reference-calorics/F4/",0)!=0 && id.rfind("reference-calorics/F5/",0)!=0) continue;
        const auto group=id.rfind("reference-calorics/",0)==0?id.substr(19,2):id.substr(0,2);
        auto& c=counts["criterion."+group];
        if(r["verdict"]=="True")++c.pass;else if(r["verdict"]=="False")++c.fail;else ++c.observations;
        if(r["verdict"].empty() || r["rule"]=="abs(V_row-V_anchored)>10*target")continue;
        double ratio=number(r,"error_over_target");
        if(!std::isfinite(ratio))ratio=std::abs(number(r,"model")-number(r,"reference"))/limit(r);
        if(std::isfinite(ratio)&&ratio>c.maximum){c.maximum=ratio;c.maximum_id=id;}
    }
    for(auto [id,c]:counts) summary<<id<<'\t'<<c.pass<<'\t'<<c.fail<<'\t'<<c.observations<<'\t'<<c.maximum<<'\t'<<c.maximum_id<<'\n';
    for(auto [id,c]:counts)std::cout<<id<<" pass="<<c.pass<<" fail="<<c.fail<<" observations="<<c.observations<<" max_error_over_target="<<c.maximum<<" maximum_id="<<c.maximum_id<<"\n";
    bool original=false;
    for(auto r:selected) {
        const double ref=number(r,"Richardson"),target=number(r,"target"),refinement=number(r,"refinement"),error=std::abs(exact-ref);
        const auto id=r["comparison_id"];const bool is_original=id.find("/original")!=std::string::npos;
        critical<<(is_original?"original":"reversed")<<'\t'<<id<<'\t'<<exact<<'\t'<<ref<<'\t'<<refinement<<'\t'<<target<<'\t'<<error<<'\t'<<error/target<<'\t'<<4*refinement/target<<'\n';
        if(is_original)original=error<=target && refinement<=target/4;
    }
    std::cout<<"original_order_pass="<<original<<" selected_rows="<<selected.size()<<"\n";
    std::ofstream moved(out+"/native-row-movements.tsv");moved<<std::setprecision(17);
    moved<<"source\tid\toccurrence\tunit\tpre_scaling\tbranch\tmain\tintegrated\tchange_from_branch\tchange_from_main\tchange_scaling\tnative_reference\tnative_error\tnative_target\terror_over_target\tnative_verdict\tcause\tmethod\trule\treference_comparison_id\tcriterion\n";
    std::map<std::string,int> causes;int changed_branch=0,changed_main=0,failed=0,missing=0,unchanged=0;
    for(auto r:read(out+"/movement-inputs.tsv")) {
        const double base=number(r,"pre_scaling"),branch=number(r,"branch"),main=number(r,"main"),value=number(r,"new");
        const bool density=std::isfinite(branch)&&value!=branch;
        const bool main_changed=std::isfinite(main)&&value!=main;
        const bool scaling=std::isfinite(base)&&std::isfinite(branch)&&base!=branch;
        if(density)++changed_branch;if(main_changed)++changed_main;
        if(!density&&!main_changed) {++unchanged;continue;}
        const std::string cause=density?(scaling?"both":"#222 density"):(scaling?"#225 scaling":"unclassified");
        ++causes[cause];
        r["model"]=r["new"];r["target"]=r["custom_target"];
        const double reference=number(r,"reference"),target=limit(r),error=r["metric_error"].empty()?std::abs(value-reference):number(r,"metric_error");
        if(r["native_verdict"]!="True"||(r["rule"]!="abs(V_row-V_anchored)>10*target"&&std::isfinite(target)&&std::isfinite(error)&&error>target))++failed;
        if(!std::isfinite(reference)||r["method"].empty())++missing;
        moved<<r["source"]<<'\t'<<r["id"]<<'\t'<<r["occurrence"]<<'\t'<<r["unit"]<<'\t'<<base<<'\t'<<branch<<'\t'<<main<<'\t'<<value<<'\t'<<value-branch<<'\t'<<value-main<<'\t'<<branch-base<<'\t'<<reference<<'\t'<<error<<'\t'<<target<<'\t'<<error/target<<'\t'<<r["native_verdict"]<<'\t'<<cause<<'\t'<<r["method"]<<'\t'<<r["rule"]<<'\t'<<r["reference_comparison_id"]<<'\t'<<r["criterion"]<<'\n';
    }
    std::ofstream classification(out+"/classification-summary.tsv");classification<<"cause\tcount\n";
    for(auto [cause,count]:causes) {classification<<cause<<'\t'<<count<<'\n';std::cout<<cause<<" count="<<count<<'\n';}
    std::cout<<"changed_vs_branch="<<changed_branch<<" changed_vs_main="<<changed_main<<" unchanged_in_both="<<unchanged<<" failed="<<failed<<" missing_reference="<<missing<<'\n';
    double max_moved=0;std::string moved_id;
    for(auto r:read(out+"/native-row-movements.tsv")) {
        if(r["rule"]=="abs(V_row-V_anchored)>10*target")continue;
        const auto ratio=number(r,"error_over_target");
        if(std::isfinite(ratio)&&ratio>max_moved){max_moved=ratio;moved_id=r["id"];}
    }
    std::ofstream movement_summary(out+"/movement-summary.tsv");
    movement_summary<<std::setprecision(17)<<"max_error_over_target\tmaximum_id\tfailed\tmissing_reference\n"
        <<max_moved<<'\t'<<moved_id<<'\t'<<failed<<'\t'<<missing<<'\n';
    return failed||missing||causes["unclassified"]||!original?1:0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
