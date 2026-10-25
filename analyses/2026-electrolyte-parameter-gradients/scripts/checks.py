"""F5 metrics and F6 independently perturbed held-coordinate charts."""
import math
from epcsaft import Parameters, regression as r
from declarations import parameters, coordinates, holds, row
from run import native_start, manufactured, prediction_columns


def group(name, rows, metric, maximum=None, mean=None, baseline=False, report_only=False):
    out = r.CheckGroup(); out.name=name;out.rows=rows;out.metric=metric
    out.maximum=maximum;out.mean=mean;out.baseline=baseline;out.report_only=report_only
    return out


def metrics(retain):
    p=parameters(); cs=coordinates(p); hs=holds(p); u=native_start(p,cs)
    chart=r.evaluate_coordinates(p,cs,u,hs)
    stored=Parameters(chart.parameters,p.component_ids,p.site_ids)
    templates=[row(p,'mean_ionic_activity_coefficient',m=m) for m in (.1,1,.5,2)]
    templates += [row(p,'solvation_gibbs_energy') for _ in range(2)]
    templates += [row(p,'bulk_relative_permittivity',x={'water':.96,'sodium-cation':.02,'chloride-anion':.02}) for _ in range(2)]
    data=manufactured(stored,cs,u,hs,templates,[.01,-.03,-.005,.015,2,-4,.1,-.2])
    objective=manufactured(stored,cs,u,hs,[row(p,'mean_ionic_activity_coefficient',m=1)])
    invalid=row(p,'solvation_gibbs_energy',form='log_ratio',observed=1)
    data=objective+data+[invalid]
    log=r.CheckMetric.AbsoluteLog10Ratio; diff=r.CheckMetric.AbsoluteDifference
    groups=[group('fail',[1,2],log,.02,.025),group('pass',[3,4],log,.02,.02),
        group('energy',[5,6],diff,5,4),group('epsilon',[7,8],diff,.25,.1),
        group('baseline',[1,2],log,0,0,True),group('unavailable-report',[9],diff,10,report_only=True),
        group('unavailable-guard',[9],diff,10)]
    budget=r.FitControls(maximum_iterations=100,maximum_elapsed_time_seconds=600.0)
    result=r.fit(stored,cs,data,weights=[1]+[0]*9,holds=hs,checks=groups,baseline=p,controls=budget)
    assert result.status=='converged',result.message
    retain('F5','objective-cost',result.final_cost,0,1e-24)
    expected=[(.03,.02,2,'fail',[2]),(.015,.01,4,'pass',[]),(4,3,6,'pass',[]),(.2,.15,8,'fail',[])]
    for c,(maximum,mean,argmax,status,failed) in enumerate(expected):
        check=result.checks[c]
        retain('F5',f'max-{c}',check.maximum,maximum,1e-10); retain('F5',f'mean-{c}',check.mean,mean,1e-10)
        assert (check.argmax,check.status,check.failed_rows)==(argmax,status,failed), (check.name,check.status,check.failed_rows)
    baseline=r.evaluate_checks(p,data,[group('source',[1,2],log,0,0)])[0]; check=result.checks[4]
    retain('F5','baseline-max',check.baseline_maximum,baseline.maximum,1e-10)
    retain('F5','baseline-mean',check.baseline_mean,baseline.mean,1e-10)
    assert result.checks[-1].status==result.checks[-2].status=='unavailable'
    assert result.checks[-1].failed_rows==[9] and result.predictions[9]<0
    # Available negative energy retains a signed difference despite its unavailable logarithmic residual.
    retain('F5','available-negative-energy',result.checks[-1].errors[0],result.predictions[9]-1,1e-10)
    report=r.fit(stored,cs,data,weights=[1]+[0]*9,holds=hs,checks=[groups[-2]],controls=budget)
    assert report.status=='converged' and report.acceptance=='pass' and math.isnan(report.residuals[9])
    return {'status':result.status,'acceptance':result.acceptance,'objective_cost':result.final_cost,
        'groups':[{'name':v.name,'maximum':v.maximum,'mean':v.mean,'status':v.status,'argmax':v.argmax,'failed_rows':v.failed_rows} for v in result.checks],
        'report_only_acceptance':report.acceptance}


def maps(retain):
    p=parameters();cs=coordinates(p);hs=holds(p);source=native_start(p,cs)
    truth=[v+(.05 if abs(v+.05)<=.95 else -.05) if c%2==0 else v+(.1 if abs(v+.1)<=.95 else -.1) for c,v in enumerate(source[:6])]+[source[-1]-1]
    # Nonunit scales and tied C transport must enter each complete optimizer column.
    scales=(.2,.4,.5,.25,.3,.6,2.)
    unit = coordinates(p)
    for c,scale in zip(cs,scales):c.scale=scale
    uncertainty=[]
    for point,theta in enumerate((source,truth)):
        physical = r.evaluate_coordinates(p,unit,theta,hs)
        u=native_start(Parameters(physical.parameters,p.component_ids,p.site_ids),cs); center=r.evaluate_coordinates(p,cs,u,hs);n=len(cs)
        uncertainty += [h.uncertainty for h in center.holds]
        templates=[row(p,kind,m=m,form='difference') for kind in ('mean_ionic_activity_coefficient','osmotic_coefficient') for m in (.1,1,4.5)]
        observations=manufactured(p,cs,u,hs,templates)
        exact=prediction_columns(p,cs,u,hs,observations)
        for c in range(n):
            hstep=1e-4*max(abs(u[c]),1)
            for factor in (10,1,.1):
                step=hstep*factor;plus,minus=list(u),list(u);plus[c]+=step;minus[c]-=step
                high=r.evaluate_coordinates(p,cs,plus,hs);low=r.evaluate_coordinates(p,cs,minus,hs)
                uncertainty += [h.uncertainty for chart in (high,low) for h in chart.holds]
                for ion in range(3):
                    col=center.map[(len(center.physical)-3+ion)*n+c]
                    cd=(high.physical[-3+ion]-low.physical[-3+ion])/(2*step)
                    tol=2e-5*max(abs(center.map[(len(center.physical)-3+j)*n+c]) for j in range(3))+1e-7
                    if factor==1:retain('F6',f'root-{point}-{ion}-{c}',col,cd,tol)
                    retain('F6-step',f'root-{point}-{ion}-{c}-{factor}',col,cd,tol,decisive=False)
                if factor==1:
                    yp=manufactured(p,cs,plus,hs,templates);ym=manufactured(p,cs,minus,hs,templates)
                    column=[exact[i*n+c] for i in range(len(templates))]
                    for i,value in enumerate(column):
                        cd=(yp[i].observed-ym[i].observed)/(2*step)
                        retain('F6',f'complete-{point}-{i}-{c}',value,cd,2e-5*max(abs(v) for v in column)+1e-7)
    retain('F6','root-uncertainty',max(uncertainty),0,1e-10)
    return {'maximum_root_uncertainty_angstrom':max(uncertainty),'scales':scales}
