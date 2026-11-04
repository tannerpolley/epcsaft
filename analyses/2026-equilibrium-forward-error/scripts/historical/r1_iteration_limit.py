from r1_common import *
from tests.packets import parameters_json
m=Mixture(Parameters.from_json(parameters_json('gross-2001-methane-ethane',3),components=('methane','ethane')))
p=Problem(T=300.,P=1e5,feed=Amounts({'methane':.5,'ethane':.5}),phases=[Phase('L',kind='liquid')],row_tolerance=2.)
opt=OUT/'r1-iteration-limit.opt'
opt.write_text('max_iter 0\noutput_file '+str(OUT/'r1-iteration-limit-ipopt.log')+'\nfile_print_level 12\n')
os.environ['EPCSAFT_IPOPT_OPTIONS_FILE']=str(opt)
with trace(diagnostics=True) as c:r=solve_equilibrium(m,p)
os.environ.pop('EPCSAFT_IPOPT_OPTIONS_FILE')
save('r1-iteration-limit',dict(result=summary(r),events=events(c)))
print(summary(r))
