g3 1 1 0	# problem unknown
 6 4 1 0 4 	# vars, constraints, objectives, ranges, eqns
 2 1 0 0 0 0	# nonlinear constrs, objs; ccons: lin, nonlin, nd, nzlb
 0 0	# network constraints: nonlinear, linear
 2 1 1 	# nonlinear vars in constraints, objectives, both
 0 1 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0 	# discrete variables: binary, integer, nonlinear (b,c,o)
 9 2 	# nonzeros in Jacobian, obj. gradient
 24 13	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
F0 1 -1 epcsaft_neutral_pressure
C0	#eos_pressure
o2	#*
n5e-06
f0 6	#pressure
n400.0
o2	#*
n100.0
v1	#fs.unit.r
h5:2.002
h6:3.6184
h6:208.11
h8:0.044096
C1	#fs.unit.volume_balance
o16	#-
o2	#*
v1	#fs.unit.r
v0	#fs.unit.q
C2	#fs.unit.flow_balance
n0
C3	#fs.unit.pressure_balance
n0
O0 0	#objective
o2	#*
n0.5
o5	#^
v0	#fs.unit.q
n2
x6	# initial guess
0 0.2	#fs.unit.q
1 0.2	#fs.unit.r
2 0.04	#fs.unit.f_out
3 0.5	#fs.unit.p_out
4 0.04	#fs.unit.f_in
5 0.5	#fs.unit.p_in
r	#4 ranges (rhs's)
4 0.0	#eos_pressure
4 0.0	#fs.unit.volume_balance
4 0.0	#fs.unit.flow_balance
4 0.0	#fs.unit.pressure_balance
b	#6 bounds (on variables)
0 0.1 1.0	#fs.unit.q
0 0.1 1.0	#fs.unit.r
0 0.0 1.0	#fs.unit.f_out
0 5e-06 1.0	#fs.unit.p_out
0 0.0 1.0	#fs.unit.f_in
0 5e-06 1.0	#fs.unit.p_in
k5	#intermediate Jacobian column lengths
1
3
5
7
8
J0 2	#eos_pressure
1 0
3 -1.0
J1 3	#fs.unit.volume_balance
0 0
1 0
2 1
J2 2	#fs.unit.flow_balance
2 -1
4 1
J3 2	#fs.unit.pressure_balance
3 -1
5 1
G0 2	#objective
0 0
2 -1
