syms x1 x2 x3 u e1 real;
Z = [x1 x2 x3]';
v = [x1 x2 x3 u]';
h1=-x1^2 - 1;
h2=x1^4 + x1^2;
A = [h1 h2 x1*x3; 0 0 1; 0 0 -1];
B = [0; 0; 1];
C = [0 1 0];
N = [0 1 0 0; 0 0 0 1]; % N represents [C 0; 0 I] matrix
I = [1 0 0; 0 1 0; 0 0 1];
solver_opt.solver = 'sedumi';

%%
prog = sosprogram([x1 x2 x3 u]);
% Here P represents the matrix inv(P) in the paper. That means we are
% directly using inv(P) as the varible matrix
[prog,P] = sospolymatrixvar(prog, monomials(Z,0), [3 3], 'symmetric');
% W represents [R S; S' Q] matrix
[prog, W] = sospolymatrixvar(prog, monomials(Z,0), [2 2], 'symmetric');

%%
T1 = [A'*P+P*A P*B; B'*P' 0];
exp1 = Z'*(P)*Z;
exp2 = v'*( N'*W*N - T1 )*v;
prog = sosineq(prog,exp1);
prog = sosineq(prog,exp2);

%%
prog = sossolve(prog, solver_opt);