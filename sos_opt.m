%% Calculate matrix P and K 
syms x1 x2 x3 z1 z2 z3 e1 real;
x = [x1; x2; x3];
z = [z1; z2; z3];
I = [1 0 0; 0 1 0; 0 0 1];
h1 = -x1^2-1;
h2 = x1^4+x1^2;
A = [h1 h2 x1*x3; 0 0 1; 0 0 -1];
B = [0 0 1]';
solver_opt.solver = 'sedumi';

%%
prog1 = sosprogram([x1 x2 x3 z1 z2 z3]);
prog1 = sosdecvar(prog1,[e1]);
[prog1,e2] = sossosvar(prog1,[x]);
[prog1,P] = sospolymatrixvar(prog1,monomials(x,0),[3,3],'symmetric');
[prog1,R] = sospolymatrixvar(prog1,monomials(x,0),[1,1]); % Here R represents inv(R)

%%
K = R*B'*I';
mat2 = P*A'+A*P+K'*B'+B*K;
exp1 = z'*(P - e1*I)*z;
exp2 = -z'*(mat2 + e2*I)*z;
prog1 = sosineq(prog1,e1);
prog1 = sosineq(prog1,exp1);
prog1 = sosineq(prog1,exp2);

%%
prog1 = sossolve(prog1, solver_opt);

%% Calculate Q
Ps = sosgetsol(prog1,P);
Rs = sosgetsol(prog1,R);
Ks = Rs*B'*I';
Q = -(inv(Ps)*A + A'*inv(Ps) + inv(Ps)*B*Ks*inv(Ps));