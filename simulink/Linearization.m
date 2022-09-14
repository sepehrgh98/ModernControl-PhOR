a1 = 6.11e-3;
a2 = -2.89e-3;
a3 = -4.24e-3;
a4 = 3.01e-3;
a5 = 2.05e-3;                              
a6 = 1.92e-3;
a7 = 1.60e-3;
a8 = -8.32e-3;
bw = 30;
cw = 1500;
%Non Linear :
syms n1 n2 n3 n4 u1 u2 d1 d2 s t
M = [a1 + a2*cos(2*n2) + a3*sin(2*n2) + a4*cos(n2) + a5*sin(n2) 0; 0 a6];
M_inv = inv(M);
V = [-2*a1*n3*n4*sin(2*n2) + 2*a3*n3*n4*cos(2*n2) + a4*n3*n4*cos(n2) - a5*n3*n4*sin(n2) ; 
    2*a2*(n3)^2*cos(n2)*sin(n2) - a3*(n3)^2*cos(n2) - 0.5*a4*(n3)^2*cos(n2) + 0.5*a5*(n3)^2*sin(n2) + a7*sin(n2) + a8*cos(n2)];
f = bw*[n3;n4] + cw*[n1;n2];
u = [u1 ; u2];
d = [d1 ; d2];
C = [1 0 0 0 ; 0 1 0 0 ];
x1_dot = [n3 ; n4];
x2_dot = (-M_inv)*V - (M_inv)*f + (M_inv)*u + (M_inv)*d;
n_dot = [x1_dot ; x2_dot];

%Linearization:
SOLV = solve(n_dot == 0 ,[n1 n2 n3 n4 u1 u2 d1 d2]);
EQ = [SOLV.n1  SOLV.n2  SOLV.n3  SOLV.n4 SOLV.u1 SOLV.u2 SOLV.d1 SOLV.d2];
var = [n1 n2 n3 n4 u1 u2 d1 d2];
JAC_A = jacobian(n_dot,[n1 n2 n3 n4]);
A = subs(JAC_A,var,EQ);
A= double(A);
JAC_B = jacobian(n_dot,u);
B = subs(JAC_B, var,EQ);
B = vpa(B,3);
B= double(B);
JAC_Bd = jacobian(n_dot,d);
Bd = subs(JAC_Bd,var,EQ);

%Jordan matrix:
[T_A, J_A] = jordan(A);
vpa(T_A);
vpa(J_A);

%State Transsmission Matrix:
Phi = ilaplace(inv(s*eye(4)-A));

%Zero Input 
x0 = [0.5; 0.5 ; 0 ; 0 ];
X_zir = Phi * x0;
vpa(X_zir , 3);

%Checking Controllability & Observability:
Controllability = ctrb(A,B);
if rank(Controllability) == 4
    disp('Controllable !! ')
end
Observability = obsv(A,C);
if rank(Observability) == 4
    disp('Observable !! ')
end

%State Feedback
syms s
pc=[-4-0.86j -4+0.86j -80 -90]
k = place(A,B,pc);

%static pre compensator
Gcl=inv(C*inv(-A+B*k)*B);

%Observer
po = [-20 -25 -400 -360];
l = place(A',C',po)';

%Lyapanov Stability first method
stability_check=eig(A)

%Lyapanov Stability second method
q=eye(4);
pm=lyap(A,q);
pm_check=eig(pm);

%WTF
syms n1 n2 n3 n4 x1 x2 x1_dot x2_dot
M=[a1 + a2*cos(2*n2) + a3*sin(2*n2) + a4*cos(n2) + a5*sin(n2) 0; 0 a6];
M_inv=inv(M)
V=[-2*a1*n3*n4*sin(2*n2)+2*a3*n3*n4*cos(2*n2)+a4*n3*n4*cos(n2)-a5*n3*n4*sin(n2);
    2*a2*(n3)^2*cos(n2)*sin(n2)-a3*(n3)^2*cos(n2)-0.5*a4*(n3)^2*cos(n2)+0.5*a5*(n3)]
f=bw*[n3;n4]+cw*[n1;n2];
u=[u1;u2];
d=[d1;d2];
x1=[n1;n2];
x2=[n3;n4];
x1_dot=[n3;n4];
x2_dot=-M_inv*V - M_inv*f+M_inv*u+M_inv*d;
v_dot=2*x1.'*x1_dot+2*x2.'*x2_dot;

%COST_FUNCTION
syms p11 p12 p13 p14 p21 p22 p23 p24 p31 p32 p33 p34 p41 p42 p43 p44
R = eye(2);
Q = 1000*[1 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
P = [p11 p12 p13 p14 ;
    p21 p22 p23 p24;
    p31 p32 p33 p34;
    p41 p42 p43 p44];
E = A'*P + P*A - P*B*inv(R)*B'*P == -Q ;
sP = solve(E , P);
i = 1;
P_ = [sP.p11(i) sP.p12(i) sP.p13(i) sP.p14(i) ;
      sP.p21(i) sP.p22(i) sP.p23(i) sP.p24(i);
      sP.p31(i) sP.p32(i) sP.p33(i) sP.p34(i);
      sP.p41(i) sP.p42(i) sP.p43(i) sP.p44(i)];
P_ = vpa(P_,3);
k1= inv(R)*B'*P_;
k1 = vpa(k,3);
k1 = double(k1)