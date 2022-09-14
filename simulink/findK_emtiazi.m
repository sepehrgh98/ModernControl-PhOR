clear all
close all
clc

A = [1 2 0; 0 -1 3; 0 1 -1];
b = [1 3 1]';
C = [0 1 0];
d = 0;
[numOL,denOL] = ss2tf(A,b,C,d);
h = tf(numOL,denOL)
roots(denOL);
Wc = ctrb(A,b);
disp('rank of controllability mtrix = ')
l=rank(Wc)
%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp(' place in MATLAB')
pc=[-1,-2,-3];
K0 = place(A, b, pc)
toc
%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp(' Bass-Gura Formula')
alpha = [1 6 11 6];
a = denOL;
Omega = [1 a(2) a(3);0 1 a(2);0 0 1];
K1 = (alpha(2:4)-a(2:4))*inv(Omega)*inv(Wc)
toc
%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp(' Ackermann Formula')
alpha = [1 6 11 6];
alpha_of_A = zeros(3,3);
for i=1:4
alpha_of_A = alpha_of_A + alpha(i)*A^(4-i);
end
K2 = [0 0 1]*inv(Wc)*alpha_of_A
toc
%%%%%%%%%%%%%%%%%%%%%%%%%



