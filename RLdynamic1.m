
function xout = RLdynamic1(t,x)
global dphi1x;
global dphi2x;
global W1;
global W2;
global W3;
global W4;

t 

A = 2;
B1 = 1;
B2 = 3;

R11 = 0.64;
R12 = 1;
R21 = 1;
R22 = 1;

Q1 = 9;
Q2 = 9;

% global t;
a1 = 0.05;
a2 = 0.05;
a3 = 0.05;
a4 = 0.05;

% t=0;
x1 = x(1);
X = [x(1)]';

W1 = [x(2)]'
W2 = [x(3)]'
W3 = [x(4)]'
W4 = [x(5)]'

phi1x = [x1^2]';
phi2x = [x1^2]';
dphi1x = [2*x1];
dphi2x = [2*x1];

n = size(A,1);
num = n*(n+1)/2;

f = A*X;
g = B1;
k = B2;

u3 = -0.5*inv(R11)*g'*dphi1x'*W3;
d4 = -0.5*inv(R22)*k'*dphi2x'*W4;

delta3 =  dphi1x*f + dphi1x*g*u3 +dphi1x*k*d4;
delta4 =  dphi2x*f + dphi2x*g*u3 +dphi2x*k*d4;

delta3_bar = delta3/(delta3'*delta3+1);
delta4_bar = delta4/(delta4'*delta4+1);

ms1 = delta3'*delta3+1;
ms2 = delta4'*delta4+1;

m1 = delta3/ms1^2;
m2 = delta4/ms2^2;

% check the reasonabilty of selected parameters by RL  
% q1 = rand(n,n);
% q2 = rand(n,n);

q1 = 9*eye(n);
q2 = 9*eye(n);

% F1 = rand(num,1);
% F3 = rand(num,1);

F1 = ones(num,1);
F3 = ones(num,1);

% F2 = rand(num,num);
% F4 = rand(num,num);

F2 = 10*eye(num);
F4 = 10*eye(num);
% check the reasonabilty of selected parameters by RL  

D1_bar = dphi1x*g*inv(R11)'*g'*dphi1x';
D2_bar = dphi2x*k*inv(R22)'*k'*dphi2x';

E1_bar = dphi2x*g*inv(R11)'*g'*dphi1x';
E2_bar = dphi1x*k*inv(R22)'*k'*dphi2x';

m42 = -0.5*F1 - D1_bar*W1/(8*ms1);
m24 = m42';

m53 = -0.5*F3 - D2_bar*W2/(8*ms2);
m35 = m53';

m43 = -dphi1x*g*inv(R11)'*R21*inv(R11)*g'*dphi1x'*W1/(8*ms2) - E1_bar*W2/(4*ms2);
m34 = m43';

m52 = -dphi2x*k*inv(R22)'*R12*inv(R22)'*k'*dphi2x'*W2/(8*ms1) - E2_bar*W1/(4*ms1);
m25 = m52';

m44 = F2 - (D1_bar*W1*m1'+m1*W1'*D1_bar + dphi1x*g*inv(R11)'*R21*inv(R11)*g'*dphi1x'*W2*m2'+...
    +m2*W2'*dphi1x*g*inv(R11)'*R21*inv(R11)'*g'*dphi1x')/8;
m55 = F4 - (D2_bar*W2*m2'+m2*W2'*D2_bar + dphi2x*k*inv(R22)'*R12*inv(R22)*k'*dphi2x'*W1*m1'+...
    +m1*W1'*dphi2x*g*inv(R22)'*R12*inv(R22)'*k'*dphi2x')/8;

M23 = [m24, m25;m34,m35];

M33 = [m44,zeros(num,num);zeros(num,num),m55];

M32 = M23'; 
I2 = eye(2);
D22 = I2 - M23*inv(M33)*M32;
[v22,e22] = eig(D22);
if   sum(diag(e22)>0) == size(D22,1)
    disp('Schur complement for I2 >0 ')
end

D33 = M33 - M32*inv(I2)*M23;
[v33,e33] = eig(D33);
if    sum(diag(e33)>0) == size(D33,1)
    disp('Schur complement for 	M33 >0 ')
end

Y1 = X'*Q1*X + u3'*R11*u3 + d4'*R12*d4;
Y2 = X'*Q2*X + u3'*R21*u3 + d4'*R22*d4;

W1_hat = -a1*m1*(delta3'*W1+Y1);
W2_hat = -a2*m2*(delta4'*W2+Y2);

W3_hat = -a3*((F2*W3-F1*delta3_bar'*W1)-1/4*(dphi1x*g*inv(R11)'*R21*inv(R11)*g'*dphi1x'*W3*m2'*W2+D1_bar*W3*m1'*W1));
W4_hat = -a4*((F4*W4-F3*delta4_bar'*W2)-1/4*(dphi2x*k*inv(R22)'*R12*inv(R22)*k'*dphi2x'*W4*m1'*W1+D2_bar*W4*m2'*W2));

% u3new = u3;
% d4new = d4;
if t <= 450
    u3new =((u3) + exp(-0.07*t)*15*(sin(t)^2*cos(t) + sin(2*t)^2*cos(0.1*t) + sin(-1.2*t)^2*cos(0.5*t) + sin(t)^5 + sin(1.12*t)^2 + cos(2.4*t)*sin(2.4*t)^3));
    d4new =((d4) + exp(-0.07*t)*15*(sin(t)^2*cos(t) + sin(2*t)^2*cos(0.1*t) + sin(-1.2*t)^2*cos(0.5*t) + sin(t)^5 + sin(1.12*t)^2 + cos(2.4*t)*sin(2.4*t)^3));
else
    u3new = u3;
    d4new = d4;
end
u3new
d4new

xout = [f + g*u3new + k*d4new; W1_hat ;W2_hat; W3_hat;W4_hat];
