clc
clear all
close all

% online solution of two area LFC
q1 = zeros(1,9);
q1(1) = 1;
q1(7) = 1;
q2 = zeros(1,9);
q2(4) = 1;
Q1 = diag(q1);
Q2 = diag(q2);
% R1 = diag([5,0]);
% R2 = diag([0,10]);
R11 = 5;
R12 = 0;
R21 = 0;
R22 = 10;
Tgi = 0.08;
Tti = 0.3;
Kpi = 12;
Tpi = 20;
T12 = 0.2545;
ri = 2.4;
A = zeros(9,9);
A(1,1) = -1/Tpi;
A(1,2) = Kpi/Tpi;
A(1,7) = -Kpi/Tpi;
A(2,2) = -1/Tti;
A(2,3) = 1/Tti;
A(3,1) = -1/(ri*Tgi);
A(3,3) = -1/Tgi;
A(3,8) = 1/Tgi;
A(4,4) =  -1/Tpi;
A(4,5) = Kpi/Tpi;
A(4,7) = Kpi/Tpi;
A(5,5) = -1/Tti;
A(5,6) = 1/Tti;
A(6,4) = -1/(ri*Tgi);
A(6,6) = -1/Tgi;
A(6,9) = 1/Tgi;
A(7,1) = T12;
A(7,4) = -T12;
B = zeros(9,2);
B(8,1) = 1;
B(9,2) = 1;
A
B

% Approach 1
% alpha = 1;
% qq = zeros(1,9);
% qq(1) = 1;
% qq(4) = 1;
% qq(7) = 1;
% Q = diag(qq);
% rr = zeros(1,2);
% rr(1) = 5;
% rr(2) = 10;
% R = diag(rr);
% [K,S,e] = care(A,B,Q,R)


% global t;
% 
% R = 1;
% Q = eye(3);

% P1 = sdpvar(9,9);
% P2 = sdpvar(9,9);
B1 = B(:,1);
B2 = B(:,2);
% Ac = A - 0.5*(B1*inv(R11)*B1'*P1+B2*inv(R22)*B2'*P2);
% cons = [P1*Ac + Ac'*P1 + Q1 + 1/4*(P1*B1*inv(R11)'*R11*inv(R11)*B1'*P1+P2*B2*inv(R22)'*R12*inv(R22)*B2'*P2) == 0;
% P2*Ac + Ac'*P2 + Q2 + 1/4*(P1*B1*inv(R11)'*R21*inv(R11)*B1'*P1+P2*B2*inv(R22)'*R22*inv(R22)*B2'*P2) == 0;]
% solvesdp(cons)
% 
% % cons = [P1*A + A'*P1+Q1-P1*B*inv(R1)*b1'*P1-P1*B*inv(R2)*b2'*P2 - P2'*b2*inv(R2)*b2'*P1 == 0;
% % P2*A + A'*P2+Q2-P2*B*inv(R1)*b1'*P1-P2*B*inv(R2)*b2'*P2 - P1'*b1*inv(R1)*b1'*P2 == 0;]
% % cons = [P1*A + A'*P1+Q1-P1*b1*inv(R1)*b1'*P1-P1*b1*inv(R2)*b2'*P2 - P2'*b2*inv(R2)*b2'*P1 == 0;
% % P2*A + A'*P2+Q2-P2*b2*inv(R1)*b1'*P1-P2*b2*inv(R2)*b2'*P2 - P1'*b1*inv(R1)*b1'*P2 == 0;]
% % solvesdp(cons)

%% Step1: Get initial solution of coupled AREs:
Co = ctrb(A,B1);% Stabilizing
if (rank(Co) == size(A,1))
    disp('Agent 1 is stabilzing!')
end
Ob = obsv(A,sqrtm(Q1)); % Detectable 
if (rank(Ob) == size(A,1))
    disp('Agent 1 is detectable!')
end

Co = ctrb(A,B2);% Stabilizing
if (rank(Co) == size(A,1))
    disp('Agent 2 is stabilzing!')
end
Ob = obsv(A,sqrtm(Q2)); % Detectable 
if (rank(Ob) == size(A,1))
    disp('Agent 2 is detectable!')
end

Z1 = B1*inv(R11)*R21*inv(R11)*B1';
Z2 = B2*inv(R22)*R12*inv(R22)*B2';
S1 = B1*inv(R11)*B1';
S2 = B2*inv(R22)*B2';
[P10] = care(A,B1,Q1,R11);% Pl0
Qaux = Q2+P10*Z1*P10';
[P20] = care(A-B1*inv(R11)*B1'*P10,B2,Qaux,R22);% P20

% Check 
[v1,e1] = eig(A-S1*P10);
real(diag(e1))
[v2,e2] = eig(Q2+P10*Z1*P10);
real(diag(e2))
[v3,e3] = eig(A-S1*P10-S2*P20);
real(diag(e3))

C1 = ctrb(A-S1*P10,B2);
if (rank(C1)==size(A-S1*P10,1))
    disp('Stabilzing!');
end
O2 = obsv(A-S1*P10,sqrtm(Q2+P10*S1*P10));
if (rank(O2)==size(A-S1*P10,1))
    disp('Detectable!')
end


%% Step2: Lyapunov iteration 
MaxIter = 20;
epsion1 = 1e-4;
epsion2 = 1e-4;

SP1 = zeros(MaxIter+1,size(A,1)^2);
SP2 = zeros(MaxIter+1,size(A,1)^2);
SP1(1,:) = reshape(P10,1,size(A,1)^2);
SP2(1,:) = reshape(P20,1,size(A,1)^2);

iter = 0;
temp1 = 10;
temp2 = 10;
Temp1 = zeros(MaxIter+1,1);
Temp2 = zeros(MaxIter+1,1);
Temp1(1) = temp1;
Temp2(1) = temp2;

while (iter <= MaxIter && (temp1 > epsion1 || temp2 > epsion2))
    iter = iter + 1
    P1i = reshape(SP1(iter,:),size(A,1),size(A,1));
    P2i = reshape(SP2(iter,:),size(A,1),size(A,1));
    Anew = A - S1*P1i-S2*P2i;
    Qnew1 = Q1 + P1i* S1*P1i+P2i*Z2*P2i;
    Qnew2 = Q2 + P1i* Z1*P1i+P2i*S2*P2i;
    P1_new = lyap(Anew',Qnew1);%
    P2_new = lyap(Anew',Qnew2);
    SP1(iter+1,:) = reshape(P1_new,1,size(A,1)^2);
    SP2(iter+1,:) = reshape(P2_new,1,size(A,1)^2);
    temp1 = max(abs(SP1(iter+1,:)- SP1(iter,:)));
    temp2 = max(abs(SP2(iter+1,:)- SP2(iter,:)));
    Temp1(iter+1) = temp1;
    Temp2(iter+1) = temp2;
end

temp1
temp2

P1_result = reshape(SP1(iter+1,:),size(A,1),size(A,1));
P2_result = reshape(SP2(iter+1,:),size(A,1),size(A,1));

K1_result = -inv(R11)*B1'*P1_result
K2_result = -inv(R22)*B2'*P2_result



