clc
clear all
close all

 A = [  -1.9440    0.5720    1.4460   -0.5760    0.7360   -0.6010   -0.7220   -0.0880    0.9770    0.3800;
        1.4400    0.3930    1.0230   -0.7110    1.2820   -0.6790    0.0100    0.5880    1.2810   -1.4140;
       -0.8810    1.0580   -1.4920    1.1130   -1.7280    0.4980    0.3130    1.5090   -1.5360   -0.2640;
       -1.1700   -1.0550   -0.0580   -0.7230   -0.9390    1.4530   -1.0870   -0.4860    1.0660    0.2350;
        0.7360   -0.5690    1.4490   -1.3830    0.1160   -0.0520    1.3870    0.6590   -1.6580   -1.4370;
        0.0140    0.6580    0.5860   -0.8500   -0.0740   -1.3350   -0.2610   -1.0210   -0.4490    1.4440;
       -0.7340    0.6210    0.4220   -0.3690   -0.3950   -0.4530    1.2280    0.2130   -1.3800    1.3070;
        0.8200   -1.7460    0.1780   -0.8600   -1.2350   -0.9020    0.3900   -0.6560   -1.6580    1.3290;
        0.8310    0.5690    1.4080    1.5000    1.3960   -0.6050    0.3870   -0.7290    1.7170    1.3090;
        0.0510   -0.2240    1.3940    0.1040   -1.7420   -0.3860   -0.0470   -0.5050   -1.1350    1.3920 ];

B1 = [  -2.0360    1.5600   -0.9070   -1.2140    0.8130    0.0440   -0.7500    0.9010    0.9130    0.0840;
        0.6370    0.4470    1.1540   -1.0910   -0.5750    0.7290    0.6900   -1.8260    0.6350   -0.2090]';

B2 = [ -1.6480    0.1710   -0.3800   -1.4650    1.8540    0.0150    0.4580    0.2550    0.2740   -0.5020;
       -0.7590    1.2560   -1.0760   -0.1010    0.7450    1.7170   -0.0910   -1.3040   -0.7630    1.3450]';

R11 = diag([1,2]);
R12 = diag([3,4]);
R21 = diag([5,6]);
R22 = diag([7,8]);
Q1 = eye(10);
Q2 = eye(10);

% P1 = sdpvar(10,10);
% P2 = sdpvar(10,10);
% 
% Ac = A - 0.5*(B1*inv(R11)*B1'*P1+B2*inv(R22)*B2'*P2);
% 
% cons = [P1*Ac + Ac'*P1 + Q1 + 1/4*(P1*B1*inv(R11)'*R11*inv(R11)*B1'*P1+P2*B2*inv(R22)'*R12*inv(R22)*B2'*P2) == 0;
% P2*Ac + Ac'*P2 + Q2 + 1/4*(P1*B1*inv(R11)'*R21*inv(R11)*B1'*P1+P2*B2*inv(R22)'*R22*inv(R22)*B2'*P2) == 0;]
% 
% solvesdp(cons)
%% Method 1:Lyapunov iterations
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
MaxIter = 10;
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

P1_result = reshape(SP1(iter+1,:),size(A,1),size(A,1))
P2_result = reshape(SP2(iter+1,:),size(A,1),size(A,1))

K1_result = -inv(R11)*B1'*P1_result
K2_result = -inv(R22)*B2'*P2_result

%% Method 2:  Synchronous policy iteration
figure (1)
plot(0:length(Temp1)-1,Temp1,'b-*')
xlabel('迭代次数')
ylabel('P_{1} 迭代误差')
% grid
figure(2)
plot(0:length(Temp2)-1,Temp2,'r-.o')
xlabel('迭代次数')
% grid
ylabel('P_{2} 迭代误差')










