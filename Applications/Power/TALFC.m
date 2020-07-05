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
R1 = 5;
R2 = 10;
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
alpha = 1;
qq = zeros(1,9);
qq(1) = 1;
qq(4) = 1;
qq(7) = 1;
Q = diag(qq);
rr = zeros(1,2);
rr(1) = 5;
rr(2) = 10;
R = diag(rr);
[K,S,e] = care(A,B,Q,R)


% global t;

R = 1;
Q = eye(3);

P1 = sdpvar(9,9);
P2 = sdpvar(9,9);
b1 = B(:,1);
b2 = B(:,2);
% cons = [P1*A + A'*P1+Q1-P1*B*inv(R1)*b1'*P1-P1*B*inv(R2)*b2'*P2 - P2'*b2*inv(R2)*b2'*P1 == 0;
% P2*A + A'*P2+Q2-P2*B*inv(R1)*b1'*P1-P2*B*inv(R2)*b2'*P2 - P1'*b1*inv(R1)*b1'*P2 == 0;]
cons = [P1*A + A'*P1+Q1-P1*b1*inv(R1)*b1'*P1-P1*b1*inv(R2)*b2'*P2 - P2'*b2*inv(R2)*b2'*P1 == 0;
P2*A + A'*P2+Q2-P2*b2*inv(R1)*b1'*P1-P2*b2*inv(R2)*b2'*P2 - P1'*b1*inv(R1)*b1'*P2 == 0;]
solvesdp(cons)

% g=[1; 3.5;0 ];
% 
% x0=[1 -1 1 ones(1,6) rand(1,6) rand(1,6)];
% % z(1,:)=x0;
% % f=[-0.376*x0(1)+x0(2) ;-1.8115*x0(1)+2*x0(3);0.725*x0(1)]; 
% % phixx=[x0(1)^2 x0(1)*x0(2) x0(1)*x0(3)   x0(2)^2 x0(2)*x0(3) x0(3)^2];
% % dphixx=[2*x0(1) 0 0 ; x0(2) x0(1) 0 ; x0(3) 0  x0(1); 0 2*x0(2) 0;0 x0(3) x0(2); 0 0 2*x0(3) ];
% % uo(1)=-0.5*inv(R)*g'*dphixx'*[x0(10) x0(11) x0(12) x0(13) x0(14) x0(15) ]';
% % V(1)=[x0(4) x0(5) x0(6) x0(7) x0(8) x0(9) ]*phixx';
% % Y(1)=[x0(4) x0(5) x0(6) x0(7) x0(8) x0(9) ]*dphixx*(f+g*uo(1))+([x0(1) x0(2) x0(3) ]*Q*[x0(1) x0(2) x0(3) ]'+uo(1)*R*uo(1)')  ;
% 
% % e2(1)= p'*([x0(6) x0(7) x0(8)] -[x0(3) x0(4) x0(5)])';   
% % for policy=2:10
% % policy 
%  
% %  
% %  options = odeset('InitialStep',1,...
% %                'MaxStep',500,'RelTol',1e-5,'AbsTol',1e-4,'OutputFcn',@odeplot);
%    options = odeset('OutputFcn',@odeplot);
% %  [t,x]=ode23(@sys,[tstart tfinal],x0,options);
%    [t,x]= ode23('dynamicsnew3',[0 600],x0,options);
% %    x01=[1 1 -2 x(length(x),1:12)];
% %    [t,x]= ode23('dynamicsnew3',[0 500],x01,options);
% %     
% %    H=[x(length(t),6) x(length(t),7) x(length(t),8)]';
% % x0=[x(length(t),1) x(length(t),2) x(length(t),3)  x(length(t),4) x(length(t),5) x(length(t),6) x(length(t),7)  x(length(t),8) x(length(t),9) x(length(t),10) x(length(t),11) x(length(t),12) x(length(t),13) x(length(t),14) x(length(t),15) ]; %reset and initializing of the weights
% % dphixo=[2*x(length(t),1) 0 0; x(length(t),2) x(length(t),1)  0; x(length(t),3) 0  x(length(t),1);0 2*x(length(t),2) 0;0 x(length(t),3) x(length(t),2) ;0 0 2*x(length(t),3) ];
% % 
% %  uo(policy)=-0.5*inv(R)*g'*dphixo'*[x(length(t),10) x(length(t),11) x(length(t),12) x(length(t),13) x(length(t),14) x(length(t),15)  ]';
% %  f=[-0.376*x(length(t),1)+x(length(t),2) ;-1.8115*x(length(t),1)+2*x(length(t),3);0.725*x(length(t),1)]; 
% % %  V(policy)=[x(length(t),3) x(length(t),4) x(length(t),5)]*dphixo*(f+g*uo(policy));
% % phixo=[x(length(t),1)^2 x(length(t),1)*x(length(t),2) x(length(t),1)*x(length(t),3)  x(length(t),2)^2 x(length(t),2)*x(length(t),3) x(length(t),3)^2 ]';     
% % V(policy)=[x(length(t),4) x(length(t),5) x(length(t),6) x(length(t),7) x(length(t),8) x(length(t),9)  ]*phixo;
% % 
% % z(policy,:)=x(length(t),:);
% % Y(policy)=[x(length(t),4) x(length(t),5) x(length(t),6) x(length(t),7) x(length(t),8) x(length(t),9)]*dphixo*(f+g*uo(policy))+([x(length(t),1) x(length(t),2) x(length(t),3) ]*Q*[x(length(t),1) x(length(t),2) x(length(t),3) ]'+uo(policy)*R*uo(policy)')  ;
% % 
% % 
% % 
% % p=(-0.5*inv(R)*g'*dphixo')';
% % 
% % e2(policy)= p'*([x(length(t),10) x(length(t),11) x(length(t),12) x(length(t),13) x(length(t),14) x(length(t),15) ] -[x(length(t),4) x(length(t),5) x(length(t),6) x(length(t),7) x(length(t),8) x(length(t),9) ])';
% % 
% % end
% % 
% 
% % figure (1);
% % plot(t,x(:,1:3));
% % figure (2);
% % plot(t,x(:,4:9));
% %%%%%%%%%%%%%%%%%5
% figure (1);
% plot(t,x(:,1:3));
% title ('System States');
% xlabel ('Time (s)');
% figure (2);
% plot(t,x(:,4:9));
% title ('Parameters of the critic NN');
% xlabel ('Time (s)');
% legend ('W_{c1}','W_{c2}', 'W_{c3}','W_{c4}','W_{c5}', 'W_{c6}');
%  figure (3);
% plot(t,x(:,10:15)); 
% title ('Parameters of the actor NN');
% xlabel ('Time (s)');
% legend ('W_{a1}','W_{a2}', 'W_{a3}','W_{a4}','W_{a5}', 'W_{a6}');
%  figure (4);
% plot(t,x(:,16:21)); 
% title ('Parameters of the disturbance NN');
% xlabel ('Time (s)');
% legend ('W_{d1}','W_{d2}', 'W_{d3}','W_{d4}','W_{d5}', 'W_{d6}');
% %%%%%%%%%%%%%%%%%%
% 
% 
% 
% B=[0 0 1]';
% A=[-1.01887 0.90506 -0.00215; 0.82225 -1.07741 -0.17555; 0 0 -1];
% D=[1 0 0]';
% % 
% % A=[-0.376 1 0;-1.8115 0 2;0.725 0 0];
% % B=[1;3.5;0];
% Q = eye(3);
% 
% g = 5;
% B_new = [D , B];
% m1 = size(D,2);
% m2 = size(B,2);
% R = [-g^2*eye(m1) zeros(m1,m2) ; zeros(m2,m1) eye(m2)];
% 
% PTheor= care(A,B_new,Q,R);
% 
% Papprox=[x(length(x),4)    x(length(x),5)/2  x(length(x),6)/2
%      x(length(x),5)/2  x(length(x),7) x(length(x),8)/2  
%      x(length(x),6)/2 x(length(x),8)/2 x(length(x),9) 
%      ];







