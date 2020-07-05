
function xout=dynamicsnew3(t,x)
global R;
global g;
global dphix;
global H;
global W;
% global t;
a = 20;
b = 1;
c = 1;
R = 1;
gamma = 5;
Q=[1 0 0 ; 0 1 0 ; 0 0 1];
% t=0;
x1 = x(1);
x2 = x(2);
x3 = x(3);
W = [x(4) x(5) x(6) x(7) x(8) x(9) ]';
H = [x(10) x(11) x(12) x(13) x(14) x(15) ]';
WD = [x(16) x(17) x(18) x(19) x(20) x(21) ]';
% u=[x(9)];


phix = [x1^2 x1*x2 x1*x3 x2^2 x2*x3 x3^2]';
dphix = [2*x1 0 0; x2 x1 0; x3 0 x1; 0 2*x2 0; 0 x3 x2; 0 0 2*x3];
f = [-1.01887*x1+0.90506*x2-0.00215*x3;0.82225*x1-1.07741*x2-0.17555*x3;-1*x3];
g = [0; 0;1];
k = [1; 0;0];
u = -0.5*inv(R)*g'*dphix'*H;
d = 0.5*(1/(gamma^2))*k'*dphix'*WD;
s = dphix*f + dphix*g*u +dphix*k*d;
Y =(-[x1 x2 x3]*Q*[x1 x2 x3]'-u*R*u'+(gamma^2)*d*d');
e = W'*s-Y;

%     W(:,k+1)=W(:,k)-e(k)*a*(s(:,k)/(s(:,k)'*s(:,k)+1));
% W(:,k+1)=W(:,k)-a*s(:,k)*(e(k));
%


p = (-0.5*inv(R)*g'*dphix')';

e2 = p'*(H -W);

E2 = H -W;
Win =-a*(s./(s'*s+1)^2)*e';
%      Hin=-b*tanh(E2);
%     Hin=0.25*dphix*g*inv(R)*g'*dphix'*(W-H);
%     Hin=0.25*b*inv(R)*g'*dphix'*(W-H)*dphix*g*inv(R)'; 
F = 10;%Positive definite matrix  
%  D1=dphix*g*inv(R)*g'*dphix';
% Hin=0.5*b*D1*(W-H)+0.25*b*D1*H*(s'./(s'*s+1)^2)*W;
%  F1=10*ones(1,length(H));

F1 = 1*eye(length(W));

sbar = (s/(s'*s+1));
   
D1 = dphix*g*inv(R)*g'*dphix';
E1 =(1/(gamma^2))*dphix*k*k'*dphix';
F2 = 10*eye(length(H));
 
Hin = b*(F2*W-F2*H) + 0.25*b*D1*H*(s'./(s'*s+1)^2)*W;
WDin = c*(F2*W-F2*WD)-0.25*c*E1*WD*(s'./(s'*s+1)^2)*W;
if t <= 550
    unew =((u) + exp(-0.007*t)*37*(sin(t)^2*cos(t) + sin(2*t)^2*cos(0.1*t) + sin(-1.2*t)^2*cos(0.5*t) + sin(t)^5 + sin(1.12*t)^2 + cos(2.4*t)*sin(2.4*t)^3));
    dnew =((u) + exp(-0.007*t)*37*(sin(t)^2*cos(t) + sin(2*t)^2*cos(0.1*t) + sin(-1.2*t)^2*cos(0.5*t) + sin(t)^5 + sin(1.12*t)^2 + cos(2.4*t)*sin(2.4*t)^3));
    % unew=((u)+2*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*
    % t)+sin(t)^5));
else
    unew = u;
    dnew = u;
end

xout = [f + g*unew + k*dnew; Win ; Hin; WDin];
