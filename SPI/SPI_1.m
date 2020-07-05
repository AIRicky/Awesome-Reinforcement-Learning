clear all
close all
clc

% global t;

% P10 =  1.1330;
% P20 =  1.1869;

P10 =  1;
P20 =  1;

W10 =  Operator_TM(P10);
W20 =  Operator_TM(P20);

x0=[0 W10 W20 rand rand];
 
%  options = odeset('InitialStep',1,...
%                'MaxStep',500,'RelTol',1e-5,'AbsTol',1e-4,'OutputFcn',@odeplot);
   options = odeset('OutputFcn',@odeplot);
%  [t,x]=ode23(@sys,[tstart tfinal],x0,options);
   [t,x]= ode23('RLdynamic1',[0 500],x0,options);

figure (1);
plot(t,x(:,1));
title ('System States');
xlabel ('Time (s)');
figure (2);
plot(t,x(:,2));
title ('Parameters of the critic NN of player 1');
xlabel ('Time (s)');
legend ('W1_{c1}');
 figure (3);
plot(t,x(:,3)); 
title ('Parameters of the critic NN of player 2');
xlabel ('Time (s)');
legend ('W2_{c1}')
figure (4);
plot(t,x(:,4)); 
title ('Parameters of the Action NN of Player 1');
xlabel ('Time (s)');
legend ('W3_{a1}');
figure (5);
plot(t,x(:,5)); 
title ('Parameters of the Action NN of Player 2');
xlabel ('Time (s)');
legend ('W4_{a1}');

P1_th = 1.1330
P1_bar = x(2);
P1approx =  Operator_ITM(P1_bar)

P2_th = 1.1869
P2_bar = x(3);
P2approx =  Operator_ITM(P2_bar)
