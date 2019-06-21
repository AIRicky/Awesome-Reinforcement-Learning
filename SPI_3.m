clear all
close all
clc

% global t;

P10 = [ 0.4550    0.7570    0.3704;
        0.7570    2.3832    1.2097;
        0.3704    1.2097    0.8877];
P20 = [ 0.8666    1.9790    1.0731
         1.9790    6.7884    3.5192
         1.0731    3.5192    2.1728];
W10 =  Operator_TM(P10);
W20 =  Operator_TM(P20);

x0=[1 0.5 1  W10/2  W20/2  5*rand(1,6) 5*rand(1,6)];
 
%  options = odeset('InitialStep',1,...
%                'MaxStep',500,'RelTol',1e-5,'AbsTol',1e-4,'OutputFcn',@odeplot);
   options = odeset('OutputFcn',@odeplot);
%  [t,x]=ode23(@sys,[tstart tfinal],x0,options);
   [t,x]= ode23('RLdynamic3',[0 1],x0,options);

figure (1);
plot(t,x(:,1:3));
title ('System States');
xlabel ('Time (s)');
figure (2);
plot(t,x(:,4:9));
title ('Parameters of the critic NN of player 1');
xlabel ('Time (s)');
legend ('W1_{c1}','W1_{c2}', 'W1_{c3}','W1_{c4}','W1_{c5}', 'W1_{c6}');
 figure (3);
plot(t,x(:,10:15)); 
title ('Parameters of the critic NN of player 2');
xlabel ('Time (s)');
legend ('W2_{c1}','W2_{c2}', 'W2_{c3}','W2_{c4}','W2_{c5}', 'W2_{c6}');
figure (4);
plot(t,x(:,16:21)); 
title ('Parameters of the Action NN of Player 1');
xlabel ('Time (s)');
legend ('W3_{a1}','W3_{a2}', 'W3_{a3}','W3_{a4}','W3_{a5}', 'W3_{a6}');
figure (5);
plot(t,x(:,22:27)); 
title ('Parameters of the Action NN of Player 2');
xlabel ('Time (s)');
legend ('W4_{a1}','W4_{a2}', 'W4_{a3}','W4_{a4}','W4_{a5}', 'W4_{a6}');

P1_th = [ 0.4550    0.7570    0.3704;
        0.7570    2.3832    1.2097;
        0.3704    1.2097    0.8877]
P1_bar = x(4:9);
P1approx =  Operator_ITM(P1_bar)

P2_th = [ 0.8666    1.9790    1.0731
         1.9790    6.7884    3.5192
         1.0731    3.5192    2.1728]
    P2_bar = x(10:15);
P2approx =  Operator_ITM(P2_bar)
