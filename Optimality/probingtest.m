% Noise test
clc
clear all
close all
Ti = 0.01; % 采样周期
Tstep = 0.001;% 采样周期内积分计算时间间隔
ii = 1:50000;
Magu = 5;
prob = zeros(1,length(ii));
tt = zeros(1,length(ii));
for i = 1: length(ii)
    t = Tstep*i;
    tt(i) = t;
    Magu = 5;
    Damp = exp(-0.05*t);
    prob(i)=Damp*Magu*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)...
        +sin(t)^5+sin(1.12*t)^2+cos(2.4*t));
end
figure
plot(tt,prob)
xlabel('time(s)');
ylabel('Probing signal')