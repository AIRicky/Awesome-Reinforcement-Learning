clc
clear all
close all

disp('-----在线无模型求解线性连续时间系统微分博弈均衡------ ')
% System dynamics
A = [0 1;-1 -3];
B1 = [0 0.6]';
B2 = [1 4]';
Q = diag(ones(1,2));
R = 1;
gammar = 11;

% 矩阵维数
n = size(A,1);
m = size(B1,2);
q = size(B2,2);
num = n*(n+1)/2;

% 理论值
B = [B2 , B1];
m1 = size(B2,2);
m2 = size(B1,2);
RX = [-gammar^2*eye(m1) zeros(m1,m2) ; zeros(m2,m1) eye(m2)];
P_star = care(A,B,Q,RX);
K_star = inv(R)*B1'*P_star;
L_star = 1/gammar^2*B2'*P_star;

%% IRL学习求解
Ti = 0.02; % 采样周期
Tstep = 0.001;% 采样周期内积分计算时间间隔
M = Ti/Tstep; % 用于获取数值积分所需的数据 M个

N = 25; 
Nmin = n*(n+1)/2+m*n+n*q; % 最小采样数据
if N > Nmin
    disp('    采样数据个数满足要求');
end

stopflag = 1;  % 学习flag
epsion = 1e-2; % 学习误差上界
MaxIter = 20;  % 最大学习迭代次数

In = diag(ones(1,n));
P0 = 0;

x0 = [0.2,0.2]';
%  x0 = [0.01,0.01]';

K = zeros(MaxIter,n); 
K(1,:) = [0,0];
L = zeros(MaxIter,n);
L(1,:) = [0,0];
Pbar_Hist = zeros(MaxIter-1,num);

% 初始控制稳定性检查
[v,e]= eig(A-B1*K(1,:));
if sum(diag(real(e))< 0) == n
    disp('    初始控制有效');
end

MaxCount = N*MaxIter;
% Tf = MaxCount*Ti;

% Define history matrix
X_hist = zeros(n,MaxCount);
U_hist = zeros(1,MaxCount);
W_hist = zeros(1,MaxCount);
E1_hist = zeros(1,MaxCount);
E2_hist = zeros(1,MaxCount);
% K_hist = zeros(MaxIter,n);
% L_hist = zeros(MaxIter,n);

iter = 0;
newSt = x0;

tcount = 1;

while(stopflag) % first loop: Iteration loop
   
    iter = iter +1;

    Phi = zeros(num+2*n,N); 
    Theta = zeros(N,1);
    Result = zeros(num+2*n,1);
    
    % second loop: collect N data sets;
    for jj = 1:N    
      
        u_Int = zeros(1,M);
        w_Int = zeros(1,M);
        x_Int = zeros(n,M);
        r_Int = zeros(1,M);
        phi_Int = zeros(num+2*n,M);
        theta_Int = zeros(1,M);
%         e1_Int = (1e-2)*normrnd(0,0.1,1,M); % Gauss white noise
%         e2_Int = (1e-2)*normrnd(0,0.05,1,M); % Gauss white noise  
        e1_Int = zeros(1,M);
        e2_Int = zeros(1,M);
 
        X_hist(:,tcount) = newSt;
        U_hist(tcount) = -K(iter,:)*newSt;
        W_hist(tcount) = L(iter,:)*newSt;
%         E1_hist(tcount) = e1_Int(1);
%         E2_hist(tcount) = e2_Int(1); 
        tcount = tcount + 1;
        
        % begin third loop: calculate the integral term
        for ii = 1:M
            x_Int(:,ii) = newSt;
            u_Int(ii) = -K(iter,:)*newSt;
            w_Int(ii) = L(iter,:)*newSt;
%             e1_Int(ii) = 0.01*sin(Ti*tcount +Tstep*ii);
%             e2_Int(ii) = 0.005*sin(10*(Ti*tcount +Tstep*ii));
                t = Ti*tcount +Tstep*ii;
%                 e1_Int(ii) =0.05*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t));
%                 e2_Int(ii) =0.005*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t));
        
%             e1_Int(ii) = exp(-10*t)*0.25*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t));
%             e2_Int(ii) = exp(-10*t)*0.25*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t));
            if  iter == 1
%                 e1_Int(ii) = 0.01*sin(Ti*tcount +Tstep*ii);
%                 e2_Int(ii) = 0.005*sin(10*(Ti*tcount +Tstep*ii));
                  e1_Int(ii) = 0;
                  e2_Int(ii) = 0;
            else
                Damp = exp(-10*t);
%                   Damp = 1;
                Magu = 5;
                Magw = 5;
                e1_Int(ii) = Damp*Magu*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t));
                e2_Int(ii) = Damp*Magw*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t));

            end
%             e1_Int(ii) = exp(-10*t)*0.25*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t));
%             e2_Int(ii) = exp(-10*t)*0.25*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t));
            
            r_Int(ii) = x_Int(:,ii)'*Q*x_Int(:,ii)+ u_Int(ii)'*R*u_Int(ii)-gammar^2*w_Int(ii)'*w_Int(ii);% theta_int
            phi_Int(num+1:num+n,ii) = 2*(kron(x_Int(:,ii),e1_Int(ii)))'*(kron(In,R));
            phi_Int(num+n+1:num+2*n,ii) = 2*gammar^2*(kron(x_Int(:,ii),e2_Int(ii)))';
          
            u = u_Int(ii);
            w = w_Int(ii);
            e1 = e1_Int(ii);
            e2 = e2_Int(ii);
            temp = newSt;
            [T,Xf] = ode45('RLsystem2',[0 Tstep],newSt,[],u,w,e1,e2);
%             [T,Xf] = ode45('RLsystem',[0 Tstep],newSt,[],u,w,0,0);
            a = size(Xf);
            newSt = Xf(a(1),:)';
        end
        E1_hist(tcount) = e1_Int(1);
        E2_hist(tcount) = e2_Int(1); 
        % end third loop
        Theta(jj,1) =  RLIntegral(r_Int,Tstep);
        Phi(1:num,jj) = Operator_Tv(x_Int(:,1)) - Operator_Tv(x_Int(:,M));% x_Int(M) maybe newSt
        temp1 = zeros(n,1);
        temp2 = zeros(n,1);
        for lii = 1:n
            temp1(lii) = RLIntegral(phi_Int(num+lii,:),Tstep);
            temp2(lii) = RLIntegral(phi_Int(num+n+lii,:),Tstep);
        end
        Phi(num+1:num+n,jj) = temp1;
        Phi(num+n+1:num+2*n,jj) = temp2;
    end
    
    % Update K and L
    Phi
%     Result = inv(Phi*Phi')*Phi*Theta;
%     Result = pinv(Phi*Phi')*Phi*Theta;
    Result = pinv(Phi')*Theta; 
    Pbar(iter,:) = Result(1:num);
    K(iter+1,:) = Result(num+1:num+n);
    L(iter+1,:) = Result(num+n+1:num+2*n);
    
    % 迭代过程系统稳定性检查
    iter 
    [v,e]= eig(A-B1*K(iter+1,:))
    if sum(diag(real(e))< 0) == n
        ['第',num2str(iter) '次迭代生成的控制策略能镇定系统']
    else
        ['第',num2str(iter) '次迭代生成的控制策略不能镇定系统']
    end
    
    if  iter == 1
%         StopIter = norm(Pbar(iter,:)-P0,2)
          StopIter = sum(abs(Pbar(iter,:)-P0))
    else
%         StopIter = norm(Pbar(iter,:)-Pbar(iter-1,:),2)
          StopIter = sum(abs(Pbar(iter,:)-Pbar(iter-1,:)))
    end

    if (iter+1> MaxIter) || (StopIter <= epsion)
        stopflag = 0;
    end
end

% 将学习到的策略应用于实际系统 仿真
MM = 50;
Xnew = zeros(MM*N,n);
for tj = 1:MM*N
    xn0 = newSt;
    % % [t,x]=ode23(@sys,[tstart tfinal],x0,options);
    e1 = 0;
    e2 = 0;
    u = K(iter+1,:)*newSt;
    w = L(iter+1,:)*newSt;
    [T,Xf] = ode45('RLsystem2',[0 Ti],newSt,[],u,w,e1,e2);
    a = size(Xf);
    newSt = Xf(a(1),:)';
    Xnew(tj,:) = newSt;
end

% figure(5)
NData = [temp';Xnew]';
% ttt = (tcount-1):(tcount-1+MM*N);
% plot(ttt*Ti,NData(1,:),ttt*Ti,NData(2,:))
% xlabel('Time(s)')
% ylabel('System trajectory')
% legend('x1','x2')

iter
P_learning = Operator_ITM(Pbar(end,:))
P_star
K_learning1 = inv(R)*B1'*P_learning
K
K_learning = K(iter+1,:)
K_star 
L_learning1 = 1/gammar^2*B2'*P_learning
L
L_Learning = L(iter+1,:)
L_star

% figure(1)
% plot(Pbar)
% xlabel('迭代次数')

figure(1)
Data = [X_hist(:,1:tcount-1),NData];
plot((0:tcount-1+MM*N)*Ti,Data(1,:),(0:tcount-1+MM*N)*Ti,Data(2,:))
% % ylim([-1,1])
xlabel('Time(s)')
ylabel('System trajectory')
legend('x1','x2')

figure(2)
% tt = 0:MaxIter-1;
tk = 0:(iter);
NPbar = [zeros(1,num);Pbar];
plot(Ti*N*tk,NPbar(:,1),'b-*')
hold on
plot(Ti*N*tk,NPbar(:,2),'r-o')
plot(Ti*N*tk,NPbar(:,3),'g-+')
% plot([zeros(1,num);Pbar],'-*')
% xbase = 1:size([zeros(1,num);Pbar],1);
% set(gca,'XTickLabel',xbase*Ti*N);

xlabel('time(s)')
ylabel('convergence of critic network')

figure(3)
% tp = 1: MaxIter;
% plot(tt*N*Ti,K)
tk = 0:(iter);
plot(Ti*N*tk,K(1:length(tk),1),'b-*')
hold on
plot(Ti*N*tk,K(1:length(tk),2),'r-o')
% plot(Ti*N*tk,K(:,3),'g-+')
% plot(K,'-*')
% set(gca,'XTickLabel',(1:size(K,1)-1)*Ti*N);
xlabel('time(s)')
ylabel('convergence of action network')

figure(4)
% tp = 1: MaxIter;
% plot(tt*N*Ti,L)
% plot(L,'-*')
tk = 0:(iter);
plot(Ti*N*tk,L(1:length(tk),1),'b-*')
hold on
plot(Ti*N*tk,L(1:length(tk),2),'r-o')
% set(gca,'XTickLabel',(1:size(L,1)-1)*Ti*N);
xlabel('time(s)')
ylabel('convergence of disturbance network')

% figure(5)
% % tp = 1: MaxIter;
% % plot(tt*N*Ti,L)
% % plot(L,'-*')
% tk = 0:(iter);
% plot(Ti*N*tk,U_hist(1:length(tk)),'b-*')
% hold on
% plot(Ti*N*tk,W_hist(1:length(tk)),'r-o')
% % set(gca,'XTickLabel',(1:size(L,1)-1)*Ti*N);
% xlabel('time(s)')
% ylabel('convergence of disturbance network')
% legend('U','W')

