function y = RLIntegral(x,th)
% x:　离散点的数目
% th:  离散点的时间间隔
lm = length(x);
y = 0;
for lj = 1:lm-1
    y = y + (x(lj)+x(lj+1))*th/2;
end
