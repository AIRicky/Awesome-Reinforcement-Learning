function y = RLIntegral(x,th)
% x:����ɢ�����Ŀ
% th:  ��ɢ���ʱ����
lm = length(x);
y = 0;
for lj = 1:lm-1
    y = y + (x(lj)+x(lj+1))*th/2;
end
