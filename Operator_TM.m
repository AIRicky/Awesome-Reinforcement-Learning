% Tm operator for symmetric matrix
function P_bar =Operator_TM(P)
n = size(P,1);
P_bar = zeros(1,n*(n+1)/2);
temp = 1;
for i = 1:n
    for j = i :n
        if j == i
            P_bar(temp) =  P(i,j);
        else
            P_bar(temp) = 2*P(i,j);
        end
        temp = temp + 1;
    end
   
end
