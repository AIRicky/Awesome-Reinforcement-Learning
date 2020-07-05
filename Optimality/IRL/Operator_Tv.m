% Tv operator
function x_hat =Operator_Tv(x)
n = length(x);
x_hat = zeros(1,n*(n+1)/2);
temp = 1;
for i = 1:n
    for j = i:n
        x_hat(temp) = x(i)*x(j);
        temp = temp +1;
    end
 
end
