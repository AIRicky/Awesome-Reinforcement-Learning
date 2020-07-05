% ITm operator for symmetric matrix
function P = Operator_ITM(P_bar)
N = length(P_bar);
n = floor(sqrt(2*N));
P = zeros(n,n);
temp = 1;
for i = 1:n
    for j = i :n
        if j == i
          P(i,j) = P_bar(temp);
        else
          P(i,j) = P_bar(temp)/2;
          P(j,i) = P(i,j);
        end
        temp = temp + 1;
    end
end
