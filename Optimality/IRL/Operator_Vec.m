function VecM = Operator_Vec(M)
[m,n] = size(M);
VecM = reshape(M,m*n,1);