function dx = RLsystem2(t,x,flag,u,w,e1,e2)
A = [0 1;-1 -3];
B1 = [0 0.6]';
B2 = [1 4]';
dx = A*x + B1*(u+e1) + B2*(w+e2);