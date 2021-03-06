function dx = RLsystem(t,x,flag,u,w,e1,e2)
A = [-0.0665 8 0 0;
    0 -3.663 3.663 0;
   -6.86 0 -13.736 -13.736;
   0.6 0 0 0];
B1 = [0 0 13.736 0]';
B2 = [-8 0 0 0]';
dx = A*x + B1*(u+e1) + B2*(w+e2);