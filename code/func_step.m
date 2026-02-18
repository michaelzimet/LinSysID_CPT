function dy = func_lin(t,y)
global A B;

x1 = y(1);
x2 = y(2);
u = 0;

dx1 = A(1,1)*x1 + A(1,2)*x2 + B(1)*u;
dx2 = A(2,1)*x1 + A(2,2)*x2 +  B(2)*u;

dy = [dx1;dx2];
end