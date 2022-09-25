f=@(x) sin(x)-x/10;
fsolve(@(x) sin(x)-x/10,1.57)
x=0:0.01:2*pi;
figure(1)
plot(x,sin(x),"Red",x,x/10)
grid on






