function M=shockmach(gamma,P02,P1)
a=@(x) (P02/P1)-(((0.5*(gamma+1)*x^2)^(gamma/(gamma-1)))/((2*gamma*x^2)/(gamma+1)-(gamma-1)/(gamma+1))^(1/(gamma-1)));
M=fsolve(a,1.5);
end
