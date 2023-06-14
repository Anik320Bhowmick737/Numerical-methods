function theta_max=Theta_max(M)
c=2.4*(2.4*(M/2)^4+0.2*M^2+1);
c=c^0.5;
c=c-1+2.4*(M/2)^2;
c=c/(1.4*M^2);
a=c;
theta_max=atand((2*cot(asin(a^0.5))*((M^2)*a-1))/(2+(M^2)*(2.4-2*a)));
end
