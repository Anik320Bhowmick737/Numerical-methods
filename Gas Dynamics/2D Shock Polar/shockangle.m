function [weak,strong]=shockangle(M,Theta)
if Theta<Theta_max(M)
a=@(x) 2*cot(x)*(((M*sin(x)).^2)-1)/(2+M^2*(1.4+cos(2*x)))-tan(deg2rad(Theta));
weak=0;
strong=0;
weak=rad2deg(fsolve(a,asin(1/M)));
strong=rad2deg(fsolve(a,deg2rad(90)));
end
if Theta==Theta_max(M)
a=@(x) 2*cot(x)*(((M*sin(x)).^2)-1)/(2+M^2*(1.4+cos(2*x)))-tan(deg2rad(Theta));    
weak=rad2deg(fsolve(a,1.1));
strong=rad2deg(fsolve(a,1.1));
end
if Theta>Theta_max(M)
    fprintf("Detached");
end
end
