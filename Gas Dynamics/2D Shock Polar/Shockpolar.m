function [deflection_angle,Shockangle,Pressure_ratio]=Shockpolar(M)
deflection_angle=linspace(0,Theta_max(M),500);
Pressure_ratio=zeros(length(deflection_angle),2);
Shockangle=zeros(length(deflection_angle),2);
for i=1:length(deflection_angle)
    [weak,strong]=shockangle(M,deflection_angle(i));
    Shockangle(i,1)=weak;
    Shockangle(i,2)=strong;
    Pressure_ratio(i,1)=P2_P1(M,weak);
    Pressure_ratio(i,2)=P2_P1(M,strong);
end
shockpolarplot(deflection_angle,Pressure_ratio,M,"b",0);
end



