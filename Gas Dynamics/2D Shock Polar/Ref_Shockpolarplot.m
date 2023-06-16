function Ref_Shockpolarplot(M,theta)
[weakphi,strongphi]=shockangle(M,theta);
disp(weakphi);
p2byp1M=P2_P1(M,weakphi);
gamma=1.4;
M2=(1/(sin(deg2rad(abs(weakphi-theta)))))*(((gamma-1)*(M*sin(deg2rad(weakphi)))^2+2)/(2*gamma*(M*sin(deg2rad(weakphi)))^2-(gamma-1)))^0.5;
tmax=Theta_max(M2);
l=1000;
angle=linspace(0,tmax,l);
weak=zeros(l,1);
strong=zeros(l,1);
for i=1:l
    [weak(i,1),strong(i,1)]=shockangle(M2,angle(i));
end
weakP2byP1=zeros(l,1);
strongP2byP1=zeros(l,1);
for i=1:l
    weakP2byP1(i,1)=P2_P1(M2,weak(i,1));
    strongP2byP1(i,1)=P2_P1(M2,strong(i,1));
end
hold on;
plot(theta+angle,p2byp1M*weakP2byP1,theta+angle,p2byp1M*strongP2byP1,theta-angle,p2byp1M*weakP2byP1,theta-angle,p2byp1M*strongP2byP1,color='red');
legend(" M2 = "+string(M2));
end