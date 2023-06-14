function Shockpolarplot(M)
tmax=Theta_max(M);
l=1000;
angle=linspace(0,tmax,l);
weak=zeros(l,1);
strong=zeros(l,1);
for i=1:l
    [weak(i,1),strong(i,1)]=shockangle(M,angle(i));
end
weakP2byP1=zeros(l,1);
strongP2byP1=zeros(l,1);
for i=1:l
    weakP2byP1(i,1)=P2_P1(M,weak(i,1));
    strongP2byP1(i,1)=P2_P1(M,strong(i,1));
end
hold on;
plot(angle,weakP2byP1,angle,strongP2byP1,-angle,weakP2byP1,-angle,strongP2byP1,color='blue');
end