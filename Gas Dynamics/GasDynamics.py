import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
## Isentropic
class Isentropic_relations():
    def __init__(isen,gamma,Mach=False,TT0=False,PP0=False,rhorho0=False,AAstar_sup=False,AAstar_sub=False,TTstar=False,PPstar=False,rhorhostar=False):
        isen.gamma=gamma
        isen.Mach=Mach
        isen.TT0=TT0
        isen.PP0=PP0
        isen.rhorho0=rhorho0
        isen.AAstar_sup=AAstar_sup
        isen.AAstar_sub=AAstar_sub
        isen.TTstar=TTstar
        isen.PPstar=PPstar
        isen.rhorhostar=rhorhostar
    def G(isen,M=False):
        return (1+0.5*(isen.gamma-1)*(M)**2)       
    def fromMach(isen):
        if isen.Mach==False:
            print("Give Valid input")
        else:
            g=isen.G(M=isen.Mach)  
            g0=isen.G(M=1)  
            AbyAstar=(1/isen.Mach)*((2/(isen.gamma+1))*g)**(0.5*(isen.gamma+1)/(isen.gamma-1))
            TbyT0=1/g
            PbyP0=(TbyT0)**((isen.gamma)/(isen.gamma-1))
            TbyTstar=g0/g
            PbyPstar=TbyTstar**((isen.gamma)/(isen.gamma-1))
            rhobyrho0=PbyP0/TbyT0
            rhobyrhostar=PbyPstar/TbyTstar
            if isen.Mach>=1:
                Mach_angle=np.arcsin(1/isen.Mach)
                PM_angle=(((isen.gamma+1)/(isen.gamma-1))**0.5)*(np.arctan((((isen.gamma-1)/(isen.gamma+1))*((isen.Mach)**2-1))**0.5))-np.arctan(((isen.Mach)**2-1)**0.5)
                Data={"Mach":isen.Mach,"A/A*":AbyAstar,"T/T0":TbyT0,"P/P0":PbyP0,"rho/rho0":rhobyrho0,"T/T*":TbyTstar,"P/P*":PbyPstar,"rho/rho*":rhobyrhostar,"Mach_angle":Mach_angle,"PM_angle":PM_angle}
                Listed_Data=list(Data.values())   
            else:
                Data={"Mach":isen.Mach,"A/A*":AbyAstar,"T/T0":TbyT0,"P/P0":PbyP0,"rho/rho0":rhobyrho0,"T/T*":TbyTstar,"P/P*":PbyPstar,"rho/rho*":rhobyrhostar,"Mach_angle":None,"PM_angle":None} 
                Listed_Data=list(Data.values())   
            return Data,Listed_Data
    def fromAAstar_sup(isen):
        if isen.AAstar_sup<=1:
            print('A/A* must be greater than 1')
        else:
            def func(x):
                f=(1/x)*((2/(isen.gamma+1))*(1+0.5*(isen.gamma-1)*x**2))**(0.5*(isen.gamma+1)/(isen.gamma-1))-isen.AAstar_sup
                return f
            isen.Mach=fsolve(func,1.2)[0]
            Data=isen.fromMach()
            return Data
    def fromAAstar_sub(isen):
        if isen.AA_star_sub<=1:
            print('A/A* must be greater than 1')
        else:    
            def func(x):
                f=(1/x)*((2/(isen.gamma+1))*(1+0.5*(isen.gamma-1)*x**2))**(0.5*(isen.gamma+1)/(isen.gamma-1))-isen.AAstar_sub
                return f
            isen.Mach=fsolve(func,0.3)[0]
            Data=isen.fromMach()
            return Data    
    def fromTT0(isen):
        if isen.TT0==False or isen.TT0>=1:
            print("Give Valid input")
        else:    
            isen.Mach=((1/isen.TT0-1)*(2/(isen.gamma-1)))**0.5
            Data=isen.fromMach()
            return Data
    def fromPP0(isen):
        if isen.PP0==False or isen.TT0>=1:
            print("Give Valid input")
        else:    
            isen.TT0=(isen.PP0)**((isen.gamma-1)/(isen.gamma))
            Data = isen.fromTT0()
            return Data
    def fromrhorho0(isen):
        if isen.rhorho0==False:
            print("Give Valid input")
        else:    
            isen.PP0=(isen.rhorho0)**(isen.gamma)  
            Data=isen.fromPP0()   
            return Data   
    def fromTTstar(isen):
        if isen.TTstar==False:
            print("Enter valid number")
        else:  
            isen.Mach=(((((1+isen.gamma)/2)*(1/(isen.TTstar)))-1)*(2/(isen.gamma-1)))**0.5
            Data=isen.fromMach()
            return Data
    def fromPPstar(isen):
        if isen.PPstar==False:
            print("Enter valid number")
        else:
            isen.TTstar=(isen.PPstar)**((isen.gamma-1)/(isen.gamma))
            Data = isen.fromTTstar()
            return Data
    def fromrhorhostar(isen):
        if isen.rhorhostar==False:
            print("Enter valid number")
        else:
            isen.PPstar=(isen.rhorhostar)**(isen.gamma)
            Data = isen.fromPPstar()
            return Data 
##Normal Shock                  
class NormalShock_relations():
    def __init__(Norm,gamma,M1=False,M2=False,P02byP01=False,P1byP02=False,P2byP1=False,rho2byrho1=False,T2byT1=False):
        Norm.gamma=gamma
        Norm.M1=M1
        Norm.M2=M2
        Norm.P02byP01=P02byP01
        Norm.P1byP02=P1byP02
        Norm.P2byP1=P2byP1
        Norm.rho2byrho1=rho2byrho1
        Norm.T2byT1=T2byT1
    def fromM1(Norm):
        if Norm.M1 == False or Norm.M1 <=1:
            print("Give Valid Input")
        else:
            Norm.M2 = (((Norm.M1)**2 + 2/(Norm.gamma-1))/(((2*Norm.gamma)/(Norm.gamma-1))*(Norm.M1)**2-1))**0.5
            Norm.P2byP1 = (1+Norm.gamma*(Norm.M1**2))/(1+Norm.gamma*(Norm.M2)**2)
            Norm.P02byP01 = (Isentropic_relations(Norm.gamma,Mach=Norm.M1).fromMach()[0]['P/P0'])/(Isentropic_relations(Norm.gamma,Mach=Norm.M2).fromMach()[0]['P/P0'])*Norm.P2byP1 
            Norm.P1byP02 = (Isentropic_relations(Norm.gamma,Mach=Norm.M1).fromMach()[0]['P/P0'])/(Norm.P02byP01)
            Norm.T2byT1 = (1+0.5*(Norm.gamma-1)*(Norm.M1**2))/(1+0.5*(Norm.gamma-1)*(Norm.M2)**2)
            Norm.rho2byrho1 = (Norm.P2byP1/Norm.T2byT1)
            Data = {"M1":Norm.M1,"M2":Norm.M2,"P02/P01":Norm.P02byP01,"P1/P02":Norm.P1byP02,"P2/P1":Norm.P2byP1,"rho2/rho1":Norm.rho2byrho1,"T2/T1": Norm.T2byT1}
            return Data
    def fromM2(Norm):
        if Norm.M2 == False or Norm.M2 >=1:
            print("Give Valid Input from M2")
        else:
            Norm.M1=(((Norm.M2)**2 + 2/(Norm.gamma-1))/(((2*Norm.gamma)/(Norm.gamma-1))*(Norm.M2)**2-1))**0.5
            Data=Norm.fromM1()
            return Data
    def fromP02byP01(Norm):
        if Norm.P02byP01==False or Norm.P02byP01>=1:
            print("Give Valid Input from P02byP01")
        else:
            def func(x):
                t=(Norm.gamma/(Norm.gamma-1))
                f=((((Norm.gamma+1)*(x**2))/(2+(Norm.gamma-1)*x**2))**t)*((Norm.gamma+1)/(2*Norm.gamma*(x**2)-(Norm.gamma-1)))**(1/(Norm.gamma-1))-Norm.P02byP01
                return f
            Norm.M1 = fsolve(func,2)[0]  
            Data=Norm.fromM1()
            return Data
    def fromP1byP02(Norm):
        if Norm.P1byP02==False or Norm.P1byP02>=0.52828178:
            print("Give Valid Input from P1byP02")
        else:
            def func(x):
                f=((((2*Norm.gamma)/(Norm.gamma+1))*(x**2)-((Norm.gamma-1)/(Norm.gamma+1)))**(1/(Norm.gamma-1)))*(0.5*(Norm.gamma+1)*(x**2))**(Norm.gamma/(1-Norm.gamma))-Norm.P1byP02
                return f
            Norm.M1 = fsolve(func,1.5)[0]  
            Data=Norm.fromM1()
            return Data
    def fromP2byP1(Norm):
        if Norm.P2byP1==False or Norm.P2byP1<=1:
            print("Give Valid Input from P2byP1") 
        else:
            def func(x):
                y = (((x)**2 + 2/(Norm.gamma-1))/(((2*Norm.gamma)/(Norm.gamma-1))*(x)**2-1))   
                f=(1+Norm.gamma*(x**2))/(1+Norm.gamma*y)-Norm.P2byP1 
                return f 
            Norm.M1 = fsolve(func,5)[0]
            Data=Norm.fromM1()
            return Data       
    def fromT2byT1(Norm):
        if Norm.T2byT1==False or Norm.T2byT1<=1:
            print("Give Valid Input from T2byT1") 
        else:
            def func(x):
                y = (((x)**2 + 2/(Norm.gamma-1))/(((2*Norm.gamma)/(Norm.gamma-1))*(x)**2-1))   
                f=(1+0.5*(Norm.gamma-1)*(x**2))/(1+0.5*(Norm.gamma-1)*y)-Norm.T2byT1
                return f 
            Norm.M1 = fsolve(func,5)[0]
            Data=Norm.fromM1()
            return Data    
    def fromrho2byrho1(Norm):
        if Norm.rho2byrho1==False or Norm.rho2byrho1<=1 or Norm.rho2byrho1>=6:
            print("Give Valid Input from rho2byrho1") 
        else:
            def func(x):
                f=(((Norm.gamma+1)*(x**2))/((Norm.gamma-1)*(x**2)+2))-Norm.rho2byrho1
                return f 
            Norm.M1 = fsolve(func,3)[0]
            Data=Norm.fromM1()
            return Data            
##Oblique Shock
class ObliqueShock_relations():
    def __init__(Obl,gamma,M1=False,turn_angle=False,wave_angle=False,M1n=False):
        Obl.gamma=gamma
        Obl.M1=M1
        Obl.turn_angle=np.deg2rad(turn_angle)
        Obl.wave_angle_at_max_deflection=np.arcsin(((1/(Obl.gamma*(Obl.M1)**2))*(0.25*(Obl.gamma+1)*(Obl.M1)**2-1+((Obl.gamma+1)*((Obl.gamma+1)*(0.5*Obl.M1)**4+(0.5*(Obl.gamma-1))*(Obl.M1)**2+1))**0.5))**0.5)
        Obl.max_delfection_angle=np.arctan((1/(np.tan(Obl.wave_angle_at_max_deflection)))*(((((Obl.M1)*(np.sin(Obl.wave_angle_at_max_deflection)))**2)-1)/(0.5*(Obl.gamma+1)*(Obl.M1)**2-((((Obl.M1)*(np.sin(Obl.wave_angle_at_max_deflection)))**2)-1))))
        Obl.wave_angle=np.deg2rad(wave_angle)
        Obl.M1n=M1n
    def fromM1_turn_angle(Obl):
        def shock_angle(x):
            t=np.tan(x)
            f= (1/t)*(((((Obl.M1)*(np.sin(x)))**2)-1)/(0.5*(Obl.gamma+1)*(Obl.M1)**2-((((Obl.M1)*(np.sin(x)))**2)-1))) - np.tan(Obl.turn_angle)
            return f
        Obl.theta_w=fsolve(shock_angle,0.05)[0]
        Obl.theta_s=fsolve(shock_angle,0.5*np.pi-0.2)[0]
        Obl.Mn1w=Obl.M1*np.sin(Obl.theta_w)
        Obl.Mn1s=Obl.M1*np.sin(Obl.theta_s)
        w=NormalShock_relations(Obl.gamma,M1=Obl.Mn1w).fromM1()
        M2w=w['M2']/(np.sin(Obl.theta_w-Obl.turn_angle))
        s=NormalShock_relations(Obl.gamma,M1=Obl.Mn1s).fromM1()
        M2s=s['M2']/(np.sin(Obl.theta_s-Obl.turn_angle))
        Data_w={"M1":Obl.M1,"Weak Shock Angle":np.rad2deg(Obl.theta_w),"P2/P1":w['P2/P1'],"T2/T1":w['T2/T1'],"rho2/rho1":w['rho2/rho1'],'M1n':Obl.Mn1w,'M2n':w['M2'],'M2':M2w}
        Data_s={"M1":Obl.M1,"Strong Shock Angle":np.rad2deg(Obl.theta_s),"P2/P1":s['P2/P1'],"T2/T1":s['T2/T1'],"rho2/rho1":s['rho2/rho1'],'M1n':Obl.Mn1s,'M2n':s['M2'],'M2':M2s}
        Data_aboutM1={ "Wave_angle_at_max_deflection":np.rad2deg(Obl.wave_angle_at_max_deflection),"Max_delfection_angle":np.rad2deg(Obl.max_delfection_angle)}
        return Data_w,Data_s,Data_aboutM1


        
