# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_FS_LIP='''
dV/dt=1/C_FS*(-J-Isyn-Igap-Iran-Iapp-IL-INa-IK) : volt
J : amp * meter ** -2
Isyn=IsynRS_LIP_sup+IsynFS_LIP_sup+IsynSI_LIP_sup+IsynRS_LIP_gran+IsynFS_LIP_gran+IsynIB_LIP+IsynSI_LIP_deep+Isyn_FEF+Isyn_mdPul : amp * meter ** -2
IsynRS_LIP_sup : amp * meter ** -2
IsynFS_LIP_sup : amp * meter ** -2
IsynSI_LIP_sup : amp * meter ** -2
IsynRS_LIP_gran : amp * meter ** -2
IsynFS_LIP_gran : amp * meter ** -2
IsynIB_LIP : amp * meter ** -2
IsynSI_LIP_deep : amp * meter ** -2
Isyn_FEF : amp * meter ** -2
Isyn_mdPul : amp * meter ** -2
Igap : amp * meter ** -2
IL=gL_FS*(V-VL_FS) : amp * meter ** -2
INa=gNa_FS*m0**3*h*(V-VNa_FS) : amp * meter ** -2
    m0=1/(1+exp((-V-38*mV)/10/mV)) : 1
    dh/dt=1/tauh*(hinf-h) : 1
    hinf=1/(1+exp((V+58.3*mV)/6.7/mV)) : 1
    tauh=0.225*ms+1.125*ms/(1+exp((V+37*mV)/15/mV)) : second
IK=gK_FS*m**4*(V-VK_FS) : amp * meter ** -2
    dm/dt=1/taum*(minf-m) : 1
    minf=1/(1+exp((-V-27*mV)/11.5/mV)) : 1
    taum=0.25*ms+4.35*ms*exp(-abs(V+10*mV)/10/mV) : second
    
Iran=sig_ranFS*randn(): amp * meter ** -2 (constant over dt)

Iapp=sinp*ginp_FS*(V-Vrev_inp) : amp * meter ** -2
    dsinp/dt=-sinp/taudinp + (1-sinp)/taurinp*0.5*(1+tanh(Vinp/10/mV)) : 1
    dVinp/dt=1/tauinp*(Vlow-Vinp) : volt
    ginp_FS = ginp_FS_good* int(sin(2*pi*t*8*Hz)>=0) + ginp_FS_bad* int(sin(2*pi*t*8*Hz)<0) : siemens * meter **-2 
    ginp_FS_good : siemens * meter **-2
    ginp_FS_bad : siemens * meter **-2  
'''


##Constants :
C_FS = 0.9* ufarad * cm ** -2

gL_FS=1 * msiemens * cm **-2
VL_FS=-65*mV

gNa_FS=200 * msiemens * cm **-2
VNa_FS=50*mV

gK_FS=20 * msiemens * cm **-2
VK_FS=-100*mV

sig_ranFS=0.05* mamp * cm **-2
sig_ranFS=0.05* mamp * cm **-2*0.5

if __name__=='__main__' :
    start_scope()
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp=0* msiemens * cm **-2
    
    FS=NeuronGroup(1,eq_FS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS.V = '-110*mvolt+10*rand()*mvolt'
    FS.h = '0+0.05*rand()'
    FS.m = '0+0.05*rand()'
    FS.J='5 * uA * cmeter ** -2'
    
    V1=StateMonitor(FS,'V',record=[0])
    
#    I1=StateMonitor(FS,'IL',record=[0])
#    I2=StateMonitor(FS,'INa',record=[0])
#    I3=StateMonitor(FS,'IK',record=[0])
    
    run(1*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('FS cell')
    
#    figure()
#    plot(I1.t/second,I1.IL[0],label='L')
#    plot(I1.t/second,I2.INa[0],label='Na')
#    plot(I1.t/second,I3.IK[0],label='K')
#    plot(I1.t/second,I4.IAR[0],label='AR')
#    title('Synaptic currents')
#    legend()