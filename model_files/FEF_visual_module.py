#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 08:58:41 2020

@author: amelie
"""
#This is used for simulations that change the time constants of all interneurons of one type at the same time.

from brian2 import *

from scipy import signal
from model_files.cells.RS_FEF import *
from model_files.cells.FS_FEF import *
from model_files.cells.SI_FEF import *
from model_files.cells.VIP_FEF_vis import *

def save_raster(name,raster_i,raster_t,path):
    raster_file=open(path+'/raster_'+name+'_i.txt','w')
    for elem in raster_i:
        raster_file.write(str(elem)+',')
    raster_file.close()
    raster_file=open(path+'/raster_'+name+'_t.txt','w')
    for elem in raster_t:
        raster_file.write(str(elem)+',')
    raster_file.close()
    return

def generate_visual_neurons(t_SI,t_FS,theta_phase,N_FS,N_RS,runtime,target_on,target_time):
    #print(target_on, theta_phase)
        
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms
    
    
    ##Define neuron groups
    RS=NeuronGroup(N_RS,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    RS.V = '-70*mvolt+10*rand()*mvolt'
    RS.h = '0+0.05*rand()'
    RS.m = '0+0.05*rand()'
    RS.mAR = '0.035+0.025*rand()'
    RS.J='55 * uA * cmeter ** -2'  

    FS=NeuronGroup(N_FS,eq_FS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS.V = '-70*mvolt+10*rand()*mvolt'
    FS.h = '0+0.05*rand()'
    FS.m = '0+0.05*rand()'
    FS.J='10 * uA * cmeter ** -2' #article=code=35
    
    SI=NeuronGroup(N_FS,eq_SI_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    SI.V = '-70*mvolt+10*rand()*mvolt'
    SI.h = '0+0.05*rand()'
    SI.m = '0+0.05*rand()' 
    SI.J='50 * uA * cmeter ** -2'  
    
    VIP=NeuronGroup(N_FS,eq_VIP_vis,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    VIP.V = '-90*mvolt+10*rand()*mvolt'
    VIP.Iapp='5 * uA * cmeter ** -2'    
    

    ##Synapses
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''
    
    def generate_syn(source,target,syntype,connection_pattern,g_i,taur_i,taud_i,V_i):
        S=Synapses(source,target,model=syntype+eq_syn,method='exact')
        if connection_pattern=='':
            S.connect()
        else :
            S.connect(condition=connection_pattern, skip_if_invalid=True)
        S.g_i=g_i
        S.taur_i=taur_i
        S.taud_i=taud_i
        S.V_i=V_i  
        return S
    
    
    #From E (granular layer) cells
    S_RSFS=generate_syn(RS,FS,'IsynRS_FEF_VM','i//10==j//10',0.2* msiemens * cm **-2,0.125*ms,1*ms,0*mV) 
    S_RSRS=generate_syn(RS,RS,'IsynRS_FEF_VM','i//10==j//10',0.2* msiemens * cm **-2,0.125*ms,1*ms,0*mV) 
    S_FSRS=generate_syn(FS,RS,'IsynSI_FEF_VM','i//10==j//10',0.2* msiemens * cm **-2,0.25*ms,t_FS,-80*mV) 
    S_FSFS=generate_syn(FS,FS,'IsynSI_FEF_VM','i//10==j//10',0.2* msiemens * cm **-2,0.25*ms,t_FS,-80*mV) 

    S_VIPSI=generate_syn(VIP,SI,'IsynSI_FEF_VM','i//10==j//10',0.7* msiemens * cm **-2,0.25*ms,20*ms,-80*mV) 

    S_SIVIP=generate_syn(SI,VIP,'IsynSI_FEF_VM','',0.01* msiemens * cm **-2,0.25*ms,t_SI,-80*mV) 

    S_SIRS=generate_syn(SI,RS,'IsynSI2_FEF_VM','i//10==j//10',1.5* msiemens * cm **-2,0.25*ms,t_SI,-80*mV) 
    
    
    def generate_spike_timing(N,f,start_time,end_time=runtime):
        list_time_and_i=[]
        for i in range(N):
            list_time=[(start_time,i)]
            next_spike=list_time[-1][0]+(1+0*rand())/f  #0.01
            while next_spike<end_time:
                list_time.append((next_spike,i))
                next_spike=list_time[-1][0]+(1+0*rand())/f
            list_time_and_i+=list_time
        return array(list_time_and_i)
    
    #Adding inputs
    #Background input     
    RS.ginp_RS=7.5* msiemens * cm **-2
    SI.ginp_SI=7.5* msiemens * cm **-2
    VIP.ginp_VIP_good=2.5* msiemens * cm **-2
    VIP.ginp_VIP_bad=2.5* msiemens * cm **-2
    if theta_phase=='good':
        fLIP=50*Hz
    else :
        fLIP=13*Hz
        
    gamma_background=generate_spike_timing(N_FS,fLIP,0*ms,end_time=3000*ms)
    
    if theta_phase=='mixed':
        t0=0*ms
        t1=125*ms
        gamma_background=generate_spike_timing(N_FS,fLIP,t0,end_time=t1)
        while t0+125*ms<runtime:
            fLIP=50*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==50*Hz)
            t0,t1=t0+125*ms,t1+125*ms
            gamma_background=vstack((gamma_background,generate_spike_timing(N_FS,fLIP,t0,end_time=t1)))

    Poisson_background = SpikeGeneratorGroup(N_FS, gamma_background[:,1], gamma_background[:,0]*second)
    S_in_bg_RS=Synapses(Poisson_background,RS,on_pre='Vinp=Vhigh')
    S_in_bg_RS.connect(j='i')
    S_in_bg_SI=Synapses(Poisson_background,SI,on_pre='Vinp=Vhigh')
    S_in_bg_SI.connect(j='i')
    S_in_bg_VIP=Synapses(Poisson_background,VIP,on_pre='Vinp=Vhigh')
    S_in_bg_VIP.connect(j='i')
        

    #Target input  
    gamma_target=generate_spike_timing(10,50*Hz,target_time,end_time=target_time+100*ms)
    if target_on :
        Poisson_target = SpikeGeneratorGroup(10, gamma_target[:,1], gamma_target[:,0]*second)
    else :
        Poisson_target = SpikeGeneratorGroup(10, [], []*second)
    S_in_target_VIP=Synapses(Poisson_target,VIP,on_pre='Vinp2=Vhigh')
    S_in_target_VIP.connect(j='i')
    S_in_target_SI=Synapses(Poisson_target,SI,on_pre='Vinp2=Vhigh')
    S_in_target_SI.connect(j='i')
    SI.ginp_SI2=2.5* msiemens * cm **-2
    RS.ginp_RS2=2.5* msiemens * cm **-2
    VIP.ginp_VIP2=3* msiemens * cm **-2

    
    #Define monitors and run network :
    R5=SpikeMonitor(RS,record=True)
    R6=SpikeMonitor(FS,record=True)
    R7=SpikeMonitor(VIP,record=True)
    R8=SpikeMonitor(SI,record=True)
    
    all_neurons=RS,FS,VIP,SI,Poisson_background,Poisson_target
    all_synapses=S_FSRS,S_RSFS,S_RSRS,S_FSFS,S_VIPSI,S_SIVIP,S_SIRS,S_in_bg_RS,S_in_bg_SI,S_in_bg_VIP,S_in_target_VIP,S_in_target_SI
    all_monitors=R5,R6,R7,R8
    
    return all_neurons,all_synapses,all_monitors

def sim_FEF_v_alone(simu,path,plot_raster=False):
    start_scope()   
    
    target_time,N_simu,t_SI,t_FS,theta_phase,g_LIP_FEF_v,target_on,runtime=simu[0],simu[1],simu[2],simu[3],simu[4],simu[5],simu[6],simu[7]
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    if not plot_raster :
        new_path=path+"/results_"+str(N_simu)
        os.mkdir(new_path)
    else :
        new_path=path
    
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
    
    N_FS,N_RS=20,20
    
    all_neurons,all_synapses,all_monitors=generate_visual_neurons(t_SI,t_FS,theta_phase,N_FS,N_RS,runtime,target_on,target_time)    
    
    net=Network()
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    
    
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp       
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2  
    
    taurinp3=2*ms
    taudinp3=40*ms
    tauinp3=taudinp3  
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)

    R5,R6,R7,R8=all_monitors

    save_raster('FEF RS v',R5.i,R5.t,new_path)
    save_raster('FEF FS v',R6.i,R6.t,new_path)
    save_raster('FEF VIP v',R7.i,R7.t,new_path)
    save_raster('FEF SOM v',R8.i,R8.t,new_path)

    
    if plot_raster :    
        figure(figsize=(4,5))
        plot(R5.t,R5.i+0,'r.',label='RS')
        plot(R6.t,R6.i+20,'b.',label='FS')
        plot(R7.t,R7.i+40,'k.',label='VIP')
        plot(R8.t,R8.i+60,'g.',label='SOM')
        xlim(0.2,runtime/second)
        xlabel('Time (s)')
        ylabel('Neuron index')

    return





if __name__=='__main__':
    close('all')
    start_scope()    
    
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    FLee=(0.05*mS/cm**2)/(0.4*uS/cm**2)*0.5
    theta_phase='mixed' #'good' or 'bad' or 'mixed'
#    theta_phase='good' #'good' or 'bad' or 'mixed'
#    theta_phase='bad' #'good' or 'bad' or 'mixed'
    runtime=1*second
    target_on=True
    target_time=650*msecond
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
    
    N_FS,N_RS=20,20
    all_neurons,all_synapses,all_monitors=generate_visual_neurons(theta_phase,N_FS,N_RS,runtime,target_on,target_time)    
    
    net=Network()
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp       
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2   
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)

    R5,R6,R7,R8=all_monitors
    
    figure(figsize=(4,5))
    plot(R5.t,R5.i+0,'r.',label='RS')
    plot(R6.t,R6.i+20,'b.',label='FS')
    plot(R7.t,R7.i+40,'k.',label='VIP')
    plot(R8.t,R8.i+60,'g.',label='SOM')
    xlim(0.2,runtime/second)
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    clear_cache('cython')