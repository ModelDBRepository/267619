#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 13:56:08 2020

@author: amelie
"""
#This is used for simulations that change the time constants of all interneurons of one type at the same time.

from brian2 import *

from scipy import signal

from model_files.LIP_full import *
from model_files.FEF_full import *

from itertools import *
from joblib import Parallel, delayed
import multiprocessing


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

def generate_syn(source,target,syntype,connection_pattern,g_i,taur_i,taud_i,V_i):
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''

    S=Synapses(source,target,model=syntype+eq_syn,method='exact')
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(j=connection_pattern, skip_if_invalid=True)
    S.g_i=g_i
    S.taur_i=taur_i
    S.taud_i=taud_i
    S.V_i=V_i  
    return S

def FEF_and_LIP(simu,path,plot_raster=False):
    prefs.codegen.target = 'numpy' 
    target_time,N_simu,t_SI,t_FS,theta_phase,g_LIP_FEF_v,target_on,runtime=simu[0],simu[1],simu[2],simu[3],simu[4],simu[5],simu[6],simu[7]
    
    if not plot_raster :
        new_path=path+"/results_"+str(N_simu)
        os.mkdir(new_path)
    else :
        new_path=path
    
    start_scope()
    close('all')
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2
    
    taurinp3=2*ms
    taudinp3=40*ms
    tauinp3=taudinp3

    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    N_SI,N_RS_gran,N_SI_gran=20,20,20
    N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_RS_vm,N_SI_vm=[20]*6

    
    all_SIdFSg=[2*msiemens * cm **-2] #1
    all_FSgRSg=[1* msiemens * cm **-2]
    all_RSgFSg=[1*msiemens * cm **-2]
    all_RSgRSg=[0.5*msiemens * cm **-2]
    all_FSgFSg=[0.3* msiemens * cm **-2]
    all_RSgRSs=[2*msiemens * cm **-2]
    all_RSgFSs=[0.1*msiemens * cm **-2]
    all_FSgRSs=[0.1* msiemens * cm **-2]
    all_J_RSg=['15 * uA * cmeter ** -2']
    all_J_FSg=['-5 * uA * cmeter ** -2']
    all_thal=[10* msiemens * cm **-2]
    thal=all_thal[0]
    
    all_syn_cond=list(product(all_SIdFSg,all_FSgRSg,all_RSgFSg,all_RSgRSg,all_FSgFSg,all_RSgRSs,all_RSgFSs,all_FSgRSs))
    all_J=list(product(all_J_RSg,all_J_FSg))
    syn_cond=all_syn_cond[0]
    J=all_J[0]
    
    print('Network setup')
    
    net=Network()
    
    all_neurons_LIP, all_synapses_LIP, all_gap_junctions_LIP, all_monitors_LIP=make_full_network(syn_cond,J,thal,t_SI,t_FS,theta_phase)
    V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd,R5,R6,R7,V5,V6,V7=all_monitors_LIP
    RS_sup_LIP,IB_LIP,SI_deep_LIP=all_neurons_LIP[0],all_neurons_LIP[5],all_neurons_LIP[9]
    RS_gran_LIP,FS_gran_LIP=all_neurons_LIP[7],all_neurons_LIP[8]
    
    all_neurons_FEF,all_synapses_FEF,all_monitors_FEF=create_FEF_full2(N_RS_vis,N_FS_vis,N_RS_mot,N_RS_vm,N_SI_vm,t_SI,t_FS,theta_phase,target_on,runtime,target_time)
    R8,R9,V_RS,V_SOM,inp_mon_FEF,R11,R12,R13,R14,mon_RS=all_monitors_FEF
    RSvm_FEF,SIvm_FEF,RSv_FEF,SIv_FEF,VIPv_FEF=all_neurons_FEF[0],all_neurons_FEF[1],all_neurons_FEF[4],all_neurons_FEF[7],all_neurons_FEF[6]
    
    IB_LIP.ginp_IB=0* msiemens * cm **-2 #the input to RS_sup_LIP is provided with synapses from FEF 
    SI_deep_LIP.ginp_SI=0* msiemens * cm **-2
#    RSvm_FEF.ginp_RS=0* msiemens * cm **-2
    SIvm_FEF.ginp_SI=0* msiemens * cm **-2 ####
    RSv_FEF.ginp_RS=0* msiemens * cm **-2
    SIv_FEF.ginp_SI=0* msiemens * cm **-2
    VIPv_FEF.ginp_VIP_good=0* msiemens * cm **-2
    VIPv_FEF.ginp_VIP_bad=0* msiemens * cm **-2
    
    if theta_phase=='good':
        RS_gran_LIP.ginp_RS_good=15* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_good=15* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_bad=15* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_bad=15* msiemens * cm **-2
    if theta_phase=='mixed':
        RS_gran_LIP.ginp_RS_good=2.5* msiemens * cm **-2
        RSvm_FEF.ginp_RS2_good=2.5* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_good=2.5* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_bad=5* msiemens * cm **-2
        RSvm_FEF.ginp_RS2_bad=5* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_bad=5* msiemens * cm **-2

    
    net.add(all_neurons_FEF)
    net.add(all_synapses_FEF)
    net.add(all_monitors_FEF)    
    
    net.add(all_neurons_LIP)
    net.add(all_synapses_LIP)
    net.add(all_gap_junctions_LIP)
    net.add(all_monitors_LIP)
    
    S_FEF_IB_LIP=generate_syn(RSvm_FEF,IB_LIP,'Isyn_FEF','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_FEF_SIdeep_LIP=generate_syn(RSvm_FEF,SI_deep_LIP,'Isyn_FEF','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_LIP_RS_FEF=generate_syn(RS_sup_LIP,RSvm_FEF,'Isyn_LIP','',0.009*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.009*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   

    # S_LIP_RS_FEF=generate_syn(RS_sup_LIP,RSvm_FEF,'Isyn_LIP','',0.01*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    # S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.01*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   

    S_LIP_RSv_FEF=generate_syn(RS_sup_LIP,RSv_FEF,'Isyn_LIP','',g_LIP_FEF_v,0.125*ms,1*ms,0*mV)   
    S_LIP_SIv_FEF=generate_syn(RS_sup_LIP,SIv_FEF,'Isyn_LIP','',0.025*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_VIPv_FEF=generate_syn(RS_sup_LIP,VIPv_FEF,'Isyn_LIP','',0.005*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
   
    RSv_FEF.ginp_RS2=2.5* msiemens * cm **-2
    SIv_FEF.ginp_SI2=2.5* msiemens * cm **-2
    VIPv_FEF.ginp_VIP2=2.5* msiemens * cm **-2
        
    net.add(S_FEF_IB_LIP)
    net.add(S_FEF_SIdeep_LIP)
    net.add(S_LIP_RS_FEF)
    net.add(S_LIP_FS_FEF)
    net.add(S_LIP_RSv_FEF)
    net.add(S_LIP_SIv_FEF)
    net.add(S_LIP_VIPv_FEF)
    
    print('Compiling with cython')
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2
    
    taurinp3=2*ms
    taudinp3=40*ms
    tauinp3=taudinp3
    
    noise_good=0* uA * cmeter ** -2
    noise_level=-30* uA * cmeter ** -2
    # noise_level=-40* uA * cmeter ** -2
    if theta_phase=='mixed':
        t0,t1=125*ms,250*ms
        i0,i1=int(t0//defaultclock.dt)+1,int(t1//defaultclock.dt)+1
        noise_array=ones((200000,20))* noise_good
        noise_array[i0:i1,:]=noise_level* rand(12500,20)
        while t0+250*ms<runtime:
            t0,t1=t0+250*ms,t1+250*ms
            i0,i1=int(t0//defaultclock.dt)+1,int(t1//defaultclock.dt)+1
            noise_array[i0:i1,:]=noise_level* rand(12500,20)
            
    elif theta_phase=='good':
        noise_array=ones((200000,20))* noise_good
    elif theta_phase=='bad':
        noise_array=ones((200000,20))* noise_level
#    print(noise_array)
    noise=TimedArray(noise_array,dt=defaultclock.dt)
    
#    defaultclock.dt = 0.01*ms
    net.run(runtime,report='text',report_period=300*second)
    
    print(RSvm_FEF.ginp_RS)
    print(RSvm_FEF.ginp_RS2)
    print(RSvm_FEF.Iinp1[:])
    print(RSvm_FEF.Iinp2[:])
    print(RSvm_FEF.J[:])
    
    save_raster('LIP RS',R1.i,R1.t,new_path)
    save_raster('LIP FS',R2.i,R2.t,new_path)
    save_raster('LIP SI',R3.i,R3.t,new_path)
    save_raster('LIP IB',R4.i,R4.t,new_path)
    save_raster('LIP RS gran',R5.i,R5.t,new_path)
    save_raster('LIP FS gran',R6.i,R6.t,new_path)
    save_raster('LIP SI deep',R7.i,R7.t,new_path)
    save_raster('FEF RS vm',R8.i,R8.t,new_path)
    save_raster('FEF SOM vm',R9.i,R9.t,new_path)
    save_raster('FEF RS v',R11.i,R11.t,new_path)
    save_raster('FEF FS v',R12.i,R12.t,new_path)
    save_raster('FEF VIP v',R13.i,R13.t,new_path)
    save_raster('FEF SI v',R14.i,R14.t,new_path)
    save_raster('FEF RS m',mon_RS.i,mon_RS.t,new_path)

    if plot_raster:
        #LIP Plot
        figure()
        plot(R1.t,R1.i+140,'r.',label='RS cells')
        plot(R2.t,R2.i+120,'b.',label='FS cells')
        plot(R3.t,R3.i+100,'g.',label='SOM cells')
        plot(R5.t,R5.i+70,'r.')
        plot(R6.t,R6.i+50,'b.')
        plot(R4.t,R4.i+20,'m.',label='IB cells')
        plot(R7.t,R7.i,'g.')
        xlim(0.2,runtime/second)
        legend(loc='upper left')
        xlabel('Time (s)')
        ylabel('Neuron index')
        title('LIP raster')
        
        #FEF Plots    
        figure(figsize=(10,4))
        subplot(131)
        title('Visual Neurons')
        plot(R11.t,R11.i+0,'r.',label='RS cells')
        plot(R12.t,R12.i+20,'b.',label='FS cells')
        plot(R13.t,R13.i+40,'k.',label='VIP cells')
        plot(R14.t,R14.i+60,'g.',label='SOM cells')
        xlim(0.2,runtime/second)
        legend(loc='upper left')   
        xlabel('Time (s)')
        ylabel('Neuron index')
        
        subplot(132)
        title('Visual-Motor Neurons')
        plot(R8.t,R8.i+60,'r.',label='RS cells')
        plot(R9.t,R9.i+40,'g.',label='SOM cells')
        xlim(0.2,runtime/second)
        legend(loc='upper left') 
        xlabel('Time (s)')
        ylabel('Neuron index')
        
        subplot(133)
        title('Decision cells')
        plot(mon_RS.t,mon_RS.i+0,'r.',label='RS cells')
        xlim(0.2,runtime/second)
        legend(loc='upper left') 
        xlabel('Time (s)')
        ylabel('Neuron index')

        tight_layout()
        suptitle('FEF raster')
        
        # LIP and FEF Plots
        figure(figsize=(4,9))
    #    subplot(411)
        up=200
        plot(R1.t,R1.i+140+up,'r.',label='RS cells')
        plot(R2.t,R2.i+120+up,'b.',label='FS cells')
        plot(R3.t,R3.i+100+up,'g.',label='SOM cells')
        plot([0.2,runtime/second],[95+up,95+up],'k--')
        plot(R5.t,R5.i+70+up,'r.')
        plot(R6.t,R6.i+50+up,'b.')
        plot([0.2,runtime/second],[45+up,45+up],'k--')
        plot(R4.t,R4.i+20+up,'m.',label='IB cells')
        plot(R7.t,R7.i+up,'g.')
        xlim(0.2,runtime/second)
        plot([0.2,runtime/second],[up-10,up-10],'k')
        
        up=40
        plot(R11.t,R11.i+0+up,'r.')
        plot(R12.t,R12.i+20+up,'b.')
        plot(R13.t,R13.i+40+up,'k.',label='VIP cells')
        plot(R14.t,R14.i+60+up,'g.')
        xlim(0.2,runtime/second)
        plot([0.2,runtime/second],[up-10,up-10],'k')

        up=140
        plot(R8.t,R8.i+20+up,'r.')
        plot(R9.t,R9.i+0+up,'g.')
        xlim(0.2,runtime/second)
        plot([0.2,runtime/second],[up-10,up-10],'k')

        plot(mon_RS.t,mon_RS.i+0,'r.')
        xlim(0.2,runtime/second)

        xlabel('Time (s)',fontsize=12)
        ylabel('Neuron index',fontsize=12)
        xticks(fontsize=12)
        yticks(fontsize=12)
        
        tight_layout()    
        subplots_adjust(bottom=0.1)   
        title('Full network raster')
        
        # figure()
        # plot(inp_mon_FEF.t,inp_mon_FEF.Iinp2[0],'r.',label='RS cells')
        
        figure()
        plot(inp_mon_FEF.t,inp_mon_FEF.Isyn[0])
        
#        show()
    
 #   clear_cache('cython')


if __name__=='__main__':
    path=""
    if os.name == 'nt':
        path=os.path.join(ntpath.dirname(os.path.abspath(__file__)),"results_"+str(datetime.datetime.now()).replace(':','-'))
    else :
        path="./results_"+str(datetime.datetime.now())
    
    os.mkdir(path)
    
    N=50
    liste_target_time=[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]
 
    liste_simus=[]
    for t in liste_target_time:
        liste_simus+=[t]*N
    
    liste_simus=[[liste_simus[i],i+750] for i in range(len(liste_simus))]
    
    simus_pas_faites=list(range(700))
    
    order=[50*i for i in range(len(liste_target_time))]
    for ind in range(1,50):
        order+=[50*i+ind for i in range(len(liste_target_time))]
    
    liste_simus=[liste_simus[i] for i in order if i in simus_pas_faites]

    liste_simus.reverse()
    
    print(liste_simus)
    
    print('Number of simulations: '+str(len(liste_simus)))
    
    for simu in liste_simus:
        FEF_and_LIP(simu,path)
        #clear_cache('cython')
    
#    clear_cache('cython')