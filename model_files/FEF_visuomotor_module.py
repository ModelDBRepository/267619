# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""

#This is used for simulations that change the time constants of all interneurons of one type at the same time.

from brian2 import *

from scipy import signal
from model_files.cells.RS_FEF_VM import *
from model_files.cells.FS_FEF import *
from model_files.cells.SI_FEF_more_h import *
from model_files.cells.VIP_FEF import *

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

def zeros_ones_monitor(spikemon,record_dt,runtime):
    L=int(runtime/record_dt)
    zeros_ones=[0]*L
    for time in spikemon.t:
        zeros_ones[int(time/record_dt)]+=1
    return zeros_ones

def generate_deepSI_and_gran_layers(t_SI,t_FS,theta_phase,N_RS,N_SOM,runtime):
    
    if theta_phase=='bad':
        LIP_input=3* msiemens * cm **-2
        
    if theta_phase=='good' or theta_phase=='mixed':
        LIP_input=5* msiemens * cm **-2
        # LIP_input=6* msiemens * cm **-2
        

    
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms
    
    #Single column network
    
    ##Define neuron groups
    RS=NeuronGroup(N_RS,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    RS.V = '-70*mvolt+10*rand()*mvolt'
    RS.h = '0+0.05*rand()'
    RS.m = '0+0.05*rand()'
    RS.mAR = '0.035+0.025*rand()'
    RS.J_fixed='45 * uA * cmeter ** -2'
    # RS.J_fixed='80 * uA * cmeter ** -2'
    # RS.J_fixed='75 * uA * cmeter ** -2'
    
    SOM=NeuronGroup(N_SOM,eq_SI_FEF_more_h,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    SOM.V = '-110*mvolt+10*rand()*mvolt'
    SOM.h = '0+0.05*rand()'
    SOM.m = '0+0.05*rand()'
    # SOM.J='40 * uA * cmeter ** -2'
    SOM.J='40 * uA * cmeter ** -2'
    
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
            S.connect(j=connection_pattern, skip_if_invalid=True)
        S.g_i=g_i
        S.taur_i=taur_i
        S.taud_i=taud_i
        S.V_i=V_i  
        return S
    
    
    #From RS cells
    S_RSRS=generate_syn(RS,RS,'IsynRS_FEF_VM','',0.6*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.4
    S_RSSOM=generate_syn(RS,SOM,'IsynRS_FEF_VM','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.6

    #From SOM cells,
    S_SOMRS=generate_syn(SOM,RS,'IsynSI_FEF_VM','',0.8*msiemens * cm **-2,0.25*ms,t_SI,-80*mV) #0.35
    
    def generate_spike_timing(N,f,start_time,end_time=runtime):
        list_time_and_i=[]
        for i in range(N):
            list_time=[(start_time,i)]
            next_spike=list_time[-1][0]+(1+0.01*rand())/f
            while next_spike<end_time:
                list_time.append((next_spike,i))
                next_spike=list_time[-1][0]+(1+0.01*rand())/f
            list_time_and_i+=list_time
        return array(list_time_and_i)

    #Defining inputs
    #mdPul input
    RS.ginp_RS2_good=2.5* msiemens * cm **-2
    RS.ginp_RS2_bad=5* msiemens * cm **-2
    # RS.ginp_RS2_good=5* msiemens * cm **-2
    # RS.ginp_RS2_bad=5* msiemens * cm **-2
    # RS.ginp_RS2_good=0* msiemens * cm **-2
    # RS.ginp_RS2_bad=0* msiemens * cm **-2
    fmdPul=13*Hz
    # fmdPul=30*Hz
    spikes_mdPul=generate_spike_timing(N_RS,fmdPul,0*ms,end_time=3000*ms)
    #Theta=4Hz
    theta_frequency=4*Hz
    if theta_phase=='mixed':
        t0=0*ms
        # t1=125*ms
        t1=0.5/theta_frequency
        spikes_mdPul=generate_spike_timing(N_SOM,fmdPul,t0,end_time=t1)
        while t0+1/theta_frequency<runtime:
            t0,t1=t0+1/theta_frequency,t1+1/theta_frequency
            spikes_mdPul=vstack((spikes_mdPul,generate_spike_timing(N_SOM,fmdPul,t0,end_time=t1)))
            
    # print(spikes_mdPul[:,0:10])
    
    G_in_mdPul = SpikeGeneratorGroup(N_RS, spikes_mdPul[:,1], spikes_mdPul[:,0]*second)
    
    S_in_mdPul=Synapses(G_in_mdPul,RS,on_pre='Vinp2=Vhigh')
    S_in_mdPul.connect(j='i')
         
    
    #LIP input
    SOM.ginp_SI=LIP_input
    RS.ginp_RS=LIP_input

    if theta_phase=='good' or theta_phase=='mixed':
        fLIP=50*Hz
        # fLIP=150*Hz
    else :
        fLIP=13*Hz

    LIPspikes=generate_spike_timing(N_SOM,fLIP,0*ms,end_time=2100*ms)

    if theta_phase=='mixed':
        t0=0*ms
        t1=0.5/theta_frequency
        fLIP=50*Hz
        # fLIP=150*Hz
        LIPspikes=generate_spike_timing(N_SOM,fLIP,t0,end_time=t1)
        while t0+1/theta_frequency<runtime:
            t0,t1=t0+0.5/theta_frequency,t1+0.5/theta_frequency
            fLIP=50*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==50*Hz)
            # fLIP=150*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==150*Hz)
            LIPspikes=vstack((LIPspikes,generate_spike_timing(N_SOM,fLIP,t0,end_time=t1)))
             
    G_in_LIP = SpikeGeneratorGroup(N_SOM, LIPspikes[:,1], LIPspikes[:,0]*second)
    S_LIP_in=Synapses(G_in_LIP,SOM,on_pre='Vinp=Vhigh')
    S_LIP_in.connect(j='i')
    
    S_LIP_in2=Synapses(G_in_LIP,RS,on_pre='Vinp=Vhigh')
    S_LIP_in2.connect(j='i')
    
    #Define monitors and run network :
    R5=SpikeMonitor(RS,record=True)
    R6=SpikeMonitor(SOM,record=True)
    
    inpmon=StateMonitor(RS,'ginp_RS2',record=True)
    # inpmon=StateMonitor(RS,'Isyn',record=True)
    
    V_RS=StateMonitor(RS,'V',record=True)
    V_SOM=StateMonitor(SOM,'V',record=True)
    
    all_neurons=RS,SOM,G_in_mdPul,G_in_LIP
    all_synapses=S_RSRS,S_RSSOM,S_SOMRS,S_in_mdPul,S_LIP_in,S_LIP_in2
    all_monitors=R5,R6,V_RS,V_SOM,inpmon
    
    return all_neurons,all_synapses,all_monitors

def sim_FEF_vm_alone(simu,path,plot_raster=False):
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
#    taurinp=2*ms
#    taudinp=10*ms
#    tauinp=taudinp   
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
    
    N_RS,N_SOM=20,20
    all_neurons,all_synapses,all_monitors=generate_deepSI_and_gran_layers(t_SI,t_FS,theta_phase,N_RS,N_SOM,runtime)    
    
    net=Network()
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    
#    taurinp=2*ms
#    taudinp=10*ms
#    tauinp=taudinp
    
#    taurinp2=0.1*ms
#    taudinp2=0.5*ms
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2   
    
    taurinp3=2*ms
    taudinp3=40*ms
    tauinp3=taudinp3   
    
    noise_good=0* uA * cmeter ** -2
    noise_level=-30* uA * cmeter ** -2
    # noise_level=0* uA * cmeter ** -2
    theta_frequency=4*Hz
    if theta_phase=='mixed':
        t0,t1=0.5/theta_frequency,1/theta_frequency
        i0,i1=int(t0//defaultclock.dt)+1,int(t1//defaultclock.dt)+1
        noise_array=ones((200000,20))* noise_good
        noise_array[i0:i1,:]=noise_level* rand(int(0.5/theta_frequency*100000),20)
        while t0+1/theta_frequency<runtime:
            t0,t1=t0+1/theta_frequency,t1+1/theta_frequency
            i0,i1=int(t0//defaultclock.dt)+1,int(t1//defaultclock.dt)+1
            noise_array[i0:i1,:]=noise_level* rand(int(0.5/theta_frequency*100000),20)
#    print(noise_array)
            
    elif theta_phase=='bad':
        noise_array=ones((200000,20))* noise_level
    elif theta_phase=='good':
        noise_array=ones((200000,20))* noise_good
    noise=TimedArray(noise_array,dt=defaultclock.dt)
    
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

    # taurinp3=2.5*ms
    # taudinp3=40*ms
    # tauinp3=taudinp3
    
    net.run(runtime,report='text',report_period=300*second)

    R5,R6,V_RS,V_SOM,inpmon=all_monitors
    
    RS=all_neurons[0]
    # print(RS.ginp_RS)
    # print(RS.ginp_RS2)
    # print(RS.Iinp1[:])
    # print(RS.Iinp2[:])
    # print(RS.J[:])

    save_raster('FEF RS vm',R5.i,R5.t,new_path)
    save_raster('FEF SOM vm',R6.i,R6.t,new_path)

    
    if plot_raster :
        figure()
        plot(R5.t,R5.i+0,'r.',label='RS cells')
        plot(R6.t,R6.i+20,'g.',label='SOM cells')
        xlim(0,runtime/second)
    #    legend(loc='upper left')
        xlabel('Time (s)')
        ylabel('Neuron index')
        ylim(-1,41)
        
        N_RS_spikes=zeros_ones_monitor(R5,defaultclock.dt,runtime)
        f, t, Sxx = signal.spectrogram(array(N_RS_spikes), 100000*Hz,nperseg=20000,noverlap=15000)
        figure()
        pcolormesh(t, f, Sxx)#, cmap=)
        colorbar(format='%.1e')
        ylabel('Frequency (Hz)',fontsize=12)
        xlabel('Time (s)',fontsize=12)
        xticks(fontsize=12)
        yticks(fontsize=12)
        ylim(0,45)
        title('Power ($V^2$)',fontsize=12)
        
        figure()
        plot(inpmon.t,inpmon.ginp_RS2[0])
        
        # figure()
        # plot(inpmon.t,inpmon.Iapp[0])
        mean_V_RS=1/20*sum(V_RS.V,axis=0)
        mean_V_file = open(new_path+"/mean_V_RS.txt", "w")
        savetxt(mean_V_file, mean_V_RS)
        mean_V_file.close()
        
    return

if __name__=='__main__':
    close('all')
    start_scope()    
    
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms

    theta_phase='mixed' #'good' or 'bad' or 'mixed'
    runtime=2*second
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp 
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
    
    N_RS,N_SOM=20,20
    all_neurons,all_synapses,all_monitors=generate_deepSI_and_gran_layers(t_SI,t_FS,theta_phase,N_RS_gran,N_SOM,runtime)    
    
    net=Network()
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2   
    
    taurinp3=2*ms
    taudinp3=40*ms
    tauinp3=taudinp3   
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)

    R5,R6,V_RS,V_SOM,inpmon=all_monitors
    
    figure()
    plot(inpmon.t,inpmon.Iinp2[0])
    xlabel('Time (s)')
    ylabel('Iinp2')
    tight_layout()
    
    figure()
    plot(R5.t,R5.i+20,'r.',label='RS')
    plot(R6.t,R6.i+40,'g.',label='SOM')
    xlim(0,runtime/second)
#    legend(loc='upper left')
    xlabel('Time (s)')
    ylabel('Neuron index')
    ylim(-1,61)
    
#    figure()
#    plot(V_RS.t,V_RS.V[0])
    
    min_t=int(50*ms*100000*Hz)
    LFP_V_RS=1/20*sum(V_RS.V,axis=0)[min_t:]
    LFP_V_SOM=1/20*sum(V_SOM.V,axis=0)[min_t:]
    
    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_SOM=signal.periodogram(LFP_V_SOM, 100000,'flattop', scaling='spectrum')
    
    figure()
    subplot(221)
    plot((V_RS.t/second)[min_t:],LFP_V_RS)
    ylabel('LFP')
    title('gran RS cell')
    subplot(223)
    plot((V_SOM.t/second)[min_t:],LFP_V_SOM)
    ylabel('LFP')
    title('gran FS cell')
    
    subplot(222)
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('gran RS cell')
    subplot(224)
    plot(f,Spectrum_LFP_V_SOM)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('gran FS cell')
    
    figure()
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum')
    xlabel('Frequency (Hz)')
    xlim(0,50)
    

    f, t, Sxx = signal.spectrogram(LFP_V_RS, 100000*Hz,nperseg=20000,noverlap=15000)
    figure()
    pcolormesh(t, f, Sxx)#, cmap=)
    colorbar(format='%.1e')
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    ylim(0,45)
    title('Power ($V^2$)')
    
    clear_cache('cython')