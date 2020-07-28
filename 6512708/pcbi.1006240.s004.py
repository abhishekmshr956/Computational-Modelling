# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:12:29 2017

@author: aborst
"""

import numpy as np
import Borst2018Library as bs
import matplotlib.pyplot as plt

init=1
nofexamples=5
saveswitch=0

# plot params

plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

mylabelsize=8
mytitlesize=10
mylegendsize=8
mylw=2

mycolor=np.array(['green','orange','brown'])
mylabel=np.array(['full model','PDE only','NDS only'])

def setmyaxes(myxpos,myypos,myxsize,myysize,polar):
    if polar==1:    
        ax=plt.axes([myxpos,myypos,myxsize,myysize],projection='polar')
    else:
        ax=plt.axes([myxpos,myypos,myxsize,myysize])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    
# other params

if init==1:

    maxtime=1000
    deltat=10
    myt=np.linspace(0,maxtime*deltat*0.001,maxtime)
    motionnoisestim=np.zeros((6,200,200,1000))
    myzscore=np.zeros((6,6))
    
    PD_tf_tuning=np.zeros((3,7))
    ND_tf_tuning=np.zeros((3,7))
    dir_tuning=np.zeros((3,13))
    photonnoise_dep=np.zeros((3,6))
    motionnoise_dep=np.zeros((3,6))

def loadstimulus(example,i):
    
    mydir='Coherent Noise Stimuli/example'+np.str(example)+'/'
    
    if i==0: fname='randnoise00coh.npy'
    if i==1: fname='randnoise20coh.npy'
    if i==2: fname='randnoise40coh.npy'
    if i==3: fname='randnoise60coh.npy'
    if i==4: fname='randnoise80coh.npy'
    if i==5: fname='randnoise100coh.npy'
    
    print mydir+fname
    
    mystim=np.load(mydir+fname)
    
    return mystim
    
def load_motionnoise_stimuli(example):
    for i in range(6):
        motionnoisestim[i]=loadstimulus(example,i)
        
def calc_motionnoise_stimuli():
    
    for i in range(6):
        perccoher=i*20.0
        blurrindex=1
        motionnoisestim[i]=bs.calc_PDND_motion_noise(maxtime,perccoher,blurrindex)

def save_motionnoise_stimuli(example):
    
    mydir='Coherent Noise Stimuli/example'+np.str(example)+'/'
    
    calc_motionnoise_stimuli()
    
    for i in range(6):
        if i==0: fname='randnoise00coh.npy'
        if i==1: fname='randnoise20coh.npy'
        if i==2: fname='randnoise40coh.npy'
        if i==3: fname='randnoise60coh.npy'
        if i==4: fname='randnoise80coh.npy'
        if i==5: fname='randnoise100coh.npy'
        
        fname=mydir+fname
        print fname
        np.save(fname,motionnoisestim[i])
        
def detzscore(resp):
    
    PDresp=resp[070:450]
    NDresp=resp[570:950]

    PDmean=np.mean(PDresp)
    PDvar=np.mean((PDresp-PDmean)**2)

    NDmean=np.mean(NDresp)
    NDvar=np.mean((NDresp-NDmean)**2)

    PDNDzscore=np.abs((PDmean-NDmean))/np.sqrt(PDvar+NDvar)

    PDNDzscore=np.int(100*PDNDzscore)/100.0
    
    return PDNDzscore
        
# -------- response calculations ------------
        
def calc_dir_tuning():
    
    myvelo=np.zeros(300)
    myvelo[50:300]=1.0
    
    for w in range(13):
        print
        img_rot=w*30
        print 'angle=', img_rot
        movie=bs.calc_sinegrating(myvelo,img_rot,5,1.0)
        
        for i in range(3):
        
            output=bs.NewEMD(stimulus=movie,deltat=10,det_switch=3+i,ret_switch=0)
            dir_tuning[i,w]=np.mean(output[150:300])-np.mean(output[0:50])
    
    for i in range(3):
        dir_tuning[i,:]=bs.normalize(dir_tuning[i,:])
    
def calc_tf_tuning():
    
    maxvelo=np.array([0.1,0.2,0.5,1.0,2.0,5.0,10.0])
    velofct=np.zeros(300)
    velofct[50:300]=1.0

    for i in range(7):
        movie=bs.calc_sinegrating(1*velofct*maxvelo[i],0,5,1)
        for j in range(3):
            output=bs.NewEMD(stimulus=movie,deltat=10,det_switch=3+j,ret_switch=0)
            PD_tf_tuning[j,i]=np.mean(output[150:300])-np.mean(output[0:50])
                
    for i in range(7):
        movie=bs.calc_sinegrating((-1)*velofct*maxvelo[i],0,5,1)
        for j in range(3):
            output=bs.NewEMD(stimulus=movie,deltat=10,det_switch=3+j,ret_switch=0)
            ND_tf_tuning[j,i]=np.mean(output[150:300])-np.mean(output[0:50])
        
    for i in range(3):
        ND_tf_tuning[i,:]=ND_tf_tuning[i,:]/np.max(PD_tf_tuning[i,:])
        PD_tf_tuning[i,:]=PD_tf_tuning[i,:]/np.max(PD_tf_tuning[i,:])
    
def calc_motionnoise_dep():  
    
    for k in range(nofexamples):
        
        example=k+1
        
        print
        print example
        
        if init==1:
            load_motionnoise_stimuli(example)
    
        for i in range(6):
            
            mystim=motionnoisestim[i]*1.0
            
            for j in range(3):
        
                resp=bs.NewEMD(mystim,20,j+3,0,0)
                resp[0]=resp[1]
                resp=bs.normalize(resp-np.min(resp))
            
                PDNDzscore=detzscore(resp)
        
                motionnoise_dep[j,i]+=PDNDzscore
                
    motionnoise_dep[:,:]=motionnoise_dep[:,:]/(1.0*nofexamples)
   
def calc_photonnoise_dep():  
    
    tf=np.zeros(1000)
    tf[050:450]=+1.0
    tf[550:950]=-1.0
    
    meanlum=np.array([1,2,4,8,16,32])
    
    mystim=bs.calc_sinegrating(tf, img_rot=0, spat_freq=5, contrast=1.0)
    
    for i in range(6):
        
        nfac=meanlum[i]
    
        for j in range(3):
            
            resp=bs.NewEMD(mystim,20,j+3,0,noisefac=nfac)
            resp[0]=resp[1]
            resp=bs.normalize(resp-np.min(resp))
        
            PDNDzscore=detzscore(resp)
        
            photonnoise_dep[j,i]=PDNDzscore
    
def calc_all():  

    calc_tf_tuning()
    calc_dir_tuning()
    calc_photonnoise_dep()
    calc_motionnoise_dep()
    
# ---- plotting -----------------
    
def plot_tf_tuning(PDND):

    maxvelo=np.array([0.1,0.2,0.5,1.0,2.0,5.0,10.0])
    
    for i in range(3):
        if PDND==0:
            plt.plot(maxvelo, PD_tf_tuning[i,:], linewidth=mylw, color=mycolor[i],label=mylabel[i])
            plt.scatter(maxvelo, PD_tf_tuning[i,:],color=mycolor[i])
        if PDND==1:
            plt.plot(maxvelo, ND_tf_tuning[i,:], linewidth=mylw, color=mycolor[i],label=mylabel[i])
            plt.scatter(maxvelo, ND_tf_tuning[i,:],color=mycolor[i])
    
    plt.plot(maxvelo,maxvelo*0, color='black')
    plt.xlabel('temp. frequency [Hz]',fontsize=8)
    if PDND==0: plt.ylabel('response',fontsize=8)
    plt.xscale('log')
    plt.xlim(0.1,10)
    plt.ylim(-0.1,1.1)
           
def plot_dir_tuning():
        
        myalpha=np.linspace(0,360,13)
        
        for i in range(3):
            plt.polar(myalpha/360.0*2*np.pi,bs.rect(dir_tuning[i,:],0),linewidth=mylw,color=mycolor[i])
        plt.ylim(0,1.2)
        plt.yticks(np.arange(3)*0.5)
        locs, labels=plt.yticks()
        plt.yticks(locs, ('','',''))
        
def plot_photonnoise_dep():
    
    myx=np.linspace(0,5,6)
    
    for i in range(3):
        
        plt.plot(myx,photonnoise_dep[i,:]**2,linewidth=mylw,color=mycolor[i],label=mylabel[i])
        plt.scatter(myx,photonnoise_dep[i,:]**2,color=mycolor[i])
    
    plt.xlabel('log mean luminance [arb.units]',fontsize=8)
    plt.xlim(-0.5,5.5)
    plt.ylabel('signal-to-noise ratio',fontsize=8)
    plt.yscale('log')
    plt.ylim(0.1,100)
               
def plot_motionnoise_dep():
    
    myx=np.linspace(0,100,6)
    
    for i in range(3):
        
        plt.plot(myx,motionnoise_dep[i,:]**2,linewidth=mylw,color=mycolor[i],label=mylabel[i])
        plt.scatter(myx,motionnoise_dep[i,:]**2,color=mycolor[i])
    
    plt.xlabel('coherence [%]',fontsize=8)
    plt.xlim(-10,110)
    plt.yscale('log')
    plt.ylim(0.1,1000)
         
def plot_all():
    
    plt.figure(figsize=(7.5,10.))
    
    xsize=0.18
    ysize=0.12
    xoffset=0.332
    yoffset=0.37
    
    setmyaxes(xoffset+0.02,yoffset+0.35,xsize,ysize,0)
    plot_tf_tuning(0)
    setmyaxes(xoffset+0.25,yoffset+0.35,xsize,ysize,0)
    plot_tf_tuning(1)
    setmyaxes(xoffset+0.02,yoffset+0.05,xsize,ysize,0)
    plot_photonnoise_dep()
    setmyaxes(xoffset+0.25,yoffset+0.05,xsize,ysize,0)
    plot_motionnoise_dep()
    setmyaxes(xoffset+0.47,yoffset+0.18,xsize,ysize,1)
    plot_dir_tuning()
    
# -----------------------------------------
    
def go():
    
    calc_all()
    plot_all()
        
    if saveswitch==1:
        plt.savefig('Borst2018Fig4.tiff', dpi=600)
    
go()

        



        
    
