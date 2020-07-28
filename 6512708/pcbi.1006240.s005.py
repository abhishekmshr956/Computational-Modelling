# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:12:29 2017

@author: aborst
"""

import numpy as np
import Borst2018Library as bs
import matplotlib.pyplot as plt

init=1
saveswitch=0

# plot params

plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

mylabelsize=8
mytitlesize=10
mylegendsize=8
mylw=2

mycolor=np.array(['green','orange','brown'])
mylabel=np.array(['full model','NDS only','PDE only'])

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
    contrast_dep=np.zeros((3,6))

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
    
    return PDNDzscore, PDvar, NDvar

def plotresponse(EMDswitch,column,resp):   
    
    if EMDswitch==0:
        lower=0
        upper=5
        row=0
    if EMDswitch==1:
        lower=2
        upper=7
        row=2
    if EMDswitch==2:
        lower=1
        upper=6
        row=1
        
    nofbins=50
        
    resp[0]=resp[1]
    
    PDresp=resp[070:450]
    NDresp=resp[570:950]
    
    PDhist,edges=np.histogram(PDresp,bins=nofbins,range=(lower,upper), density=True)
    NDhist,edges=np.histogram(NDresp,bins=nofbins,range=(lower,upper), density=True)
    
    #---- timecourse ----    
    
    xpos=0.35+0.35*column
    ypos=0.65-0.18*row
    xsize=0.195
    ysize=0.13
    
    setmyaxes(xpos,ypos,xsize,ysize,0)
    
    plt.plot(myt,resp, color='black', label=mylabel[row])
    plt.xlabel('time [s]', fontsize=mylabelsize)
    plt.ylabel('response [mV]', fontsize=mylabelsize)
    plt.ylim(lower,upper)
    plt.legend(loc=1,fontsize=mylegendsize, frameon=False)
    plt.yticks(np.arange(1+upper-lower)+lower)
    
    #---- histo---
    
    xpos=0.55+0.35*column
    xsize=0.08
    
    setmyaxes(xpos,ypos,xsize,ysize,0)
    
    myy=np.linspace(lower,upper,nofbins)
    plt.fill_betweenx(myy,PDhist,color='blue')
    plt.fill_betweenx(myy,NDhist,color='red')

    if row==0:
        plt.legend(loc=1,fontsize=mylegendsize, frameon=False)
    else:
        plt.legend(loc=1,fontsize=mylegendsize, frameon=False)

    plt.ylim(lower,upper)
    plt.yticks(np.arange(1+upper-lower)+lower,'')
    plt.xticks(np.arange(5),'')
    
def motionnoise_detloop(cohindex):
    
    mystim=motionnoisestim[cohindex]

    returnswitch=0
    
    for i in range(3):
        EMDswitch=i      
        resp=bs.NewEMD(mystim,deltat,EMDswitch+3,returnswitch)
        PDNDzscore, PDvar, NDvar = detzscore(resp)
        print PDNDzscore, PDvar, NDvar
        plotresponse(EMDswitch,1,resp)
        
def photonnoise_detloop(nfac):
    
    returnswitch=0
    
    tf=np.zeros(1000)
    tf[050:450]=+1.0
    tf[550:950]=-1.0
    
    mystim=bs.calc_sinegrating(tf, img_rot=0, spat_freq=5, contrast=1.0)
    
    for i in range(3):
        EMDswitch=i
        resp=bs.NewEMD(mystim,deltat,EMDswitch+3,returnswitch,noisefac=nfac)
        PDNDzscore, PDvar, NDvar = detzscore(resp)
        print PDNDzscore, PDvar, NDvar
        plotresponse(EMDswitch,0,resp)
        
def go():
    if init==1:
        load_motionnoise_stimuli(2)
    plt.figure(figsize=(7.5,10.))
    photonnoise_detloop(4)
    motionnoise_detloop(4)
    if saveswitch==1:
        plt.savefig('Borst2018Fig5.tiff', dpi=600)
    
go()
    
        
        



        
    
