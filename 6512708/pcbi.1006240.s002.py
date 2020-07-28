# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:12:29 2017

@author: aborst
"""

import numpy as np
import Borst2018Library as bs
import matplotlib.pyplot as plt

saveswitch=0

# plot params

plt.figure(figsize=(7.5,10.))

plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

mylabelsize=8
mytitlesize=10
mylegendsize=6
mylw=1

def setmyaxes(myxpos,myypos,myxsize,myysize):
    ax=plt.axes([myxpos,myypos,myxsize,myysize])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

maxtime=1000

# T4 parameters
    
Eexc=+50.0
Einh=-10.0
gleak=1.0
deltat=1
memcap=10.0
myt=np.linspace(0,maxtime*0.001*deltat,maxtime)

gin=np.zeros(maxtime)
gIh=np.zeros(maxtime)  
gexc=np.zeros(maxtime)
ginh=np.zeros(maxtime)
Vm=np.zeros(maxtime)

def InhibExcitCell(inhinput,excinput):
    
    global Vm, gexc, ginh
    
    gexc=gexc*0.0
    ginh=ginh*0.0+1.0
    
    ginh[100:300]=ginh[100:300]-inhinput
    ginh=bs.rect(ginh,0)
    
    gexc[100:300]=excinput
    gexc=bs.rect(gexc,0)
    
    Vm=0.0*Vm
    
    for t in range(1,maxtime,1):
        
        Vm[t]=Eexc*gexc[t]+Einh*ginh[t]+Vm[t-1]*memcap/deltat
        gin[t]=gleak+gexc[t]+ginh[t]+memcap/deltat
        Vm[t]=Vm[t]/gin[t]
        Vm[t]=bs.rect(Vm[t],-100)
        Vm[t]=bs.ceil(Vm[t],+100)
        
    return Vm
    
def ExcitExcitCell(inp1,inp2):
    
    global Vm, gexc
    
    gexc=gexc*0.0
    
    gexc[100:300]=inp1+inp2
    gexc=bs.rect(gexc,0)
    
    Vm=0.0*Vm
    
    for t in range(1,maxtime,1):

        Vm[t]=Eexc*gexc[t]+Vm[t-1]*memcap/deltat
        gin[t]=gleak+gexc[t]+memcap/deltat
        Vm[t]=Vm[t]/gin[t]
        Vm[t]=bs.rect(Vm[t],-100)
        Vm[t]=bs.ceil(Vm[t],+100)
        
    return Vm
    
def plottc(x):
    
    plt.plot(myt,x[0],color='black', linewidth=mylw, label='Inp x')
    plt.plot(myt,x[1], color='green',linewidth=mylw, label='Inp y')
    plt.plot(myt,x[2], color='red',linewidth=mylw, label='Inp x & y')
    plt.plot(myt,x[3], color='blue',linewidth=mylw, label='Lin Exp')
    plt.legend(loc=1,frameon=False,fontsize=mylegendsize)
    
def plotNL(g,R1,R2,SeqR,LinR):
    
    plt.plot(g,R1, linewidth=mylw, color='black', label='Inp x')
    plt.plot(g,R2, linewidth=mylw, color='green', label='Inp y')
    plt.plot(g,SeqR, linewidth=mylw, color='red', label='Inp x & y')
    plt.plot(g,LinR, linewidth=mylw, color='blue', label='Lin Exp')
    plt.legend(loc=2,frameon=False,fontsize=mylegendsize)
    
def myfigure():
    
    x=np.zeros((2,4,1000))
    
    xsize=0.17
    ysize=0.12
    
    myx=np.array([0.30, 0.55, 0.80])
    myy=np.array([0.75, 0.55])
    
    # -------- left column ------------
    
    # 1st dim: Cell1, Cell2
    # 2nd dim: 1,2,12, linexp
    # 3rd dim: time
    
    x[0,0,:]=ExcitExcitCell(1,0)        
    x[0,1,:]=ExcitExcitCell(0,1)*1.05
    x[0,2,:]=ExcitExcitCell(1,1)
    
    x[1,0,:]=InhibExcitCell(1,0)        
    x[1,1,:]=InhibExcitCell(0,1)
    x[1,2,:]=InhibExcitCell(1,1)  
    
    for i in range(3):
        x[1,i,0:50]=x[1,i,50]
    
    for i in range(2):
        x[i,3,:]=x[i,0,:]+x[i,1,:]-x[i,0,90]
    
    myylims=np.array([-10,60,-10,30])
    
    for i in range(2):

        setmyaxes(myx[0],myy[i],xsize,ysize)
        plottc(x[i])
        plt.ylim(myylims[i*2],myylims[i*2+1])

        if i==0:
            plt.text(-0.3,25,'response [mV]',fontsize=mylabelsize,rotation=90,verticalalignment='center')
        if i==1:
            plt.text(-0.3,10,'response [mV]',fontsize=mylabelsize,rotation=90,verticalalignment='center')
            
        plt.xlabel('time [s]',fontsize=mylabelsize)
    
    # -------- upper row ------------
    
    g=np.linspace(0,1,11)
    R1=Eexc*g/(g+1)
    R2=Eexc*g/(g+1)+1
    
    LinR=Eexc*2.0*g/(g+1)
    SeqR=Eexc*2.0*g/(2*g+1)
    NonR=Eexc*(-2.0*g**2)/(2.0*g**2+3.0*g+1)
    
    setmyaxes(myx[1],myy[0],xsize,ysize)
    plotNL(g,R1,R2,SeqR,LinR)
    plt.text(-0.3,30,'response [mV]',fontsize=mylabelsize,rotation=90,verticalalignment='center')
    plt.xlabel('conductance [1/g leak]',fontsize=mylabelsize)
    plt.ylim(0,60)
    
    setmyaxes(myx[2],myy[0],xsize,ysize)
    plt.plot(g,NonR, color='green',linewidth=mylw)
    plt.plot(g,g*0, color='black',linewidth=mylw)
    plt.text(-0.3,-5,'nonlin response [mV]',fontsize=mylabelsize,rotation=90,verticalalignment='center')
    plt.xlabel('conductance [1/g leak]',fontsize=mylabelsize)
    plt.ylim(-20,10)
    
    #----------- lower row -------------
    
    R1=g*(-Einh)/(2.0*(2-g))
    R2=g*(2.0*Eexc-Einh)/(2.0*(2+g))

    LinR=(Eexc*(2-g)-2*Einh)*g/(4-g**2)
    SeqR=(Eexc-Einh)*g/2.0
    NonR=(Eexc*(2-g)+g*Einh)*g**2/(2*(4-g**2))
    
    setmyaxes(myx[1],myy[1],xsize,ysize)
    plotNL(g,R1,R2,SeqR,LinR)
    plt.text(-0.3,15,'response [mV]',fontsize=mylabelsize,rotation=90,verticalalignment='center')
    plt.xlabel('conductance [1/g leak]',fontsize=mylabelsize)
    plt.ylim(0,30)
    
    setmyaxes(myx[2],myy[1],xsize,ysize) 
    plt.plot(g,NonR, color='green',linewidth=mylw)
    plt.plot(g,g*0, color='black',linewidth=mylw)
    plt.xlabel('conductance [1/g leak]',fontsize=mylabelsize)
    plt.text(-0.3,4,'nonlin response [mV]',fontsize=mylabelsize,rotation=90,verticalalignment='center')
    plt.ylim(-2,10)

myfigure()

if saveswitch==1:
    plt.savefig('Borst2018Fig2.tiff', dpi=600)
        
    
