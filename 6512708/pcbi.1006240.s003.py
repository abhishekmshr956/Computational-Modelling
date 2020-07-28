# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 16:12:29 2017

@author: aborst
"""

import numpy as np
import Borst2018Library as bs
import matplotlib.pyplot as plt

maxtime=1000
deltat=10
myt=np.linspace(0,maxtime*deltat*0.001,maxtime)
noff=40
rect_level=0
mytf=1.0

saveswitch=0

# plot params

plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)

mylabelsize=8
mytitlesize=10
mylegendsize=6
mylw=1

xpos=0.40
ypos=0.80
xsize=0.50
ysize=0.10
ystep=0.12

def setmyaxes(myxpos,myypos,myxsize,myysize):

    ax=plt.axes([myxpos,myypos,myxsize,myysize])
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    

def plotresponse(myrow,resp):   
     
    xpos=0.3
    ypos=0.82-myrow*0.20
    xsize=0.6
    ysize=0.15
    
    setmyaxes(xpos,ypos,xsize,ysize,0)
    
    plt.plot(myt,resp, color='black')
    plt.xlabel('time [s]')
    plt.ylabel('response [mV]')
    
    
plt.figure(figsize=(7.5,10.))

tf=np.zeros(1000)
tf[100:400] = mytf
tf[500:800] = -mytf

stimulus=bs.calc_sinegrating(tf, img_rot=0, spat_freq=5, contrast=1.0)

lptau=50.0/deltat

R16=bs.rebin(stimulus,noff,noff,maxtime)
    
interim=bs.highpass(R16,250/deltat)
interim=interim+0.1*R16
L1=bs.rect(interim,0)   
lp=bs.lowpass(R16,lptau)
hp=L1

Eexc=+50.0
Einh=-20.0
gleak=1.0

Mi9=bs.rect(1.0-lp[:,0:noff-2,:],0)
Mi1=bs.rect(hp[:,1:noff-1,:],0)
Mi4=bs.rect(lp[:,2:noff-0,:],0)

gexc=bs.rect(Mi1,0)
ginh=bs.rect(Mi9+Mi4,0)

T4a=(Eexc*gexc+Einh*ginh)/(gexc+ginh+gleak)

T4a_rect=bs.rect(T4a,rect_level)
    
T4a_mean=np.mean(np.mean(T4a_rect,axis=0),axis=0)  

# first row

setmyaxes(xpos,ypos,xsize,ysize)

plt.plot(myt,R16[10,10,:],color='red', label='left',linewidth=mylw)
plt.plot(myt,R16[10,11,:],color='blue', label='central',linewidth=mylw)
plt.plot(myt,R16[10,12,:],color='green',label='right',linewidth=mylw)

plt.ylim(-0.5,1.5)
plt.text(-1.2,0.5,'input',fontsize=mylabelsize,rotation=90,verticalalignment='center')
plt.legend(loc=4, frameon=False, fontsize=mylegendsize)

# second row

setmyaxes(xpos,ypos-1*ystep,xsize,ysize)

plt.plot(myt,Mi9[10,10,:],color='red', label='left',linewidth=mylw)
plt.plot(myt,Mi1[10,10,:],color='blue', label='central',linewidth=mylw)
plt.plot(myt,Mi4[10,10,:],color='green', label='right',linewidth=mylw)

plt.ylim(-0.5,1.5)
plt.text(-1.2,0.5,'filtered input',fontsize=mylabelsize,rotation=90,verticalalignment='center')
plt.legend(loc=1, frameon=False, fontsize=mylegendsize)

# third row

setmyaxes(xpos,ypos-2*ystep,xsize,ysize)

plt.plot(myt,gexc[10,10,:],color='blue', label='gexc',linewidth=mylw)
plt.plot(myt,ginh[10,10,:],color='red', label='ginh',linewidth=mylw)

plt.ylim(-0.5,2)
plt.text(-1.2,0.75,'conductance [/gleak]',fontsize=mylabelsize,rotation=90,verticalalignment='center')
plt.legend(loc=1, frameon=False, fontsize=mylegendsize)

# forth row

setmyaxes(xpos,ypos-3*ystep,xsize,ysize)

plt.plot(myt,T4a[10,10,:],color='black', label='local',linewidth=mylw)

plt.ylim(-15,15)
plt.text(-1.2,0,'response [mV]',fontsize=mylabelsize,rotation=90,verticalalignment='center')
plt.legend(loc=1, frameon=False, fontsize=mylegendsize)

# fifth row

setmyaxes(xpos,ypos-4*ystep,xsize,ysize)

plt.plot(myt,T4a_mean,color='black', label='global',linewidth=mylw)

plt.ylim(-1,4)
plt.xlabel('time [s]',fontsize=mylabelsize)
plt.text(-1.2,1.5,'response [mV]',fontsize=mylabelsize,rotation=90,verticalalignment='center')
plt.legend(loc=1, frameon=False, fontsize=mylegendsize)

if saveswitch==1:
    plt.savefig('Borst2018Fig3.tiff', dpi=600)
        

    
    


        



        
    
