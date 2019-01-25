# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 09:53:55 2018

@author: chp16hcw
"""

from __future__ import division
import csv, numpy, glob, os
import matplotlib.pyplot as plt
from scipy import integrate

os.chdir("C:\\Users\\chp16hcw\\Google Drive\\PhD\\Results\\Ktest\\ATR\\Raw\\")

def findmax(InArr):
    Amax = max(InArr)
    iMax = [i for i,x in enumerate(InArr) if x == Amax]
    iMax = iMax[0]
    return Amax, iMax

def findslope(t,T):
    mc,res,_,_,_ = numpy.polyfit(t,T,1,full = True)
    return mc, res

# calculates the rsquared of a set of values given y and mean
def r2(y,yhat):
    y2 = []
    yhat2 = []
    yavg = numpy.mean(y)
    y2[:] = [(y[x]-yavg)**2 for x in range(len(y))]
    yhat2[:] = [(yhat[x]-yavg)**2 for x in range(len(yhat))]
    sst = sum(y2)
    ssr = sum(yhat2)
    rsq = ssr/sst
    return rsq

def U(tin,Tin,Tamb):
    lnTfit =[]
    Tmax, iTmax = findmax(Tin)
    iStart      = round((len(Tin) - iTmax)/2) + iTmax
    fitt        = tin[iStart:]
    Trawfit     = Tin[iStart:]
    Trawfit[:]  = [Trawfit[x]-Tamb[x] for x in range(len(Trawfit))]
    Tln         = numpy.log(Trawfit)
    eq, R       = findslope(fitt,Tln)
    lnTfit[:]   = [fitt[x]*eq[0]+eq[1] for x in range(len(fitt))]
    rsq         = r2(Tln,lnTfit)
    return eq[0], lnTfit, rsq
    
def Simp(U,T,Tamb,t):
    Tcorr, Tdif = [], []
    Tdif[:]  = [T[x]-Tamb[x] for x in range(len(T))]
    for k in range(len(T)):
        Tadj = 0
        if k > 2:
            Tadj = -U*integrate.simps(Tdif[0:k],t[0:k])
        Tcorr.append(T[k]+Tadj)
    return Tcorr

def mvAvg(Val,pt):
    num = int((pt-1)/2)
    Val2 = []
    for k in range(len(Val)):
        if k > num and k < ((len(Val)-1)-num):
            Val2.append(numpy.mean(Val[(k-num):(k+num)]))
        else:
            Val2.append(Val[k])
    return Val2
        
def dxdt(Val,t):
    dvaldt = []
    for k in range(len(Val)):
        if k > 0:
            dvaldt.append((Val[k]-Val[k-1])/(t[k]-t[k-1]))
        else:
            dvaldt.append(0)
    return dvaldt
            

#Glob finds all csvs and opens them
for file in glob.glob("*.csv"):
    
    print(file)
    
#Raw data arrays, used to hold data from thee CSV   
    t = []
    T1raw, T2raw, T3raw, Tjraw = [] , [], [], []
    mraw = []
    Tdif = []
    Hraw= []
    
#Reads the data from the current CSV and stores data in arrays
    with open(file,'rt') as csvfile:
        csvRead = csv.reader(csvfile,delimiter=',')
        csvfile.readline()
        for row in csvRead:
            t.append(float(row[0]))
            T1raw.append(float(row[1]))
            T2raw.append(float(row[2]))
            T3raw.append(float(row[3]))
            Tjraw.append(float(row[4]))
            Hraw.append(float(row[5]))
            mraw.append(float(row[6]))
                
#Get the slopes of ln(T-Tamb), U heat transfer coefficient
# fit as well as actual fit and rsq
    U1, Tfit1, rsq1 = U(t, T1raw, Tjraw)
    U2, Tfit2, rsq2 = U(t, T2raw, Tjraw)
    U3, Tfit3, rsq3 = U(t, T3raw, Tjraw)
    
    #print (rsq1, rsq2, rsq3)
#Use Simpsons law to integrate and correct for heat loss in temperature
    
    T1corr = Simp(U1,T1raw,Tjraw,t)
    T2corr = Simp(U2,T2raw,Tjraw,t)
    T3corr = Simp(U3,T3raw,Tjraw,t)
    
#fix the initial height measurement by fitting straight line
    dH1dt = []
    for l in range(len(Hraw)):
        if l > 0:
            dH1dt.append(((Hraw[l]-Hraw[l-1])/(t[l]-t[l-1]))**2)
        else:
            dH1dt.append(0)
    lim = 5000
    
    for k in range(len(dH1dt)):
        exc = [i for i,x in enumerate(dH1dt) if x > lim]
    
    mfit = (Hraw[exc[-1]+2]-Hraw[exc[0]-2])/(t[exc[-1]+2]-t[exc[0]-2])
    cfit = Hraw[exc[0]-2] - mfit * t[exc[0]-2]
    
    for z in range(exc[0]-2,exc[-1]+2):
        Hraw[z] = mfit*t[z]+cfit
    
    #Fix any anomalous readings in the mass reading
    
    k = 0
    mcorr = []
    for k in range(len(mraw)):
        if k >2 :
            maxm = max(mcorr[0:(k-1)])

            if any([mraw[k] < (0.5*maxm),
                    mraw[k] < 0.8*mraw[k-1],
                    mraw[k] > 2*mraw[k-1]]):
                mcorr.append(mcorr[k-1])
            else:
                mcorr.append(mraw[k])
        else:
            mcorr.append(mraw[k])
        
        
#Apply a 21 point moving average to the outputs
    T1corr  = mvAvg(T1corr,21)
    T2corr  = mvAvg(T2corr,21)
    T3corr  = mvAvg(T3corr,21)
    Hcorr   = mvAvg(Hraw,21)
    mcorr   = mvAvg(mcorr,21)
    
#Get dx/dt for each of the responses
    dT1dt   = dxdt(T1corr,t)
    dT2dt   = dxdt(T2corr,t)
    dT3dt   = dxdt(T3corr,t)
    dHdt    = dxdt(Hcorr,t)
    dmdt    = dxdt(mcorr,t)
    
#Figures out the File name and where to split
    name = file.split('.')
    if len(name) == 3:
        nameout = name[0]+name[1]
    else: nameout = name[0]
    
#Writing all the Data to a CSV
    CSVout = nameout+ "_Corrected.csv"    
    with open(CSVout,'w', newline='') as f: 
        head = ['t', 'T1', 'T2', 'T3', 'Tj', 'H', 'm', 'T1corr', 'T2corr',
                'T3corr', 'Hcorr', 'mcorr', 'dT1dt', 'dT2dt', 'dT3dt',
                'dHdt', 'dmdt']
        csvwrite = csv.DictWriter(f,fieldnames = head)
        csvwrite.writeheader()
        for j in range(len(t)):
            csvwrite.writerow({'t': t[j], 'T1': T1raw[j], 'T2': T2raw[j], 'T3': T3raw[j],
            'Tj': Tjraw[j],'H': Hraw[j], 'm': mraw[j], 'T1corr': T1corr[j], 
            'T2corr': T2corr[j], 'T3corr': T3corr[j], 'Hcorr': Hcorr[j],
            'mcorr': mcorr[j], 'dT1dt': dT1dt[j], 'dT2dt': dT2dt[j],
            'dT3dt': dT3dt[j], 'dHdt': dHdt[j], 'dmdt': dmdt[j]})
    
    
#Code for plotting the data if needed
    pltT = 'Y'
    if pltT == 'Y' :
        
        fig = plt.figure(figsize=(10,12))
        
        ax1 = fig.add_subplot(321)
        ax2 = fig.add_subplot(322)
        ax3 = fig.add_subplot(323)
        ax4 = fig.add_subplot(324)
        ax5 = fig.add_subplot(325)
        ax6 = fig.add_subplot(326)
        
        ax1.plot(t,T1corr)
        ax1.plot(t,T1raw)
        ax1.set_xlabel('time (s)',fontsize =20)
        ax1.set_ylabel('T1 /$^\circ$C',fontsize =20)
        ax1.set_xlim(0,800)
        ax1.set_ylim(bottom = 0)
    
        ax3.plot(t,T2corr)
        ax3.plot(t,T2raw)
        ax3.set_xlabel('time (s)',fontsize =20)
        ax3.set_ylabel('T2 /$^\circ$C',fontsize =20)
        ax3.set_xlim(0,800)
        ax3.set_ylim(bottom = 0)
    
        ax5.plot(t,T3corr)
        ax5.plot(t,T3raw)
        ax5.set_xlabel('time (s)',fontsize =20)
        ax5.set_ylabel('T3 /$^\circ$C',fontsize =20)
        ax5.set_xlim(0,800)
        ax5.set_ylim(bottom = 0)
        
        ax2.plot(t,Hcorr)
        ax2.plot(t,Hraw)
        ax2.set_xlabel('time (s)',fontsize =20)
        ax2.set_ylabel('H /mm',fontsize =20)
        ax2.set_xlim(0,800)
        ax2.set_ylim(bottom = 0)
        
        ax4.plot(t,mcorr)
        ax4.plot(t,mraw)
        ax4.set_xlabel('time (s)',fontsize =20)
        ax4.set_ylabel('m /g',fontsize =20)
        ax4.set_xlim(0,800)
        ax4.set_ylim(bottom = 0)
        
        #ax6.plot(t,dHdt)
        ax6.plot(t,dmdt)
        ax6.plot(t,dT2dt)
        ax6.set_xlabel('time (s)',fontsize =20)
        ax6.set_ylabel('test',fontsize =20)
        ax6.set_xlim(0,800)
        #ax6.set_ylim(bottom = 0)
        
        
        plt.plot()
    