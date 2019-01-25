# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 18:18:04 2017

@author: chp16hcw
"""


import csv, glob,numpy, os

os.chdir("C:\\Users\\chp16hcw\\Google Drive\\PhD\\Results\\KTest\\ATR\\Corrected\\")


filelist = glob.glob("*.csv")

CSVout = 'ATRsumm.csv'
with open(CSVout,'w', newline='') as f: 
    head = ['Sample', 'DelT1max', 'DelT2max', 'DelT3max',
            'time Max T1','time Max T2','time Max T3',
            'Atime Max T1','Atime Max T2','Atime Max T3',
            'Max H', 'H Final', 'time Max H', '% Sag', 
            'm max','m Final', 'm Loss', 'Hmax/m','Hf/m',
            'Sag/m', 'Sag % / m' ]
    csvwrite = csv.DictWriter(f,fieldnames = head)
    csvwrite.writeheader()


for file in filelist:
    
    time, DelT1, DelT2, DelT3, H, m, Tj = [], [], [], [], [], [], []
    T1, T2, T3  = [], [], []
#grabs data from the current csv   
    with open(file,'rt') as csvfile:
        csvRead = csv.reader(csvfile,delimiter=',')
        csvfile.readline()
        for row in csvRead:
            time.append(float(row[0]))
            T1.append(float(row[1]))
            T2.append(float(row[2]))
            T3.append(float(row[3]))
            DelT1.append(float(row[7]))
            DelT2.append(float(row[8]))
            DelT3.append(float(row[9]))
            Tj.append(float(row[4]))
            H.append(float(row[10]))
            m.append(float(row[11]))
    
    sample = file.split('_')[0]
    print (sample)
    DelTmax1 = numpy.amax(DelT1)-numpy.average(Tj)
    DelTmax2 = numpy.amax(DelT2)-numpy.average(Tj)
    DelTmax3 = numpy.amax(DelT3)-numpy.average(Tj)
    T1max    = numpy.amax(T1) 
    T2max    = numpy.amax(T2) 
    T3max    = numpy.amax(T3) 
    Hmax     = numpy.amax(H)
    Hfin     = numpy.average(H[-200:-50])
    Sagper   = 100 * ((Hmax-Hfin) / Hmax)
    mmax     = numpy.amax(m)
    mfin     = numpy.average(m[-200:-50])
    mloss    = mmax-mfin
    Hmaxm    = Hmax/mfin
    Hfinm    = Hfin/mfin
    sagm     = ((Hmax-Hfin) / Hmax)/mfin
    sagperm  = Sagper/mfin

    k = [i for i,x in enumerate(H) if x > 0.95* Hmax][0]
    l = [i for i,x in enumerate(DelT1) if x > 0.95*numpy.amax(DelT1)][0] 
    m = [i for i,x in enumerate(DelT2) if x > 0.95*numpy.amax(DelT2)][0]
    n = [i for i,x in enumerate(DelT3) if x > 0.95*numpy.amax(DelT3)][0]
    o = [i for i,x in enumerate(T1) if x > 0.95*numpy.amax(T1)][0]
    p = [i for i,x in enumerate(T2) if x > 0.95*numpy.amax(T2)][0]
    q = [i for i,x in enumerate(T3) if x > 0.95*numpy.amax(T3)][0]
    
    tH  = time[k]

    taT1 = time[l]
    taT2 = time[m]
    taT3 = time[n]
    tT1 = time[o]
    tT2 = time[p]
    tT3 = time[q]
    
      
    with open(CSVout,'a', newline='') as f:
        head = ['Sample', 'DelT1max', 'DelT2max', 'DelT3max',
        'time Max T1','time Max T2','time Max T3',
        'Atime Max T1','Atime Max T2','Atime Max T3',
        'Max H', 'H Final', 'time Max H', '% Sag', 
        'm max','m Final', 'm Loss', 'Hmax/m','Hf/m',
        'Sag/m', 'Sag % / m' ]
        csvwrite = csv.DictWriter(f,fieldnames = head)
        csvwrite.writerow({'Sample': sample, 'DelT1max': DelTmax1,
                           'DelT2max': DelTmax2, 'DelT3max': DelTmax3,
                           'time Max T1': tT1,'Atime Max T1':taT1,
                           'Atime Max T2':taT2,'Atime Max T3':taT3,
                           'time Max T1': tT1, 'time Max T2': tT2, 
                           'time Max T3': tT3, 'Max H': Hmax, 'H Final': Hfin,
                           'time Max H': tH, '% Sag': Sagper, 'm max': mmax, 
                           'm Final': mfin, 'm Loss': mloss, 'Hmax/m': Hmaxm, 
                           'Hf/m': Hfinm, 'Sag/m': sagm, 'Sag % / m': sagperm})