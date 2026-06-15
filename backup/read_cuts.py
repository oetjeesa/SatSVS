import os
import numpy as np
import matplotlib.pyplot as plt

# import matplotlib
# matplotlib.use("TkAgg")
# from matplotlib import pyplot as plt
# os.environ['PROJ_LIB'] = '/Users/micheltossaint/Documents/anaconda3/share/proj'

def read_cuts(in_file):
    line_no     = 0
    values      = []
    
    with open(in_file, 'r') as grd:
        for line in grd:
            line_no +=1          
            if False: #line_no >= stop_line and line_no < (stop_line + 4):
                parameters.append(line.split())
            elif line == "\n":
                break
            else:
                values.append(line.split())    
        Xmin = float(values[1][0])
        DX = float(values[1][1])
        NX = int(values[1][2])
        Xmax = Xmin+(NX-1)*DX
        
        NPHI = int(np.size(values,0)/(NX+2))
        
        Ex=np.zeros((NX,NPHI))+ 1j*np.zeros((NX,NPHI))
        Ey=np.zeros((NX,NPHI))+ 1j*np.zeros((NX,NPHI))
        PHI=np.zeros(NPHI) 
        
        X = np.linspace(Xmin, Xmax, NX)
    
        for m in range(NPHI):
            PHI[m] = float(values[1+m*(NX+2)][3])
            for k in range(2+m*(NX+2),NX+2+m*(NX+2)):  #range(2,NX+2):
                Ex[k-2-m*(NX+2)][m]=float(values[k][0]) + 1j * float(values[k][1])
                Ey[k-2-m*(NX+2)][m]=float(values[k][2]) + 1j * float(values[k][3])
        
    return Ex, Ey, X, PHI, NPHI

filename = 'pattern18.7.cut'
Ex, Ey, X, PHI, NPHI = read_cuts(filename)

fig = plt.figure(figsize=(10, 5))
for i, cut in enumerate(PHI):
    plt.plot(X, 10*np.log10(np.absolute(Ex[:,i])), '-',label=str(cut)+' [deg] azimuth cut')
plt.legend()
plt.xlabel('Theta [deg]')
plt.ylabel('Gain co-polar [dB]')
plt.title(filename+' co-polar gain [dB]')
plt.grid()
plt.show()

fig = plt.figure(figsize=(10, 5))
for i, cut in enumerate(PHI):
    plt.plot(X, 10*np.log10(np.absolute(Ey[:,i])), '-',label=str(cut)+' [deg] azimuth cut')
plt.legend()
plt.xlabel('Theta [deg]')
plt.ylabel('Gain cross-polar [dB]')
plt.title(filename+' cross-polar gain [dB]')
plt.grid()
plt.show()