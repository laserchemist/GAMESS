#!/usr/bin/env python
# JMS 04 December 2018
# 1.4 JMS 25 January 2019
#     Python script to read GAMESS MD trajectory from output file and save as xmol
#     and process data. Beta version
#     Need to add retrieval of solvent molecules
# 1.5 JMS 30 January
#     Updated to read and process version dependent output
# Some debug printing remains, treat as beta
import re
import numpy as np
import scipy
from scipy import signal
import matplotlib.pyplot as plt
import math
from scipy import fftpack
import os
import sys

class gamess(object): # Class to extract trajectory from GAMESS file
    """Class to extract data from GAMESS output"""
    # For GAMESS documentation see: https://www.msg.chem.iastate.edu/gamess/
    def __init__(self,file):
        self.description="Class to parse GAMESS trajectory output"
        self.file=file
    def version(self):
    #     *         GAMESS VERSION = 30 SEP 2018 (R3)         *
        head="GAMESS VERSION"
        index=-1
        vfind=[]
        with open( self.file, "r" ) as source:
            for line in source:
                    if line.find( head )>0:
                        index+=1
                        match = re.search(r"([\d]{1,2})\s(\w+)\s([\d]{4})", line)
                        vlist=[match.group(1),match.group(2),match.group(3)]
                        vfind.append(vlist)
                        break
        self.version=vlist
        return vlist
    def xyz(self):
        #print("In calling function",self.file)
        with open( self.file, "r" ) as source:
            coord=coordinate(source," QM ATOM COORDINATES","CARTESIAN COORDINATES (ANG)")
            col2=[] # Key to creating an empty List to fill with x,y,z coordinate not [[]]!
            index=-1
            for i in coord:
                index+=1
                col2.append([])
                for j in i:
                   col2[index].append([float(i) for i in j[2:5]]) # Convert to float element by element
        coordxyz=col2
        self.xyz=coordxyz
        return coordxyz
    def atom(self):
        #print("In calling function",self.file)
        with open( self.file, "r" ) as source:
            coord=coordinate(source," QM ATOM COORDINATES","CARTESIAN COORDINATES (ANG)")
        col2=[] # Key to creating an empty 2D List to fill with x,y,z coordinate
        index=-1
        for i in coord:
            index+=1
            col2.append([])
            for j in i:
               col2[index].append(j[0]) # Grab element label
        self.atom=col2
        return col2
    def mass(self):
        #print("In calling function",self.file)
        with open( self.file, "r" ) as source:
            coord=coordinate(source," QM ATOM COORDINATES","CARTESIAN COORDINATES (ANG)")
            col2=[] # Key to creating an empty 2D List to fill with x,y,z coordinate
            index=-1
            for i in coord:
                index+=1
                col2.append([])
                for j in i:
                    col2[index].append(float(j[1])) # Grab element label
        self.mass=col2
        return col2
    def time(self): # To retireve time step data
        head="FSEC"
        index=-1
        timefind=[]
        with open( self.file, "r" ) as source:
            for line in source:
                for line in source:
                    if line.find( head )>0:
                        index+=1
                        #  *** AT T=          7.62 FSEC, THIS RUN'S STEP NO.=     762
                        match = re.search(r"(\d+\.\d+)\sFSEC\,\D+(\d+)", line)
                        timelist=[match.group(1),match.group(2)]
                        timefind.append(timelist)
                        break
        self.time=timefind
        return timefind
    def dipole(self): # To retireve time step data
        head="FSEC"
        dipfind="(A.U.)"
        index=-1
        dipoleline=[]
        with open( self.file, "r" ) as source:
            for line in source:
                if line.find( head )>0:
                    index+=1
                    dipoleline.append([])
                    for line in source:
                         if line.find(dipfind)>0:
                             match = re.search(r"([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)", line)
                             dplist=[match.group(1),match.group(2),match.group(3),match.group(4)]
                             dipoleline[index].append(dplist)
                             break
        self.dipole=dipoleline
        return dipoleline
    def dipnew(self): # To retireve time step data
        # DIPOLE      -0.028421  0.034413  0.047308
        dipfind="DIPOLE"
        index=-1
        dipoleline=[]
        with open( self.file, "r" ) as source:
            for line in source:
                dipoleline.append([])
                for line in source:
                        if line.find(dipfind)>0:
                             match = re.search(r"([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)", line)
                             dplist=[match.group(1),match.group(2),match.group(3)]
                             dipoleline[index].append(dplist)
                             break
        self.dipnew=dipoleline
        return dipoleline
    def Epot(self): # To retireve time step data
        head="FSEC"
        find="POT  ENERGY"
        index=-1
        capline=[]
        with open( self.file, "r" ) as source:
            for line in source:
                if line.find( head )>0:
                    index+=1
                    capline.append([])
                    for line in source:
                         if line.find(find)>0:
                             #print("In Epot: ",line)
                             match = re.search(r"([+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))", line)
                             dplist=[match.group(1)]
                             capline[index].append(dplist)
                             break
        self.Epot=capline
        return capline
    def Etot(self): # To retireve time step data
        head="FSEC"
        # Variant 1:      TOT  ENERGY      =     -4.851316117E+04 KCAL/MOL    0.00000000
        # Variant 2: TOT. HAM. ENERGY =     -1.402778943E+05 KCAL/MOL   -0.05032210
        regex = r"(?=.*\bTOT)(?=.*\bENERGY\b).*"
        p = re.compile(regex)
        index=-1
        capline=[]
        with open( self.file, "r" ) as source:
            for line in source:
                if line.find( head )>0:
                    index+=1
                    capline.append([])
                    for line in source:
                        test=line
                        if re.match(p, line) is not None:
                            if "ENERGY" in test:
                                #print("In Etot: ",line)
                                match = re.search(r"([+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))", line)
                                dplist=[match.group(1)]
                                capline[index].append(dplist)
                            break
        self.Etot=capline
        return capline
    def temp(self): # To retireve time step data
        head="FSEC"
        find="TEMPERATURE"
        index=-1
        capline=[]
        with open( self.file, "r" ) as source:
            for line in source:
                if line.find( head )>0:
                    index+=1
                    capline.append([])
                    for line in source:
                         if line.find(find)>0:
                             match = re.search(r"(\d+\.\d+)", line)
                             dplist=[match.group(1)]
                             capline[index].append(dplist)
                             break
        self.temp=capline
        return capline
    def dipderiv(self): # To retireve time step data
        head="FSEC"
        dipfind="(DEBYE)"
        index=-1
        dipoleline=[]
        with open( self.file, "r" ) as source:
            for line in source:
                if line.find( head )>0:
                    index+=1
                    dipoleline.append([])
                    for line in source:
                         if line.find(dipfind)>0:
                             for line in source:
                                 match = re.search(r"([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)", line)
                                 #print("match = ",match)
                                 if match:
                                     #print("MATCH++++++",index,match.group(2))
                                     dplist=[match.group(1),match.group(2),match.group(3),match.group(4)]
                                     dipoleline[index].append(dplist)
                                     break
                             break
        self.dipderiv=dipoleline
        return dipoleline

# ++++++++++++++++++ END of Class gamess +++++++++++++

# Need function to return array of COORDINATES
def coordinate(source_file,head,foot):
    buffer= []
    coord=[]
    index=-1
    atom=-1
    for line in source_file:
        if line.startswith( head ):
            index+=1
            coord.append([])
            #print("line: ",line)
            for line in source_file:
              if re.search(r"([A-Z]{1,2})\W+(\d+\.\d+)\W+([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)\W+([-\s]\d+\.\d+)", line):
                 coord[index].append(line.rstrip('\n').split())
                 #print("Coordinate: ",index,coord[index])
              else:
                 break
    return coord
def calc_derivative(array_1D, delta_t):
    ''' The derivatives of the an array by using the
    finite differences method.
    '''
    dy = np.gradient(array_1D)
    return np.divide(dy, delta_t)
def autocorr(x):
    result = np.correlate(x, x, mode='full')
    return result[int(result.size/2):] # Something wrong here
def zero_padding(sample_data):
    """ A series of Zeros will be padded to the end of the dipole moment
    array (before FFT performed), in order to obtain a array with the
    length which is the "next power of two" of numbers.

    #### Next power of two is calculated as: 2**math.ceil(math.log(x,2))
    #### or Nfft = 2**int(math.log(len(data_array)*2-1, 2))
    """
    return int(2 ** math.ceil(math.log(len(sample_data), 2)))
def choose_window(data, name,width):
    kind = name
    print("KIND:",kind)
    if kind == "Gaussian":
        sigma = 2 * math.sqrt(2 * math.log(2))
        std = float(width)
        window_function = signal.gaussian(len(data), std/sigma, sym=False)

    elif kind == "Blackman-Harris":
        window_function = signal.blackmanharris(len(data), sym=False)

    elif kind == "Hamming":
        window_function = signal.hamming(len(data), sym=False)

    elif kind == "Hann":
        window_function = signal.hann(len(data), sym=False)

    print(window_function)
    return window_function
def calc_FFT(array_1D, window):
    """
    This function is for calculating the "intensity" of the ACF at each frequency
    by using the discrete fast Fourier transform.
    """
####
#### http://stackoverflow.com/questions/20165193/fft-normalization
####
    #window = choose_window(array_1D, "Gaussian")
    # swindow=np.sum(window)
    # WE = swindow /array_1D[0].shape
    # print("Window shape: ",WE.shape)
    #wf=np.true_divide(window, WE)

    WE = np.sum(window) / len(array_1D)
    wf = window / WE
    # convolve the blackman-harris or other window function.
    sig = array_1D * wf
    # A series of number of zeros will be padded to the end of the \
    # VACF array before FFT.
    N = zero_padding(sig)
    # Tried using Numpy FFT but fftpack works better for this application
    #yfft = np.fft.fft(sig, N, axis=0) / len(sig)
    # yfft = np.fft.fft(sig, N, axis=0)/len(sig) # 
    # Try this... Works better, somehow above shifts spectrum to much higher cm-1
    yfft=fftpack.fft(sig)/len(sig)
    print("shape of yfft {:}".format(np.shape(yfft)))
    #return np.square(np.absolute(yfft))
    return np.abs(yfft)
#++++++++++++++++++++++++++++++++++++++++++++
# MAIN PROGRAM
#++++++++++++++++++++++++++++++++++++++++++++
# Files
if len(sys.argv)>1:
    file=sys.argv[1]
else:
    # Sample local files
    file='output65'
    file='imidazole_30w_output65'
    file='HCl_output77'
    file='imidazole_water_1fs_output78'
    file='acetylene_300K_MD_output84'
    file='acetylene_MD_output84'
    file='acetylene_MD_300K_long_output84'
    file='Acetyleneoutput91'
    file='acetylene_output91'
    path='/Users/jms/local/molecules/imidazole/'
    file='output65'
    file='imidazole_30w_output65'
    file='HCl_output77'
    file='imidazole_water_1fs_output78'
    file='acetylene_300K_MD_output84'
    file='acetylene_MD_output84'
    file='acetylene_MD_300K_long_output84'
    file='Acetyleneoutput91'
    file='acetylene_output91'
    file='acetylene_300K_MD_output84'
    print("No command line argument using default file: ", file)
if len(sys.argv)>2:
    path=sys.argv[2]
else:
    path='/Users/jms/local/molecules/imidazole/'
    print("No command line argument using default path: ", path)
if len(sys.argv)>3:
    fend=sys.argv[3]
else:
    fend=".txt"
filepath=path+file+fend
filepathdat=path+file+".dat"
print("\n\nGAMESS_readv1.5.py")
print("++++++++++++++++++")
print("ver. 1.5 JMS 31 January 2019\n\n")
print("File: ",filepath)
outlines=[]
numline=[]
substr='fsec'
print("COORDINATES:")
calc1=gamess(filepath)
calcdat=gamess(filepathdat)
print(calc1.description)
calc2=gamess(filepath)
print(calc1.file)
coords=calc1.xyz()
print("Length:",len(coords[0]))
print("Length overall:",len(coords))
# Write xmol file to analyze trajectory in VMD, etc.
#
natoms=len(coords[0])
steps=len(coords)
print("+ Coordinates retrieved dimension: rows,cols ", steps,natoms)
mass=calc1.mass()
filename=filepath+".xmol"
print("Writing file to:",filename)
with open(filename, 'w') as f:
    for item in coords:
        f.write(str(natoms)+"\n")
        f.write("Output from GAMESS\n")
        atom=0
        #print("coordinate line: ",item)
        for xyz in item:
            strxyz=" ".join([str(i) for i in xyz])
            textline=str(int(mass[0][atom]))+" "+strxyz
            #print(textline)
            f.write(textline)
            f.write("\n")
            atom+=1
#Now retrieve dipole moment as a function of time from trajectory
time=calc1.time()
print("Last Time: ", time[-1][0])
version=calc1.version()
print(version,version[2])
if float(version[2])<2018:
        #Dipole
        dipole=calc1.dipole()
        print("dipole length ",len(dipole))
        dipderiv=calc1.dipderiv()
        print("dipole derivative length ",len(dipderiv))
        print("dipole derivative: ", dipderiv[0],"++++++++++") # Works
        print("Lengths ", len(dipderiv[0][:]),len(dipderiv[:]),len(dipderiv[:][:]),len(dipderiv[3][:]))
        ndderiv=np.array(dipderiv[0:-1],dtype=float) # ValueError: setting an array element with a sequence.
        print("Shape of numpy array: ",ndderiv.shape)
        nddipole=np.array(ndderiv,dtype=float)[:,0,:] #Goes from 3D to 2D
        print("Shape of nddipole numpy array: ",nddipole.shape)
else:
        print("Newer version ") # 2018 version puts dipole this in .dat file
        # Capture this text:  DIPOLE      -0.028421  0.034413  0.047308
        dipole=calcdat.dipnew()
        print("dipole length ",len(dipole))
        #nddipole=np.array(dipole[0:-1],dtype=float) # ValueError: setting an array element with a sequence.
        #print("Shape of numpy array: ",nddipole.shape)
        nddipole=np.array(nddipole,dtype=float)[1:,0,:] #Goes from 3D to 2D REMOVE first row to match dimensions
        print("Shape of nddipole numpy array: ",nddipole.shape)
        time_step=float(time[-2][1])-float(time[-3][1]) # 1.0 # Time step in FSEC [22 January]
        ddipolez=calc_derivative(nddipole[:,2], time_step)
time_step=float(time[-1][0])-float(time[-2][0]) # 1.0 # Time step in FSEC [22 January]
print("Time step [fsec]: ",time_step)
print("Initial time step retrieved: ",time[0][0])
print("Final time step retrieved: ",time[-1][0])
print("Time step retrieved: ",time_step)
ddipolez=calc_derivative(nddipole[:,2], time_step)
ddipoley=calc_derivative(nddipole[:,1], time_step)
ddipolex=calc_derivative(nddipole[:,0], time_step)
nxyz=np.asfarray(coords) # Create numpy array to carry out advanced signal processing
print("Shape of numpy array: ",nxyz.shape)
zderiv=ddipolez
az=autocorr(ddipolez)
ay=autocorr(ddipoley)
ax=autocorr(ddipolex)
sum=az #+ay+ax # Using SUM of autocorrelated dipole derivatives
psum=zero_padding(sum)
print("+ Zero padding for FFT: ",psum)
sigma = 2 * math.sqrt(2 * math.log(2))
standard=4000
std = float(standard)
wf_name="Blackman-Harris"
# "The larger the number, the narrower the line-shape.
# "(e.g. 4000 yields FWHM around 10 cm^-1 for a 20-ps traj.)",
width=2000
timelength=float(time[-1][0])
width=4000/900*timelength
window = choose_window(sum, wf_name,width)
N = zero_padding(sum)
#yfft = np.fft.fft(sig, N, axis=0) / len(sig)
testfft = np.fft.fft(sum, N, axis=0) # no window func
plt.figure(figsize=(6, 5))
plt.plot(testfft)
plt.title('Test FFT')
dpfft=calc_FFT(sum, window)
plt.figure(figsize=(6, 5))
plt.plot(dpfft)
plt.title('Test FFT with window')
# TRY fftpack
#multiply by 27 THz signal
time_step=float((float(time[-1][0])-float(time[0][0]))/len(sum)) # 1.0 # Time step in FSEC
print("time_step",time_step)
end_time=float(time[-2][0])+time_step
start_time=float(time[0][0])
time_vec = np.arange(start_time,end_time, time_step) # Time step vector in FSEC
f=1/time_step
f=0.038 # 1260 cm-1 vibration
# 1 femtosecond = 1000 THz
sigsim=(np.sin(f * 2 * np.pi * time_vec)+1.0)
f=0.0007
sigsim2=(np.sin(f * 2 * np.pi * time_vec)+1.0)
# 1 femtosecond = 1000 THz
plt.figure(figsize=(6, 5))
plt.plot(time_vec, sum, label='Original signal')
#+sigsim #*sigsim*sigsim2
#plt.plot(time_vec, sum, label='Multiplied signal')
plt.legend()
sig_fft = fftpack.fft(sum)
# And the power (sig_fft is of complex dtype)
# FT ACF
# +++++++++++++++++++++++++++++++++++++++++++
power = np.abs(sig_fft)[0:int(sum.size / 2)]
wf_name="Blackman-Harris"
window = choose_window(sum, wf_name,width)
powerBH=calc_FFT(sum, window)[0:int(sum.size / 2)]
wf_name="Hann"
window = choose_window(sum, wf_name,width)
powerHann=calc_FFT(sum, window)[0:int(sum.size / 2)]
wf_name="Hamming"
window = choose_window(sum, wf_name,width)
powerHam=calc_FFT(sum, window)[0:int(sum.size / 2)]
power2=dpfft[0:int(sum.size / 2)]
# The corresponding frequencies
c = 2.9979245899e-5 # speed of light in vacuum in [cm/FSEC]
sample_freq = fftpack.fftfreq(sum.size, d=time_step/1000)[0:int(sum.size / 2)]
sample_omega=np.fft.fftfreq(sum.size, time_step*2 * np.pi/1000)[0:int(sum.size / 2)]
wavenumber = fftpack.fftfreq(sum.size, time_step * c)[0:int(sum.size / 2)]
wavenumber1 = np.fft.fftfreq(sum.size, time_step * c)[0:int(sum.size / 2)]
kb=1.38064E-23
T=298
hbar=1.0545718E-34
# Compute Infrared intensity 
IR_intensity=power2*sample_omega*np.tanh(hbar*1.0E15*sample_omega/(kb*T))
IR_intensityBH=powerBH*sample_omega*np.tanh(hbar*1.0E15*sample_omega/(kb*T))
IR_intensityHann=powerHann*sample_omega*np.tanh(hbar*1.0E15*sample_omega/(kb*T))
IR_intensityHam=powerHam*sample_omega*np.tanh(hbar*1.0E15*sample_omega/(kb*T))
# Create a series of diagnostic plots
# Plot the FFT power
plt.figure(figsize=(6, 5))
plt.plot(sample_freq, power)
plt.xlabel('Frequency [THz]')
plt.ylabel('power')
plt.figure(figsize=(6, 5))
plt.plot(sample_omega, power)
plt.xlabel('Frequency [omega]')
plt.ylabel('power')
plt.figure(figsize=(6, 5))
plt.plot(wavenumber, power)
plt.grid(True)
plt.xlim(0, 5000)
plt.xlabel('Frequency [cm^-1]')
plt.ylabel('power')
plt.figure(figsize=(6, 5))
plt.plot(wavenumber, IR_intensityBH, '-b', label="Blackman-Harris")
plt.plot(wavenumber,IR_intensityHann,'-r', label="Hann")
plt.plot(wavenumber,IR_intensityHam,'-g', label="Hamming")
plt.xlabel('Frequency [cm^-1]')
plt.legend()
plt.ylabel('IR_intensity')
plt.grid(True)
plt.xlim(0, 7000)
plt.savefig(filename+'_IRNW'+'.png')
plt.figure(figsize=(6, 5))
plt.plot(wavenumber, IR_intensity)
plt.xlabel('Frequency [cm^-1]')
plt.ylabel('IR_intensity')
label="Window width: "+str(width)
plt.text(3000,1,label)
plt.grid(True)
plt.xlim(0, 5000)
plt.savefig(filename+'_IR'+'.png')
filename=filepath+"_spec.csv"
print("Writing spectrum file to:",filename)
# Write out Spectra and MD data
# Quick way to write out 2D array
# https://www.python-course.eu/numpy_reading_writing.php
spectrum=np.vstack((wavenumber,IR_intensity,IR_intensityBH,IR_intensityHann,IR_intensityHam)).T
np.savetxt(filename, spectrum, delimiter=',', fmt='%10.4f')
potentialE=calc1.Epot()
potE=np.asfarray(potentialE)
totalE=calc1.Etot()
totE=np.asfarray(totalE)
temp=calc1.temp()
MDtemp=np.asfarray(temp)
# These arrays get mishaped, have to fix ad hoc
print(time_vec.shape)
print(MDtemp.shape)
MDtemp=MDtemp[:,0,0]
print(MDtemp.shape)
potE=potE[:-1,0,0]
print(potE.shape)
print(totE.shape)
print(totE.ndim)
if totE.ndim>1:
    totE=totE[:,0,0]
print(totE.shape)
# Need to carefully assemble output data
outdata=np.stack((time_vec,potE,totE,MDtemp), axis=-1)
print(outdata.shape)
filename=filepath+"_dat"+".csv"
np.savetxt(filename, outdata, delimiter=',', fmt='%10.4f', newline=os.linesep)
# Find the peak frequency: we can focus on only the positive frequencies
plt.figure(figsize=(6, 5))
plt.plot(time_vec, MDtemp)
plt.xlabel('time [FSEC]')
plt.ylabel('MDtemp')
plt.figure(figsize=(6, 5))
plt.plot(time_vec, potE)
plt.xlabel('time [FSEC]')
plt.ylabel('Potential Energy')
# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(3, sharex=True)
axarr[0].plot(time_vec, MDtemp)
axarr[0].set_title('MD Analysis: MDTEMP [upper], Potential Energy [lower]')
axarr[1].plot(time_vec, potE)
axarr[2].plot(time_vec, totE)
axarr[0].set_xlabel('Time [fsec]')
plt.savefig(filename+'.png')
plt.show()
# Finish