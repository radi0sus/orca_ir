# -*- coding: utf-8 -*-
import sys                              #sys files processing
import os                               #os file processing
import re                               #regular expresions
import numpy as np                      #summation
import matplotlib.pyplot as plt         #plots
import argparse                         #argument parser
from scipy.signal import find_peaks     #peak detection

# global constants
found_ir_section=False                  #check for IR data in out
specstring_start='IR SPECTRUM'          #check orca.out from here
specstring_end='The first'              #stop reading orca.out from here
w = 15                                  #w = line width, FWHM
wn_add = 250                            #add +250 to spectra x (required for convolution)
# plot config section
high_to_low_wn = True                   #go from high to low wavenumber, normal for IR spectra, low wn to high wn if False
transm_style = True                     #show spectra in transmition style, normal for IR spectra, absorption style if False
show_grid = True                        #show grid if True
show_single_gauss = False               #show single gauss functions if True
show_single_gauss_area = False          #show single gauss functions - area plot if True
label_rel_pos_y = -15                   #-15 for transmission style, 5 for absortion style 
save_spectrum = True                    #save spectrum if True
show_spectrum = True                    #show the matplotlib window if True
label_peaks = True                      #show peak labels if True
minor_ticks = True                      #show minor ticks if True
spectrum_title = "IR spectrum"          #title
spectrum_title_weight = "bold"          #weight of the title font: 'normal' | 'bold' | 'heavy' | 'light' | 'ultrabold' | 'ultralight'
y_label = "Intensity"                   #label of the y axis
x_label = r'$\tilde{\nu}$ /cm$^{-1}$'   #label of the x axis

#global lists
modelist=list()         #mode
freqlist=list()         #frequency
intenslist=list()       #intensity absolute
gauss_sum=list()        #list for to sum of single gaussian spectra = the convoluted spectra

def gauss(a,m,x,w):
    # calculation of the Gaussian line shape
    # a = amplitude (max y, intensity)
    # x = position
    # m = maximum/meadian (stick position in x, wavenumber)
    # w = line width, FWHM
    return a*np.exp(-(np.log(2)*((m-x)/w)**2))

# parse arguments
parser = argparse.ArgumentParser(prog='orca_ir', description='Easily plot IR spectra from orca.out')
parser.add_argument("filename", help="the ORCA output file")
args = parser.parse_args()


#open a file
#check existence
try:
    with open(args.filename, "r") as input_file:
        for line in input_file:
            #start exctract text 
            if line.startswith(specstring_start):
                #found IR data in orca.out
                found_ir_section=True
                for line in input_file:
                    #stop exctract text 
                    if line.startswith(specstring_end):
                        break
                    #only recognize lines that start with number
                    #split line into 3 lists mode, frequencies, intensities
                    #line should start with a number
                    if re.search("\d:",line): 
                        modelist.append(int(line.strip().split(":")[0])) 
                        freqlist.append(float(line.strip().split()[1]))
                        intenslist.append(float(line.strip().split()[2]))
#file not found -> exit here
except IOError:
    print(f"'{args.filename}'" + " not found")
    sys.exit(1)

#no IR data in orca.out -> exit here
if found_ir_section == False:
    print(f"'{specstring_start}'" + "not found in" + f"'{args.filename}'")
    sys.exit(1)

#prepare plot
fig, ax = plt.subplots()

#plotrange in x
plt_range_x=np.arange(0,max(freqlist)+wn_add,1)

#plot single gauss function for every frequency freq
#generate summation of single gauss functions
for index, freq in enumerate(freqlist):
    #single gauss function line plot
    if show_single_gauss:
        ax.plot(plt_range_x,gauss(intenslist[index], plt_range_x, freq, w),color="black",alpha=0.5) 
    #single gauss function filled plot
    if show_single_gauss_area:
        ax.fill_between(plt_range_x,gauss(intenslist[index], plt_range_x, freq, w), color="grey",alpha=0.3)
    # sum of gauss functions
    gauss_sum.append(gauss(intenslist[index], plt_range_x, freq, w))

#y values of the gauss summation
plt_range_gauss_sum_y = np.sum(gauss_sum,axis=0)

#find peaks scipy function, change height for level of detection
peaks , _ = find_peaks(plt_range_gauss_sum_y, height = 0.5)

#plot spectra
ax.plot(plt_range_x,plt_range_gauss_sum_y,color="black",linewidth=0.8)
#plot sticks
ax.stem(freqlist,intenslist,"dimgrey",markerfmt=" ",basefmt=" ")
#optional mark peaks - uncomment in case
#ax.plot(peaks,plt_range_gauss_sum_y[peaks],"|")

#label peaks
if label_peaks:
    if not transm_style:
        label_rel_pos_y=5   #in case of absorption style spectra
        
    for index, txt in enumerate(peaks):
        ax.annotate(peaks[index],xy=(peaks[index],plt_range_gauss_sum_y[peaks[index]]),ha="center",rotation=90,size=6,
            xytext=(0,label_rel_pos_y), textcoords='offset points')
    
ax.set_xlabel(x_label)                                          #label x axis
ax.set_ylabel(y_label)                                          #label y axis
ax.set_title(spectrum_title,fontweight=spectrum_title_weight)   #title
ax.get_yaxis().set_ticks([])                                    #remove ticks from y axis
plt.tight_layout()                                              #tight layout
if minor_ticks:
    ax.minorticks_on()                                          #show minor ticks

#x,y axis manipulation
if transm_style:
    plt.ylim(max(plt_range_gauss_sum_y)+max(plt_range_gauss_sum_y)*0.1,0) # +10% for labels
else:
    plt.ylim(0,max(plt_range_gauss_sum_y)+max(plt_range_gauss_sum_y)*0.1) # +10% for labels
    
if high_to_low_wn:
    plt.xlim(max(plt_range_x),0)
else:
    plt.xlim(0,max(plt_range_x))
    
#plt.xlim(max(plt_range_x),0)
if show_grid:
    ax.grid(True,which='major', axis='x',color='black',linestyle='dotted', linewidth=0.5)


#increase figure size N times
N = 1.5
params = plt.gcf()
plSize = params.get_size_inches()
params.set_size_inches((plSize[0]*N, plSize[1]*N))

#save the plot
if save_spectrum:
    filename, file_extension = os.path.splitext(args.filename)
    plt.savefig(f"{filename}-ir.png", dpi=300)

#show the plot
if show_spectrum:
    plt.show()
