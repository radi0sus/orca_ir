# -*- coding: utf-8 -*-
'''

# orca-ir
Python 3 script for (hassle-free) plotting IR spectra from ORCA output files with 
peak dectection and annotation.
It combines the stick spectrum with the convoluted spectrum (gaussian line shape). 
The full spectrum or parts of the spectrum (via matplotlib window) can be plotted.

## External modules
 `re`, `numpy`, `matplotlib`, `scipy`  

## Quick start
 Start the script with:
`python3 orca-ir.py filename`
it will save the plot as PNG bitmap:
`filename-ir.png`


## Command-line options
- `filename` required: filename
- `-w` optional: line width of the gaussian (default is 15)
- `-s` optional: shows the `matplotlib` window
- `-n` optional: do not save the spectrum

## Script options
There are numerous ways to configure the spectrum in the script:
Check `# plot config section - configure here` in the script. 
Here, you can configure an absorption or transmittance plot for example.
You can even plot the single gaussian functions.

## Code options
Colors, line thickness, line styles, level of peak detection and 
more can be changed in the code directly.

## Special options and limitations
The spectrum always starts at zero and ends at the maximum wave number. 
If you need only a part of the spectrum, you can start the script with:
`python3 orca-ir.py filename -s`
and use the matplotlib window to zoom to an area of interest and save it.

The PNG file will be replaced everytime you start the script with the same output file. 
If you want to keep the file, you have to rename it. 

## Examples:
![Example 1](/examples/example1.png)
![Example 2](/examples/example2.png)
![Example 3](/examples/example3.png)

'''

import sys                              #sys files processing
import os                               #os file processing
import re                               #regular expressions
import numpy as np                      #summation
import matplotlib.pyplot as plt         #plots
import argparse                         #argument parser
from scipy.signal import find_peaks     #peak detection

# global constants
found_ir_section=False                  #check for IR data in out
specstring_start='IR SPECTRUM'          #check orca.out from here
specstring_end='The first'              #stop reading orca.out from here
w = 15                                  #w = line width for broadening, FWHM
wn_add = 150                            #add +150 to spectra x (required for convolution)

# plot config section - configure here
high_to_low_wn = True                   #go from high to low wave number, normal for IR spectra, low wn to high wn if False
transm_style = True                     #show spectra in transmittance style, absorption style if False
show_grid = True                        #show grid if True
show_conv_spectrum = True               #show the convoluted spectra if True (if False peak labels will not be shown)
show_sticks = True                      #show the stick spectra if True
show_single_gauss = True                #show single gauss functions if True
show_single_gauss_area = True           #show single gauss functions - area plot if True
label_rel_pos_y = -15                   #-15 for transmittance style, 5 for absorption style 
save_spectrum = True                    #save spectrum if True
show_spectrum = False                   #show the matplotlib window if True
label_peaks = True                      #show peak labels if True
minor_ticks = True                      #show minor ticks if True
linear_locator = True                   #tick locations at the beginning and end of the spectrum x-axis, evenly spaced
spectrum_title = "IR spectrum"          #title
spectrum_title_weight = "bold"          #weight of the title font: 'normal' | 'bold' | 'heavy' | 'light' | 'ultrabold' | 'ultralight'
y_label_trans = "transmittance"         #label of the y axis - Transmittance
y_label_abs = "intensity"               #label of y-axis - Absorption
x_label = r'$\tilde{\nu}$ /cm$^{-1}$'   #label of the x-axis
figure_dpi = 300                        #DPI of the picture

#global lists
modelist=list()         #mode
freqlist=list()         #frequency
intenslist=list()       #intensity absolute
gauss_sum=list()        #list for the sum of single gaussian spectra = the convoluted spectrum


def roundup(x):
    #round to next 100
    return x if x % 100 == 0 else x + 100 - x % 100

def gauss(a,m,x,w):
    # calculation of the Gaussian line shape
    # a = amplitude (max y, intensity)
    # x = position
    # m = maximum/meadian (stick position in x, wave number)
    # w = line width, FWHM
    return a*np.exp(-(np.log(2)*((m-x)/w)**2))

# parse arguments
parser = argparse.ArgumentParser(prog='orca_ir', description='Easily plot IR spectra from orca.out')

#filename is required
parser.add_argument("filename", help="the ORCA output file")

#show the matplotlib window
parser.add_argument('-s','--show',
    default=0, action='store_true',
    help='show the plot window')

#do not save the png file of the spectra
parser.add_argument('-n','--nosave',
    default=1, action='store_false',
    help='do not save the spectra')

#change line with (integer) for line broadening, default is 15
parser.add_argument('-w','--linewidth',
    type=int,
    default=15,
    help='line width for broadening')

#pare arguments
args = parser.parse_args()

#change values according to arguments
show_spectrum = args.show
save_spectrum = args.nosave   

#check if w is between 1 and 100, else reset to 15
if 1<= args.linewidth <= 100:
    w = args.linewidth
else:
    print("warning! line width exceeds range, reset to 15")
    w = 15


#open a file
#check existence
try:
    with open(args.filename, "r") as input_file:
        for line in input_file:
            #start exctract text 
            if "Program Version 5" in line:
                #thanks to the orca prgrmrs intensity now in a different column
                intens_column=3
            if "Program Version 4" in line:
                intens_column=2
            if "Program Version 3" in line:
                intens_column=2
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
                        intenslist.append(float(line.strip().split()[intens_column]))
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
        ax.plot(plt_range_x,gauss(intenslist[index], plt_range_x, freq, w),color="grey",alpha=0.5) 
    #single gauss function filled plot
    if show_single_gauss_area:
        ax.fill_between(plt_range_x,gauss(intenslist[index], plt_range_x, freq, w),color="grey",alpha=0.5)
    # sum of gauss functions
    gauss_sum.append(gauss(intenslist[index], plt_range_x, freq, w))

#y values of the gauss summation
plt_range_gauss_sum_y = np.sum(gauss_sum,axis=0)

#find peaks scipy function, change height for level of detection
peaks , _ = find_peaks(plt_range_gauss_sum_y, height = 0.1)

#plot spectra
if show_conv_spectrum:
    ax.plot(plt_range_x,plt_range_gauss_sum_y,color="black",linewidth=0.8)

#plot sticks
if show_sticks:
    ax.stem(freqlist,intenslist,"dimgrey",markerfmt=" ",basefmt=" ")

#optional mark peaks - uncomment in case
#ax.plot(peaks,plt_range_gauss_sum_y[peaks],"|")

#label peaks
#show peak labels only if the convoluted spectrum is shown (first if)
if show_conv_spectrum:  
    if label_peaks:
        if not transm_style:
            label_rel_pos_y=5   #in case of absorption style spectra
            
        for index, txt in enumerate(peaks):
            
            if transm_style:
                #corr_factor - maintain distance from label to peak in transmittance style
                #sensitive to peak label font size
                corr_factor = (4-len(str(peaks[index])))*3.75
                #if one does not care:
                #corr_factor = 0 
            else:
                # not necessary for absorption style
                corr_factor = 0
                
            ax.annotate(peaks[index],xy=(peaks[index],plt_range_gauss_sum_y[peaks[index]]),ha="center",rotation=90,size=6,
                xytext=(0,label_rel_pos_y+corr_factor), textcoords='offset points')
    
ax.set_xlabel(x_label)                                          #label x axis

if transm_style:
    ax.set_ylabel(y_label_trans)                                #label y axis
else:
    ax.set_ylabel(y_label_abs)
    
ax.set_title(spectrum_title,fontweight=spectrum_title_weight)   #title
ax.get_yaxis().set_ticks([])                                    #remove ticks from y axis
plt.tight_layout()                                              #tight layout
if minor_ticks:
    ax.minorticks_on()                                          #show minor ticks

#x,y axis manipulation for transmittance or absorption style 
if transm_style:
    plt.ylim(max(plt_range_gauss_sum_y)+max(plt_range_gauss_sum_y)*0.1,0) # +10% for labels
else:
    plt.ylim(0,max(plt_range_gauss_sum_y)+max(plt_range_gauss_sum_y)*0.1) # +10% for labels
    
if high_to_low_wn:
    plt.xlim(roundup(max(plt_range_x)),0) # round to next 100
else:
    plt.xlim(0,roundup(max(plt_range_x)))

#tick locations at the beginning and end of the spectrum x-axis, evenly spaced
if linear_locator:
    ax.xaxis.set_major_locator(plt.LinearLocator())

#show grid
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
    plt.savefig(f"{filename}-ir.png", dpi=figure_dpi)

#show the plot
if show_spectrum:
    plt.show()
