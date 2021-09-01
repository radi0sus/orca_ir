# orca-ir
Python 3 script for (hassle-free) plotting IR spectra from ORCA output files. 
It combines the stick spectrum with the convoluted spectrum (gaussian line shape). 
The full spectrum or parts of the spectrum (via matplotlib window) can be plotted.

## External Modules
 `re`, `numpy`, `matplotlib`, `scipy`  

## Quick start
 Start the script with:
`python3 orca-ir.py filename`
it will save a plot as PNG bitmap:
`filename-ir.png`


## Command-line options
- `filename` required: filename
- `-w` optional: line width of the gaussian (default is 15)
- `-s` optional: shows the `matplotlib` window
- `-n` optional: do not save the spectrum

## Script options
There are numerous ways to configure the spectrum in the script:
Check `# plot config section - configure here` in the script. 
Here, you can configure an absorption or transmittance plot. 
You can even plot the single gaussian functions.

## Code options
Colors, line thickness, line styles, level peak detection and 
more can be changed in the code directly.

## Special options and limitations
The spectrum always starts at zero and ends at the maximum wave number. 
If you need only a part of the spectrum, you can start the script with:
`python3 orca-ir.py filename -s`
and use the matplotlib window to zoom to an area of interest and save it.

The PNG file will be replaced everytime you start the script with the same output file. 
If you want to keep a file, you have to rename it. 

## Examples:
![](/examples/example1.png=250x250)
