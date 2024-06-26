# orca-ir
A Python 3 script for (hassle-free) plotting of IR spectra from [ORCA](https://orcaforum.kofo.mpg.de) output files with 
peak dectection and annotation.
It combines the stick spectrum with the convoluted spectrum (gaussian line shape). 
The full spectrum or parts of the spectrum (via matplotlib window) can be plotted.

Please also have a look at the interactive [Jupyter Notebook edition](https://github.com/radi0sus/orca_ir/blob/main/orca-ir.ipynb). 
It offers almost the same functionality without the need to maintain a local Python installation.
Please upload an `orca.out` file to the same directory and adjust the code to `ir_data = imp_data('my_orca_calc_with_ir_data.out')`.    
Note: The interactive Matplotlib window does not work very well with *colab*. Replace `%matplotlib widget` with `%matplotlib inline` 
for a better performance. However, this change removes the ability to select a specific region and save the spectrum bitmap.

 <img src='examples\jn1.png' alt='Jupyter Notebook' width=500 align='center'>   

## External modules
 `re` 
 `numpy` 
 `matplotlib`
 `scipy`  
 
## Quick start
 Start the script with:
```console
python3 orca-ir.py filename
```
it will save the plot as PNG bitmap:
`filename-ir.png`

## Command-line options
- `filename` , required: filename
- `-w`  `N` , optional: line width (in cm<sup>-1</sup>) of the gaussian (default is  `N = 15`)
- `-s` , optional: shows the `matplotlib` window
- `-n` , optional: do not save the spectrum
- `-e` , optional: export the line spectrum in a csv-like fashion; filename of the export is input filename + "-mod.dat"

## Script options
There are numerous ways to configure the spectrum in the script:
Check `# plot config section - configure here` in the script. 
Here, you can configure an absorption or transmittance plot for example.
You can even configure the script to plot of the single gaussian functions.

The delimiter for the line spectrum export can be changed by changing the value of `export_delim =`.

## Code options
Colors, line thickness, line styles, level of peak detection and 
more can be changed in the code directly.

## Remarks
The spectrum always starts at zero and ends at the maximum wave number. 
If you need only a part of the spectrum, you can start the script with:
```console
python3 orca-ir.py filename -s
```
and use the matplotlib window to zoom to an area of interest and save it.
The PNG file will be replaced everytime you start the script with the same output file. 
If you want to keep the file, you have to rename it. 

## Examples:
![show](/examples/show-use4.gif)
![Example 1](/examples/example1.png)
![Example 2](/examples/example2.png)
![Example 3](/examples/example3.png)
