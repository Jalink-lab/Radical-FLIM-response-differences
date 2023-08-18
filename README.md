# FLIM-timelapse-traces

This repository contains an ImageJ Macro to measure and display fluorescence lifetimes of individual cells in time-lapse experiments, and a few R-scripts to further visualize the data.

# ImageJ Macro
Input files:
- 2-channel .tif files representing two lifetime-components measured with TCSPC (e.g. files exported from Stellaris/SP8)
  Optionally a third channel with a nuclear marker			
- .fli files from the Lambert Instruments Frequency-Domain FLIM microscope

Brief workflow: 			  
  ► Segment cells (create labelmap) using Cellpose
  ► Measure the intensity-weighted lifetime of all labels
  ► Display the lifetime traces in a graph
  ► Display the lifetime traces in a kymograph-like image
  ► Display the average lifetime trace of all labels (solid black line)

Requires the following update sites:
- CLIJ
- CLIJ2
- CSBDeep
- ImageScience
- PTBIOP
- SCF MPI CBG
You also need a working Cellpose Python environment, and the 'Lifetime' and 'Turbo' LUTs.

# Usage
A brief manual will appear here soon
