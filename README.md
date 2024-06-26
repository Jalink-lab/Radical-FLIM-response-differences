# Radical FLIM response differences

This repository contains an ImageJ Macro to measure and display fluorescence lifetimes of individual cells in time-lapse FLIM experiments (`Analyze_FLIM_timelapse_experiments_FOR_LITE_PAPER.ijm`), and a few R-scripts to further visualize the data, as used in the paper "‘Radical’ differences between two FLIM microscopes affect interpretation of cell signaling dynamics". The R Markdown file `StatsForLITE-1.1.Rmd` contains all the statistics done in the paper, with the output in `Stats for LITE.pdf`.

`Segment_cell_and_display_FLIM_traces.ijm` is a slightly cleaner version of the imageJ macro, with exactly the same functionality.

# ImageJ Macro to segment cells and display FLIM time traces
Input files:
- 2-channel .tif files representing two lifetime-components measured with TCSPC (e.g. files exported from Stellaris/SP8)
  Optionally a third channel with a nuclear marker			
- .fli files from the Lambert Instruments Frequency-Domain FLIM microscope

Brief workflow: 			  
- Segment cells (create labelmap) using Cellpose
- Measure the intensity-weighted lifetime of all labels
- Display the lifetime traces in a graph
- Display the lifetime traces in a kymograph-like image
- Display the average lifetime trace of all labels (solid black line)

Requires the following update sites:
- CLIJ
- CLIJ2
- CSBDeep
- ImageScience
- PTBIOP
- SCF MPI CBG
You also need a working Cellpose Python environment, and the 'Lifetime' and 'Turbo' LUTs.

# Usage
A brief manual will appear here soon.
Some teasers:

Cellpose segmentation:

![image](https://github.com/Jalink-lab/FLIM-timelapse-traces/assets/66722371/60a29899-309c-4488-b134-8f8215d62fbe)

RGB overlay:

![image](https://github.com/Jalink-lab/FLIM-timelapse-traces/assets/66722371/3dde580e-4008-43ab-a76c-02e8028bb9ef)

Time traces plot of all cells:

![image](https://github.com/Jalink-lab/FLIM-timelapse-traces/assets/66722371/9446ef48-5a37-45fb-93c0-e3cedaeba5c9)
