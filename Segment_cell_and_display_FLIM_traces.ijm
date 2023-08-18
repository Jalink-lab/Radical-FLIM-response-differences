/*	Macro to measure and display fluorescence lifetimes of individual cells in time-lapse experiments.
 * 	Input:
 * 	► 2-channel .tif files representing two lifetime-components measured with TCSPC (e.g. files exported from Stellaris/SP8)
 * 	  Optionally a third channel with a nuclear marker			
 * 	► .fli files from the Lambert Instruments Frequency-Domain FLIM microscope
 *
 *	Brief workflow: 			  
 *	► Segment cells (create labelmap) using Cellpose
 *	► Measure the intensity-weighted lifetime of all labels
 *	► Display the lifetime traces in a graph
 *	► Display the lifetime traces in a kymograph-like image
 *	► Display the average lifetime trace of all labels (solid black line)
 *	
 *  Requires the following update sites:
 *  - CLIJ
 *  - CLIJ2
 *  - CSBDeep
 *  - ImageScience
 *  - PTBIOP
 *  - SCF MPI CBG
 *  You also need a working Cellpose Python environment, and the 'Turbo' LUT (https://github.com/cleterrier/ChrisLUTs/blob/master/Turbo.lut)
 * 
 *  Author: Bram van den Broek, The Netherlands Cancer Institute, b.vd.broek@nki.nl
 *	
 *	
 */

version = 1.0;

#@ String (label = "Input file microscopy type", choices={"Confocal TCSPC / TauSeparation", "Fast FLIM", "TauContrast", "Frequency Domain FLIM", "Ratio Imaging", "Intensity only"}, style="listBox") microscope
#@ File[] (label = "Input files (exported .tif, .lif or .fli files)", style="File") file_list

#@ File (label = "Output folder", style = "directory") output
#@ Boolean(label="Segment all images before analyzing (faster, but dimensions must be the same)", value=true) segmentFirstAnalyzeLater

#@ String (value="<html><br>TCSPC / TauSeparation / TauContrast settings<hr></html>", visibility="MESSAGE") confocal_message
#@ Integer(label="TCSPC fitted / TauContrast / TauSeparation first channel", value=1, description="TauContrast saves with an additional (hidden) intensity channel. Choose the channel in LAS X that shows the TauContrast lifetime.") intensityChannel
#@ Double(label="Lifetime component 1 (ns)", value=0.6, description="For fitted TCSPC or TauSeparation data", style="format:0.00") tau1
#@ Double(label="Lifetime component 2 (ns)", value=3.4, description="For fitted TCSPC or TauSeparation data", style="format:0.00") tau2
#@ Boolean(label="Correct bidirectional phase mismatch",value=false) bidirectional
#@ Boolean(label="Correct xy drift (using intensity)",value=false) correctDrift_boolean
#@ String (label = "Registration against [Drift correction]", choices={"First frame","Last frame","Previous frame"}, style="listBox") registration_against
//#@ Boolean (label = "Edge-detect before registration? [Drift correction]", value=false) edge_detect

#@ Boolean(label="Remove last frame (sometimes (partially) empty)", value=false) removeLastFrame_boolean

#@ String (value="<html><br>FD-FLIM settings<hr></html>", visibility="MESSAGE") FDFILM_message
#@ File(label = "Reference file", style = "file", default="-") reference
#@ Integer(label="Number of phases", value = 12) phases
#@ Integer(label="Frequency (MHz)", value = 40) freq
#@ Double(label="Lifetime of the reference", value = 3.93) tau_ref
#@ Double(label="Default frame interval (s) (if not found)", value = 5, style="format:#.0") default_frameInterval

#@ String (value="<html><br>Cell segmentation settings<hr></html>", visibility="MESSAGE") segmentation_message
#@ String (label = "Cellpose model", choices={"cyto","cyto2","nuclei"}, style="listBox", value="cyto2") CellposeModel
//#@ Boolean(label="Enhance contrast before cell segmentation (Gamma correction)", value=false) equalize_contrast_for_cellpose
#@ Integer(label="Start frame for segmentation (-1 for the full timelapse)", value=-1) CellposeStartFrame
#@ Integer(label="End frame for segmentation (-1 for the full timelapse)", value=-1) CellposeEndFrame

#@ Integer(label="Cell diameter (0 for automatic)", value=20) CellposeDiameter
#@ Double(label="Cell flow error threshold (Cellpose default=0.4, higher -> more cells)", style="scroll bar", value=1, min=0, max=3, stepSize=0.1) CellposeFlowThreshold
#@ Double(label="Cell probability threshold (default=0, higher -> smaller cell area)", style="scroll bar", value=0, min=-6, max=6, stepSize=0.25) CellposeProbability
#@ Integer(label="Minimum cell size (pixels)", value=100) minCellSize
#@ Integer(label="Measure intensity in additional channel (-1 if N/A)", value=-1) additionalChannel

#@ String (value="<html><br>Output display settings<hr></html>", visibility="MESSAGE") display_message
#@ String (label = "Lookup table", choices={"Lifetime", "Turbo", "Fire", "mpl-viridis", "mpl-plasma", "mpl-viridis", "phase", "glow", "Grays"}, style="listBox", value="Turbo") lut
#@ Double(label="Min. displayed lifetime(ns)", value=2.0, style="format:0.0") minLifetime
#@ Double(label="Max. displayed lifetime(ns)", value=3.4, style="format:0.0") maxLifetime
#@ Double(label="Smooth traces in graph with radius", value=2.0, min=0.0, style="format:0.0") smoothRadiusTraces
#@ Double(label="Smooth lifetime in overlay movie with radius (x and y)", value=1.0, min=0.0, style="format:0.0") smoothRadiusOverlayXY
#@ Double(label="Smooth lifetime in overlay movie with radius (time)", value=1.0, min=0.0, style="format:0.0") smoothRadiusOverlayTime

#@ String (label = "Create RGB lifetime overlay on", choices={"average intensity image", "intensity movie"}, style="radioButtonHorizontal") overlayMovie
#@ Double(label="Brightness of RGB overlay [0-5]", value=3, min=0, max=5, style="format:0.0") RGB_brightness
#@ String (label = "Calibration bar position", choices={"Upper Right","Upper Left","Lower Right","Lower Left"}, style="listBox", value="cyto2") calibrationBarPosition
#@ Boolean(label="Display grid lines in plot", value=true) displayGrid
#@ Integer(label="Output plot axis font size", value=18) axisFontSize

#@ Boolean(label="Debug mode", value=false) debugMode

//Other adjustable parameters
equalize_contrast_for_cellpose = false;	//Gamma correction before cell segmentation
edge_detect = false;					//Edge-detect before drift correction (Can yield better results)
histogramBins = 50;						//Number of bins in the histograms (for non-timelapse data)
saturatedPixels = 0.35;					//Contrast settings before applying the LUT for Cellpose - Some saturation gives better results if there are also dim cells
sigma_CC = 4;							//Blur Cross Correlation with sigma (pixels)
show_CC = false;						//Show Cross Correlation image (drift correction)
subpixel = false						//Requires the Fiji update site 'ImageScience'. If enabled, the shift positions are calculated using the center of mass of (a crop of) the cross correlation images, instead of the maximum. This allows for subpixel accuracy, but it is not always better.
displayConfidenceInterval = false;		//Display red lines in the traces plot depicting the confidence interval
confidence_interval_sigma = 3;			//Number of standard deviations in the confidence interval
cellposeVersion = 2.0;					//Version of Cellpose that you use. (The plugin is called slightly differently for version 0.6).

lifetimeChannel = intensityChannel + 1;	//The lifetime channel comes after the intensity channel
RGB_brightness = pow(10, (0.5*(RGB_brightness-2)));
var pixelWidth;

output = output + File.separator;
if(!File.exists(output)) {
	print("Output directory "+output+" does not exist. Creating it.");
	File.makeDirectory(output);
}

saveSettings();

run("Colors...", "foreground=white background=black selection=gray");
run("Conversions...", " ");
print("\\Clear");

run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();

close_windows("Lifetime_Data");
close_windows("Cell_statistics");

//Skip files that are not .fli, .tif or .lif
Array.sort(file_list);
for(i=0; i<file_list.length; i++) {
	if( (!endsWith(file_list[i], ".fli")) && (!endsWith(file_list[i], ".tif")) && (!endsWith(file_list[i], ".lif"))) {
		print(file_list[i] + " is not a valid file - skipping it.");
		file_list = Array.deleteIndex(file_list, i);
		i--;
	}
}
print(file_list.length + " files to analyze.");

//Loop over all files. Open intensity and lifetime from .tif file, and open nuclei from the .lif file, if present.
alreadySegmented = false;
if(file_list.length == 1) segmentFirstAnalyzeLater = false;
for (f = 0; f < file_list.length; f++) {
	if (microscope == "Confocal TCSPC / TauSeparation" || microscope == "TauContrast") {
		run("Bio-Formats Macro Extensions");
		Ext.setId(file_list[f]);
		Ext.getSeriesCount(nr_series);
	}
	else nr_series = 1;
	
	//Loop over all series
	for(s = 0; s < nr_series; s++) {
		run("Close All");
		if(isOpen("Lifetime_Data")) {
			selectWindow("Lifetime_Data");
			run("Close");
		}
		if(isOpen("Intensity_table_ch"+additionalChannel)) {
			selectWindow("Intensity_table_ch"+additionalChannel);
			run("Close");
		}
		if(isOpen("ROI Manager")) {
			selectWindow("ROI Manager");
			run("Close");
		}
		setBatchMode(false);
		run("ROI Manager...");
		setBatchMode(true);
		roiManager("reset");
		roiManager("Set Color", "grays");
		roiManager("Set Line Width", 0);
	
		setBatchMode(true);

		run("Bio-Formats Macro Extensions");			//Run extensions again, because only one macro extension can be loaded
		Ext.setId(file_list[f]);
		Ext.setSeries(s);								//Series start at 0 here
		Ext.getSeriesName(seriesName)
		seriesName = replace(seriesName,"\\/","-");		//replace slashes by dashes in the seriesName
		if(seriesName != File.getName(file_list[f])) {
			saveName = File.getNameWithoutExtension(file_list[f]) + " - " + seriesName;
		}
		else saveName = File.getNameWithoutExtension(file_list[f]);

		run("CLIJ2 Macro Extensions", "cl_device=");	//Run extensions again, because only one macro extension can be loaded
		Ext.CLIJ2_clear();	
	
		//open the image
		if(segmentFirstAnalyzeLater == false || alreadySegmented == true) {
//			!!WARNING!! Fitted TCSPC .tif images exported from LAS X are saved with unit 'pixels'. Bio-Formats then ignores the pixel calibration! Using the standard ImageJ opener below is a workaround, but it doesn't work for multiseries files.
//			if (microscope == "Confocal TCSPC / TauSeparation") run("Bio-Formats Importer", "open=["+file_list[f]+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+1);
			if (microscope == "Confocal TCSPC / TauSeparation") open(file_list[f]);
			else if (microscope == "TauContrast") run("Bio-Formats Importer", "open=["+file_list[f]+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+1);
			else if (microscope == "Fast FLIM") open(file_list[f]);
			else if (microscope == "Frequency Domain FLIM") openfli(file_list[f], saveName);
			else if (microscope == "Ratio Imaging") open(file_list[f]);
			else if (microscope == "Intensity only") open(file_list[f]);
		}
		//Run this only once if segmentFirstAnalyzeLater == true
		if(segmentFirstAnalyzeLater == true && alreadySegmented == false && microscope != "Frequency Domain FLIM") {
			for(file=0; file<file_list.length; file++) {
				open(file_list[file]);
				input_image = getTitle();
				getDimensions(width, height, channels, slices, frames);
				print("Opening "+File.getName(file_list[file]));
				if(removeLastFrame_boolean == true && frames>1) {
					Stack.setFrame(frames);
					run("Delete Slice", "delete=frame");
					getDimensions(width, height, channels, slices, frames);
				}
				if(correctDrift_boolean == true) correct_drift(input_image);
			}
			if(file_list.length>1) run("Concatenate...", "all_open title=hyperstack open");
			getDimensions(width, height, channels, slices, frames);
			run("Stack to Hyperstack...", "order=xyctz channels="+channels+" slices="+slices*file_list.length+" frames="+frames/file_list.length+" display=Grayscale");
			if(frames>1 && slices>1) run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
			Stack.setDisplayMode("grayscale");
			rename("Stack_all");
			input_image = getTitle();
		}
		else if(segmentFirstAnalyzeLater == false || alreadySegmented == true) {
			rename(saveName);
			input_image = getTitle();
			getDimensions(width, height, channels, slices, frames);
			print("\nAnalyzing "+saveName+"\n");

			if(removeLastFrame_boolean == true && frames>1) {
				Stack.setFrame(frames);
				run("Delete Slice", "delete=frame");
				getDimensions(width, height, channels, slices, frames);
			}

		if(correctDrift_boolean == true) correct_drift(input_image);
		}
		//Get pixel size from .tif file
		getPixelSize(unit, pixelWidth, pixelHeight);

		selectWindow(input_image);
		getDimensions(width, height, channels, slices, frames);
		frameInterval = Stack.getFrameInterval();
		if(frameInterval != 0) print("Frame Interval detected: "+frameInterval+" s");
		else {
			print("Warning: Frame Interval not found! Using manual value of "+default_frameInterval+" s.");
			frameInterval = default_frameInterval; 
		}
		//TO DO: get frame interval from metadata for FDFLIM:
		//  FLIMIMAGE: TIMESTAMPS - t0 = 30905379 4265506959
		//  FLIMIMAGE: TIMESTAMPS - t1 = 30905380 526102
		
		for(c=channels; c>=1; c--) {
			Stack.setChannel(c);
			run("Enhance Contrast", "saturated=0.35");
		}
		
		//Retrieve the intensity and lifetime stacks using the correct modality
		if(microscope == "Confocal TCSPC / TauSeparation") lifetime_and_intensity_stacks = calculate_lifetime_and_intensity_TCSPC(input_image);
		else if(microscope == "TauContrast") lifetime_and_intensity_stacks = calculate_lifetime_and_intensity_TauContrast(input_image);
		else if(microscope == "FAST FLIM") lifetime_and_intensity_stacks = calculate_lifetime_and_intensity_FASTFLIM(input_image);
		else if(microscope == "Frequency Domain FLIM") lifetime_and_intensity_stacks = calculate_lifetime_and_intensity_FDFLIM(input_image);
		else if(microscope == "Ratio Imaging") lifetime_and_intensity_stacks = calculate_ratio_and_intensity(input_image);
		else if(microscope == "Intensity only") lifetime_and_intensity_stacks = calculate_intensity(input_image);
		intensity_stack = lifetime_and_intensity_stacks[0];
		selectWindow(intensity_stack);
		setBatchMode("show");
		lifetime_stack = lifetime_and_intensity_stacks[1];
		selectWindow(lifetime_stack);
		setBatchMode("show");

		//Segment the cells
		if(alreadySegmented == false) {
			labelmaps = segment_cells_no_nuclei(intensity_stack);
			labelmap_cells = labelmaps[0];
			nr_cells = getValue("Max");
		}

		//Save all labelmaps if segmentFirstAnalyzeLater == true and break out of the loop
		if(segmentFirstAnalyzeLater == true && alreadySegmented == false) {
			File.makeDirectory(output+"labelmaps");
			selectWindow(labelmap_cells);
			rename("labelmap_cells_");
			run("Image Sequence... ", "dir=["+output+"labelmaps] format=TIFF digits=4");
			File.makeDirectory(output+"intensities");
			selectWindow("intensity_image_for_Cellpose");
			rename("intensity_");
			run("Image Sequence... ", "dir=["+output+"intensities] format=TIFF digits=4");
//				File.makeDirectory(output+"lifetimes");
//				selectWindow(lifetime_stack);
//				run("Image Sequence... ", "dir="+output+"lifetimes format=TIFF digits=4");				
			alreadySegmented = true;
			f--;
			continue;
		}
		//Open a single labelmap and intensity image from disk
		if(segmentFirstAnalyzeLater == true && alreadySegmented == true) {
			print("Analyzing "+File.getName(file_list[f]));
			open(output+"labelmaps"+File.separator+"labelmap_cells_"+IJ.pad(f,4)+".tif");
			labelmap_cells = getTitle();
			nr_cells = getValue("Max");
			setBatchMode("show");
			open(output+"intensities"+File.separator+"intensity_"+IJ.pad(f,4)+".tif");
			rename("intensity_image_for_Cellpose");
			setBatchMode("show");
			run("Duplicate...", "title=intensity");	//Necessary later to create the RGB overlay
			setBatchMode("show");
//				open(output+"lifetimes"+File.separator+"intensity_stack"+IJ.pad(f,4)+".tif");
//				lifetime_stack = getTitle();
//				setBatchMode("show");
		}
		
		if(nr_cells == 0) {
			print("No cells found in this image!");
			s++;
			break;
		}
		//Overlay labelmap with intensity image
		selectWindow("intensity_image_for_Cellpose");
		run("Add Image...", "image="+labelmap_cells+" x=0 y=0 opacity=33 zero");

		//Measure the lifetime traces
		lifetimeInfo = measure_lifetime_traces(intensity_stack, lifetime_stack, labelmap_cells, saveName);
		lifetimeTable = lifetimeInfo[0];
		kymograph = lifetimeInfo[1];

		if(additionalChannel > 0) {
			additionalChannelTable = measure_intensity_in_additional_channel(input_image, additionalChannel);
		}
		
		if(microscope != "Intensity only") {
			if(overlayMovie == "intensity movie") RGB_overlay = overlay_intensity(intensity_stack, lifetime_stack, saveName, smoothRadiusOverlayXY, smoothRadiusOverlayTime);
			else if(overlayMovie == "average intensity image") RGB_overlay = overlay_intensity("intensity", lifetime_stack, saveName, smoothRadiusOverlayXY, smoothRadiusOverlayTime);
		}

		//Save images and data tables
		run("Set Measurements...", "mean redirect=None decimal=9");	//Make sure that enough decimals are saved!
		selectWindow(lifetime_stack);
		saveAs("tiff", output + saveName + " (lifetime)");
		selectWindow(labelmap_cells);
		saveAs("tiff", output + saveName + " (labelmap_cells)");
		selectWindow("intensity_image_for_Cellpose");
		saveAs("tiff", output + saveName + " (intensity & labelmap)");

		selectWindow(kymograph);
		saveAs("tiff", output + saveName + " (kymograph)");

		if(microscope != "Intensity only") {
			selectWindow(RGB_overlay);
			if(frames>1) Stack.setFrame(1);
			updateDisplay();
			saveAs("tiff", output + saveName + " (lifetime & intensity RGB overlay)");
		}

		roiManager("deselect");
		roiManager("Remove Frame Info");
		roiManager("save", output + saveName + " (ROIs).zip");

		Ext.CLIJ2_clear();
	}
}

restoreSettings();
//END


function getLabelColor(label, nrOfLabels) {
	color1 = IJ.pad(toHex(reds[label/nrOfLabels*255]),2);
	color2 = IJ.pad(toHex(greens[label/nrOfLabels*255]),2);
	color3 = IJ.pad(toHex(blues[label/nrOfLabels*255]),2);
	labelColor = "#"+color1+color2+color3;
	return labelColor;
}


function labels_to_ROI_Manager(labelmap) {
//	selectWindow(labelmap);
	Ext.CLIJ2_push(labelmap);
//	run("Clear Results");
//	Ext.CLIJ2_statisticsOfLabelledPixels(labelmap, labelmap);	//This is in fact redundant, because the statistics are already present.
	boundingBox_X = Table.getColumn("BOUNDING_BOX_X", "Cell_statistics");
	boundingBox_Y = Table.getColumn("BOUNDING_BOX_Y", "Cell_statistics");
	boundingBox_width = Table.getColumn("BOUNDING_BOX_WIDTH", "Cell_statistics");
	boundingBox_height = Table.getColumn("BOUNDING_BOX_HEIGHT", "Cell_statistics");
	Array.getStatistics(boundingBox_width, min, boundingBoxMax_X, mean, stdDev);
	Array.getStatistics(boundingBox_height, min, boundingBoxMax_Y, mean, stdDev);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap, nr_cells);

	if(isOpen("ROI Manager")) { selectWindow("ROI Manager"); }//run("Close"); }	//This step goes faster when the ROI manager is not visible.
	for (i = 0; i < nr_cells; i++) {
		showStatus("Converting "+nr_cells+" labels to ROIs...");
		Ext.CLIJ2_crop2D(labelmap, label_cropped, boundingBox_X[i], boundingBox_Y[i], boundingBoxMax_X, boundingBoxMax_Y);
		Ext.CLIJ2_labelToMask(label_cropped, mask_label, i+1);
		Ext.CLIJ2_pullToROIManager(mask_label);
		roiManager("Select",i);
		Roi.move(boundingBox_X[i], boundingBox_Y[i]);
		roiManager("update");
//		else {	//Else the label is square and doesn't contain a zero, causing a crash because it is not added to the ROI manager.
//			//TO DO: ReplaceIntensities on the labelmap in GPU memory, closeIndexGapsInLabelMap and, ultimately, pull and overwrite the final labelmap
//		}
	}
	Ext.CLIJ2_release(mask_label);
	Ext.CLIJ2_release(label_cropped);
	roiManager("Deselect");
	roiManager("Remove Frame Info");
}


function measure_intensity_in_additional_channel(image, additionalChannel) {
	selectWindow(image);
	run("Duplicate...", "title="+saveName+"_ch"+additionalChannel+" duplicate channels="+additionalChannel);
	if(isOpen("Results")) { selectWindow("Results"); run("Close"); }	//If the Results window is open measurements take longer!
	run("Set Measurements...", "mean redirect=None decimal=3");
	roiManager("deselect");
	showStatus("Measuring intensities...");
	roiManager("Multi Measure");
	selectWindow("Results");
	Table.rename("Results", "Intensity_table_ch"+additionalChannel);
	for (i = 0; i < nr_cells; i++) Table.renameColumn("Mean(cell_"+i+1+")", "cell_"+IJ.pad(i+1,4));
	Table.update;
	Table.save(output + saveName + "_intensity_ch_"+additionalChannel+".tsv");
	close(saveName+"_ch"+additionalChannel);
	return "Intensity_table_ch"+additionalChannel;
}


function measure_lifetime_traces(intensity_stack, lifetime_stack, labelmap_cells, saveName) {

	//Create masked intensity stack in case not all pixels in the lifetime image have a value
	Ext.CLIJ2_push(intensity_stack);
	Ext.CLIJ2_push(lifetime_stack);
	intensity_stack_masked = "intensity_masked";
	Ext.CLIJ2_mask(intensity_stack, lifetime_stack, intensity_stack_masked);
	Ext.CLIJ2_release(intensity_stack);
	Ext.CLIJ2_pull(intensity_stack_masked);
	setBatchMode("show");

	if(microscope == "Frequency Domain FLIM") frames = frames/phases;
	
	run("Enhance Contrast", "saturated=0.35");
	for(i=1;i<=frames;i++) {
		if(frames>1) Stack.setFrame(i);
		changeValues(0,0,NaN);	//Set zeroes to NaN in intensity stack
	}
	
	//For FDFLIM: set values <100 to NaN in intensity stack (non-illuminated corners) - disabled, because it can give artifacts!
	/*
	if(microscope == "Frequency Domain FLIM") {
		selectWindow(intensity_stack_masked);
		for(i=1;i<=frames;i++) {
			setSlice(i);
			changeValues(0,100,NaN);
		}
	}
	*/
	run("Set Measurements...", "area mean redirect=None decimal=3");

	//Convert labelmap to ROIs
	//Ideally we would measure using CLIJ2 and the labelmap, but that doesn't allow NaNs,
	//and also is not great for time-lapses. Besides, I have this already working.
	selectWindow(labelmap_cells);
	getLut(reds, greens, blues);
	
	labels_to_ROI_Manager(labelmap_cells);
	nr_cells = roiManager("count");

	//measure and save intensities
	selectWindow(intensity_stack_masked);
	setBatchMode("hide");
	if(isOpen("Results")) { selectWindow("Results"); run("Close"); }	//If the Results window is open measurements take longer, so we close it!
	run("Set Measurements...", "mean redirect=None decimal=3");
	selectWindow(intensity_stack_masked);
	roiManager("deselect");
	showStatus("Measuring intensities...");
	roiManager("Multi Measure");
	selectWindow("Results");
	Table.rename("Results", "Intensity_table");
	Table.save(output + saveName + "_intensity.tsv");
//	renameTableHeaders("Intensity_table");

	//multiply lifetime with intensity (for normalization)
	imageCalculator("Multiply create 32-bit stack", lifetime_stack, intensity_stack_masked);
	rename("Lifetime_times_intensity");
	
	Ext.CLIJ2_release(intensity_stack_masked);
	close(intensity_stack_masked);

	//measure lifetimes times intensity
	if(isOpen("Results")) { selectWindow("Results"); run("Close"); }
	run("Set Measurements...", "mean redirect=None decimal=3");
	selectWindow("Lifetime_times_intensity");
	roiManager("deselect");
	showStatus("Measuring lifetime*intensity...");
	roiManager("Multi Measure");
	selectWindow("Results");
	Table.rename("Results", "Lifetime*Intensity_table");
//	renameTableHeaders("Lifetime*Intensity_table");

	//Normalize lifetime with cell intensity (every pixel is weighted with intensity) and create tables and plot
	showStatus("Computing lifetimes...");
	selectWindow("Lifetime*Intensity_table");
	headings = Table.headings("Lifetime*Intensity_table");
	headers = split(headings, "\t");

	Table.rename("Lifetime*Intensity_table", "Results");
	lifetimeTimesIntensityImage = "lifetimeTimesIntensityImage";
	Ext.CLIJ2_pushResultsTable(lifetimeTimesIntensityImage);
	Table.rename("Intensity_table", "Results");
	intensityImage = "intensityImage";
	Ext.CLIJ2_pushResultsTable(intensityImage);
	kymograph = "kymograph";
	Ext.CLIJ2_divideImages(lifetimeTimesIntensityImage, intensityImage, kymograph);
	Ext.CLIJ2_release(lifetimeTimesIntensityImage);
	Ext.CLIJ2_release(intensityImage);
	Ext.CLIJ2_pull(kymograph);
	setBatchMode("show");
	run("Clear Results");
	kymograph_smoothed = "kymograph smoothed";
	if(smoothRadiusTraces >0) {
		Ext.CLIJ2_mean2DBox(kymograph, kymograph_smoothed, 0, smoothRadiusTraces);
		Ext.CLIJ2_pull(kymograph_smoothed);
		setBatchMode("show");
		Ext.CLIJ2_pullToResultsTable(kymograph_smoothed);
	}
	else Ext.CLIJ2_pullToResultsTable(kymograph);
	
	//Construct plot of all traces, or a histogram if there is only one time point.
	if(microscope == "Confocal TCSPC / TauSeparation") 	y_axis = "Lifetime (ns)";
	else if(microscope == "TauContrast") 				y_axis = "Lifetime (ns)";
	else if(microscope == "Fast FLIM") 					y_axis = "Lifetime (ns)";
	else if(microscope == "Frequency Domain FLIM")		y_axis = "Lifetime (ns)";
	else if(microscope == "Ratio Imaging")				y_axis = "Ratio (ch1 / ch2)";
	else if(microscope == "Intensity only")				y_axis = "Intensity";

	if(frames > 1) {
		timeArray = Array.getSequence(frames);
		timeArray = multiplyArraywithScalar(timeArray, frameInterval);
		plotName = saveName + " (lifetime traces plot)";
		Plot.create(plotName, "time (s)", y_axis);
		if(microscope != "Intensity only") Plot.setLimits(0, frames*frameInterval, minLifetime, maxLifetime);
		Plot.setFrameSize(900, 600);
		Plot.setAxisLabelSize(axisFontSize);
		Plot.setFontSize(axisFontSize);
		//Plot.setXYLabels("Distance (microns)", "Gray Value");
	}
	else {
		plotName = saveName + " (lifetime histogram plot)";
		Plot.create(plotName, y_axis, "count");
		Plot.setFrameSize(640, 480);
		//Plot.setLimits(0, frames*frameInterval, minLifetime, maxLifetime);
	}
	for(i=0;i<nr_cells;i++) {
		showStatus("Computing lifetime... "+i+"/"+nr_cells);
		showProgress(i/nr_cells);
		lifetimeData = Table.getColumn("X"+i, "Results");
		Table.renameColumn("X"+i, "cell_"+IJ.pad(i+1,4));
		if(frames > 1) {
			//Generate lines with 'glasbey on dark' colors
			color = getLabelColor(i, nr_cells);
			Plot.setColor(color);
			Plot.add("line", timeArray, lifetimeData);
		}
		roiManager("select",i);
		roiManager("rename", "cell_"+i+1);
//		roiManager("Set Color", "#"+color1+color2+color3);	//Works, but looks very confusing
	}
	run("Select None");

	//Calculate average trace
	setBatchMode(false);
	selectWindow("kymograph");
	run(lut);
	run("Rotate 90 Degrees Left");
	run("Flip Vertically");
	if(smoothRadiusTraces > 0) {
		selectWindow("kymograph smoothed");
		run(lut);
		run("Rotate 90 Degrees Left");
		run("Flip Vertically");
	}
	setBatchMode("hide");
	run("Select All");
	run("Plots...", "minimum=0");	//Make sure the profile is horizontal
	average_lifetime_trace = getProfile();
	run("Select None");
	
	//Calculate lifetime standard deviation for all cells at all timepoints
	Ext.CLIJ2_transposeXZ(kymograph, kymograph_XZ);
	Ext.CLIJ2_standardDeviationZProjection(kymograph_XZ, kymograph_stddev_XZ);
	Ext.CLIJ2_transposeXY(kymograph_stddev_XZ, kymograph_stddev);
	Ext.CLIJ2_pull(kymograph_stddev);
	Ext.CLIJ2_release(kymograph_XZ);
	Ext.CLIJ2_release(kymograph_stddev_XZ);
	Ext.CLIJ2_release(kymograph_stddev);
	selectWindow(kymograph_stddev);
	run("Select All");
	stddev_avg_lifetime_trace = getProfile();
	confidence_interval = multiplyArraywithScalar(stddev_avg_lifetime_trace, confidence_interval_sigma);

	//Add average trace to the plot and save
	if(frames > 1) {
		Plot.setColor("black");
		Plot.setLineWidth(3);
		Plot.add("line", timeArray, average_lifetime_trace);
		Plot.setLineWidth(1);
		average_trace_confidence_max = addArrays(average_lifetime_trace, confidence_interval);
		average_trace_confidence_min = subtractArrays(average_lifetime_trace, confidence_interval);
		if(displayConfidenceInterval == true) {
			Plot.setColor("red");
			Plot.setLineWidth(2);
			Plot.add("line", timeArray, average_trace_confidence_max);
			Plot.add("line", timeArray, average_trace_confidence_min);
		}
		if(microscope != "Intensity only") Plot.setLimits(0, frames*frameInterval, minLifetime, maxLifetime);
		else Plot.setLimitsToFit();
		if(displayGrid == true) Plot.setFormatFlags("11000000111111");
		else Plot.setFormatFlags("11000000001111");
	}
	else {	//Create histogram instead of time trace plot
		selectWindow(kymograph);
		getStatistics(area, mean, min, max, std);
		getHistogram(values, counts, histogramBins, minLifetime, maxLifetime);
		Plot.add("bar", values, counts);
		Plot.setStyle(0, "#0000a0, #00a0ff, 1.0, Separated Bars");
		Plot.setJustification("right");
		Plot.setFontSize(axisFontSize);
		Plot.addText("Mean ± sd: "+d2s(mean,2)+" ± "+d2s(std,2), 0.97, 0.08);
	}
	updateDisplay();
	Plot.show();
	setBatchMode("show");
	saveAs("Tiff", output + plotName);	

	selectWindow(kymograph);
	run("Rotate 90 Degrees Right");	//Rotate back
//	saveAs("tiff", output + saveName + " (kymograph)");
//	rename(kymograph);

	lifetimeTable = "Lifetime_Data";
	Table.rename("Results", lifetimeTable);
	if(frames > 1) Table.setColumn("time (s)", timeArray, lifetimeTable);
	Table.save(output + saveName + "_lifetime.tsv");

	return newArray(lifetimeTable, kymograph);
}


function calculate_lifetime_and_intensity_TCSPC(image) {
	selectWindow(image);
	//Get rid of the extension in the name (if any)
	//name = substring(image, 0, lastIndexOf(image, "."));
	run("Duplicate...", "title=IntensityStack duplicate");
	run("32-bit");
	run("Split Channels");
	imageCalculator("Add create 32-bit stack", "C"+intensityChannel+"-IntensityStack","C"+lifetimeChannel+"-IntensityStack");
	rename("Intensity");

	//Multiply amplitudes with lifetime components to get total intensities of the components
	selectWindow("C"+intensityChannel+"-IntensityStack");
	run("Multiply...", "value="+tau1+" stack");
	selectWindow("C"+lifetimeChannel+"-IntensityStack");
	run("Multiply...", "value="+tau2+" stack");
	imageCalculator("Add create 32-bit stack", "C"+intensityChannel+"-IntensityStack","C"+lifetimeChannel+"-IntensityStack");
	rename("Total_nr_photons");
	imageCalculator("Divide create 32-bit stack", "Total_nr_photons","Intensity");
	rename(image+"_lifetime");
	close("C"+intensityChannel+"-IntensityStack");
	close("C"+lifetimeChannel+"-IntensityStack");
	close("Total_nr_photons");

	//Set NaNs to 0 - Otherwise drift correction / movement (?) can cause artifacts.
	getDimensions(width, height, channels, slices, frames);
	if(frames>1) {
		for (f = 1; f <= frames; f++) {
	    	Stack.setFrame(f);
			changeValues(NaN, NaN, 0);
		}
	}
	else changeValues(NaN, NaN, 0);
	
	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_lifetime");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_lifetime");
}


function calculate_lifetime_and_intensity_TauContrast(image) {
	selectWindow(image);
	run("Duplicate...", "title=IntensityStack duplicate channels="+intensityChannel);
	run("32-bit");
	rename("Intensity");

	selectWindow(image);
	run("Duplicate...", "title=["+image+"_lifetime] duplicate channels="+lifetimeChannel);
	run("32-bit");
	//Scale TauContrast numbers to lifetime (-1 to 24 ns range equals 0-255 range)
	run("Multiply...", "value=0.097 stack");
	run("Subtract...", "value=1 stack");
	
//	Set pixels with value -1 (created after drift correction) to NaN - Creates artifacts!
//	setThreshold(-0.999, 24);
//	run("NaN Background", "stack");
//	resetThreshold();

	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_lifetime");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_lifetime");
}


function calculate_lifetime_and_intensity_FASTFLIM(image) {
	selectWindow(image);
	run("Duplicate...", "title=IntensityStack duplicate channels="+intensityChannel);
	run("32-bit");
	rename("Intensity");

	selectWindow(image);
	run("Duplicate...", "title=["+image+"_lifetime] duplicate channels="+lifetimeChannel);
	run("32-bit");
	//Scale TauContrast numbers to lifetime (-1 to 24 ns range equals 0-255 range)
	run("Multiply...", "value=0.001 stack");
	
//	Set pixels with value -1 (created after drift correction) to NaN - Creates artifacts!
//	setThreshold(-0.999, 24);
//	run("NaN Background", "stack");
//	resetThreshold();

	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_lifetime");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_lifetime");
}


function calculate_lifetime_and_intensity_FDFLIM(image) {
	openfli(reference,"reference");
	selectWindow(image);
//	name = substring(image, 0, lastIndexOf(image, "."));
//	getDimensions(w, h, channels, slices, frames);
	frames = frames/phases;
	// loop over all time frames
	for (i = 0; i < (frames); i++) {
		selectWindow(image);
		if (i<((frames)-1)){ //last frame
			run("Make Substack...", "delete slices=1-"+phases);
		}
		rename("temp_img");
		run("fdFLIM", "image1=[temp_img] boolphimod=false image2=reference tau_ref="+tau_ref+" freq="+freq);	
		close("temp_img");
		if (i==0) rename("Lifetimes_final");
		else run("Concatenate...", "  title=Lifetimes_final open image1=Lifetimes_final image2=Lifetimes");
	}
	close("reference");
	close(image);

	//Separate intensity and lifetime
	selectWindow("Lifetimes_final");
	run("Duplicate...", "duplicate slices=3");
	rename("Intensity");
	selectWindow("Lifetimes_final");
	run("Duplicate...", "duplicate slices=1");
	rename(image+"_phase_lifetime");
	close("Lifetimes_final");

	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_phase_lifetime");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");

	return newArray(image+"_intensity", image+"_phase_lifetime");
}

function calculate_ratio_and_intensity(image) {
	run("Duplicate...", "title=IntensityStack duplicate");
	run("32-bit");
	run("Split Channels");
	imageCalculator("Add create 32-bit stack", "C1-IntensityStack","C2-IntensityStack");
	rename("Intensity");
	selectWindow("C1-IntensityStack");
	setThreshold(1,65536);
	run("NaN Background", "stack");
	selectWindow("C2-IntensityStack");
	setThreshold(1,65536);
	run("NaN Background", "stack");

	//Divide intensities to get the ratio
	imageCalculator("Divide create 32-bit stack", "C1-IntensityStack","C2-IntensityStack");
	rename(image+"_ratio");
	close("C1-IntensityStack");
	close("C2-IntensityStack");

	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_ratio");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_ratio");
}

function calculate_intensity(image) {
	selectWindow(image);
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	run("Duplicate...", "title=Intensity duplicate");
	run("32-bit");
	setThreshold(1,65536);
	run("NaN Background", "stack");

	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity_nobackground");
	run("Enhance Contrast", "saturated=0.35");	
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_intensity_nobackground");
}


// Open both sample and background from a .fli file and subtract background
function openfli(input,name) {
	run("Bio-Formats", "open=["+input+"] view=Hyperstack stack_order=XYCZT series_1");
	rename(name);
	run("Bio-Formats", "open=["+input+"] view=Hyperstack stack_order=XYCZT series_2");
	rename("BG");
	imageCalculator("Subtract 32-bit stack", name, "BG");
	close("BG");
	close(name);
	rename(name);
}


function segment_cells_no_nuclei(intensity_stack) {
	//Segment cells using Cellpose
	selectWindow(intensity_stack);
	getDimensions(width, height, channels, slices, frames);
	if(frames > 1 && slices == 1) {			//timelapse, single image
		if(CellposeStartFrame == -1 && CellposeEndFrame == -1) run("Z Project...", "projection=[Sum Slices] all");
		else if(CellposeStartFrame == -1) run("Z Project...", "stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		else if(CellposeEndFrame == -1) run("Z Project...", "start="+CellposeStartFrame+" projection=[Sum Slices] all");
		else run("Z Project...", "start="+CellposeStartFrame+" stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		rename("intensity");
	}
	else if(slices > 1 && frames > 1) {		//timelapse, multiple images
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		if(CellposeStartFrame == -1 && CellposeEndFrame == -1) run("Z Project...", "projection=[Sum Slices] all");
		else if(CellposeStartFrame == -1) run("Z Project...", "stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		else if(CellposeEndFrame == -1) run("Z Project...", "start="+CellposeStartFrame+" projection=[Sum Slices] all");
		else run("Z Project...", "start="+CellposeStartFrame+" stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		rename("intensity");
	}
	else if(slices > 1 && frames == 1) {	//no timelapse, multiple images
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		run("Duplicate...", "title=intensity duplicate");
	}										//single image
	else run("Duplicate...", "title=intensity duplicate");
	
	//Correct for bidirectional phase mismatch - currently only on the projection, for cell segmentation
	if(bidirectional == true) {
		intensity_corrected = correct_bidirectional_phase("intensity");
		close("intensity");
		selectWindow(intensity_corrected);
		rename("intensity");
	}
	intensity_corrected = "intensity";	//rename back to intensity - ugly but working for now, because we use this title later
	selectWindow(intensity_corrected);
	run("Grays");

	Ext.CLIJ2_push(intensity_corrected);	//Push to GPU for intensity measurements
	
	if(equalize_contrast_for_cellpose == true) {
		//run("Square Root");
		run("Gamma...", "value=0.50");
		//run("Enhance Local Contrast (CLAHE)", "blocksize=64 histogram=256 maximum=3 mask=*None*");	//The fancy way, but not really better
		run("Enhance Contrast", "saturated=0.35");
	}
	selectWindow(intensity_corrected);
	setBatchMode("show");

	run("Duplicate...", "title=intensity_image_for_Cellpose duplicate");	//Keep "intensity" for RGB overlay later

	//If 32-bit, first convert to 16-bit, with some contrast enhancement - better segmentation in some cases
	run("Set Measurements...", "mean redirect=None decimal=3");
	percentile_threshold(0.05);	//Threshold on the 5% lowest values
	List.setMeasurements("limit");
	mean = List.getValue("Mean");
	changeValues(NaN, NaN, mean);
	run("Enhance Contrast", "saturated="+saturatedPixels);	//TO DO: another way to reliably set the B&C settings
	run("Conversions...", "scale");
	run("16-bit");

	//Run Cellpose without nuclei channel
	selectWindow("intensity_image_for_Cellpose");
	setBatchMode("show");
	setBatchMode(false);	//Cellpose doesn't return an image in batch mode
	cellpose_success = false;
	selectWindow("intensity_image_for_Cellpose");
	while(cellpose_success == false) {
		if(cellposeVersion == 0.6) run("Cellpose Advanced", "diameter="+CellposeDiameter+" cellproba_threshold="+CellposeProbability+" flow_threshold="+CellposeFlowThreshold+" anisotropy=1.0 diam_threshold=12.0 model="+CellposeModel+" nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false");
		else if(cellposeVersion > 0.6) run("Cellpose Advanced", "diameter="+CellposeDiameter+" cellproba_threshold="+CellposeProbability+" flow_threshold="+CellposeFlowThreshold+" anisotropy=1.0 diam_threshold=12.0 model="+CellposeModel+" nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
		if(getTitle() != "intensity_image_for_Cellpose-cellpose") {
			print("Cellpose failed!\nCheck the Fiji console for more information. Waiting 60 seconds and then trying again...");
			wait(60000);
		}
		else cellpose_success = true;
	}
	rename("labelmap_cells");
	setBatchMode("show");

	labelmap_cells = "labelmap_cells";
	Ext.CLIJ2_push(labelmap_cells);
	close("labelmap_cells");
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_cells, labelmap_cells_filtered, minCellSize, 10e6);
	Ext.CLIJ2_release(labelmap_cells);
	Ext.CLIJ2_closeIndexGapsInLabelMap(labelmap_cells_filtered, labelmap_cells_filtered_gapsclosed);
	
	run("Clear Results");
	Ext.CLIJ2_statisticsOfLabelledPixels(intensity_corrected, labelmap_cells_filtered_gapsclosed);
	Table.rename("Results", "Cell_statistics");
	Table.save(output + saveName + "_Cell_statistics.tsv");
	
	Ext.CLIJ2_release(intensity_corrected);
	Ext.CLIJ2_release(labelmap_cells_filtered);
	Ext.CLIJ2_pull(labelmap_cells_filtered_gapsclosed);
	Ext.CLIJ2_release(labelmap_cells_filtered_gapsclosed);
	rename("labelmap_cells");
	run("glasbey on dark");
	resetMinAndMax();
	setBatchMode("show");
	
	return newArray("labelmap_cells");
}


function correct_drift(image) {
	print("Performing drift correction...");
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
/*	
	//Change z and t dimensions, if applicable
	Stack.getPosition(channel, slice, frame);
	getDimensions(width, height, channels, slices, frames);
	if(frames == 1 && slices == 1) exit("No time- or z-dimension found.");
	else if(frames == 1 && slices >= 1) {
		print("Swapping z and t dimensions.");
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		getDimensions(width, height, channels, slices, frames);
	}
*/
	//Prepare intensity image used for registration
	if (microscope == "Confocal TCSPC / TauSeparation" || microscope == "Ratio Imaging") {
		selectWindow(image);
		run("Duplicate...", "title=C1_Intensity duplicate channels="+parseInt(intensityChannel));
		selectWindow(image);
		run("Duplicate...", "title=C2_Intensity duplicate channels="+parseInt(intensityChannel+1));
		imageCalculator("Add create 32-bit stack", "C1_Intensity","C2_Intensity");
		rename("Intensity_for_registration");
		close("C1_Intensity");
		close("C2_Intensity");
	}
	else run("Duplicate...", "title=Intensity_for_registration duplicate channels="+parseInt(intensityChannel));

	//Apply variance filter when edge detection is enabled
	if(edge_detect) {
		run("Duplicate...", "title=[For_FFT_"+image+"] duplicate");
		run("Variance...", "radius=2 stack");
	}
	
	//Change image size if required (has to be a power of two for the FFT)
	if(!edge_detect) run("Duplicate...", "title=[For_FFT_"+image+"] duplicate");
	dim = pad_image_edges("For_FFT_"+image, "Top-Left");

	reg_image = getTitle();	//image used for registration
	
	Xpos_max = newArray(frames);
	Ypos_max = newArray(frames);
	
	concat_string = "";
	if(registration_against == "First frame") run("Duplicate...", "title=frameA duplicate range=1-1");
	if(registration_against == "Last frame") run("Duplicate...", "title=frameA duplicate range="+frames+"-"+frames);
	for(f=1;f<frames;f++) {
		showProgress(f, frames);
		if(registration_against == "Previous frame") {
			selectWindow(reg_image);
			run("Duplicate...", "title=frameA duplicate range="+f+"-"+f);
		}
		selectWindow(reg_image);
		run("Duplicate...", "title=frameB duplicate range="+f+1+"-"+f+1);
		//Calculate cross-correlation
		run("FD Math...", "image1=frameA operation=Correlate image2=frameB result=CC"+f+" do");
		concat_string+="image"+f+"=CC"+f+" ";
		if(registration_against == "Previous frame") close("FrameA");
		close("FrameB");
	}
	close("FrameA");
	if(registration_against != "Previous frame") close("frameB");
	close("For_FFT_"+image);
	
	//concatenate cross correlation images
	crossCorrelationImage = "Cross Correlation";
	run("Concatenate...", "  title=["+crossCorrelationImage+"]" + concat_string);
	run("Gaussian Blur...", "sigma="+sigma_CC+" stack");
	
	run("Clear Results");
	run("CLIJ2 Macro Extensions", "cl_device=");
	Ext.CLIJ2_clear();
	Ext.CLIJ2_push(crossCorrelationImage);
	
	//determine max X positions for all frames
	Ext.CLIJ2_maximumYProjection(crossCorrelationImage, max_Y_projection);
	Ext.CLIJ2_transposeXZ(max_Y_projection, max_Y_TransPoseXZ);
	Ext.CLIJ2_argMaximumZProjection(max_Y_TransPoseXZ, max_XY_projection, arg_max_X);
	//determine max Y positions for all frames
	Ext.CLIJ2_maximumXProjection(crossCorrelationImage, max_X_projection);
	Ext.CLIJ2_transposeYZ(max_X_projection, max_X_TransPoseYZ);
	Ext.CLIJ2_argMaximumZProjection(max_X_TransPoseYZ, max_XY_projection, arg_max_Y);
	Ext.CLIJ2_transposeXY(arg_max_Y, arg_max_Y_transposed);
	//Put X and Y coordinates in the Results table and then in arrays
	Ext.CLIJ2_combineHorizontally(arg_max_X, arg_max_Y_transposed, max_XY_positions);
	Ext.CLIJ2_pullToResultsTable(max_XY_positions);
	Xpos_max = Table.getColumn("X0", "Results");
	Ypos_max = Table.getColumn("X1", "Results");
	
	//Subpixel (center of mass) determination
	if(subpixel) {
		run("Clear Results");
		if(show_CC) {
			crossCorrelationCropImage = "Normalized Cross Correlation Crop";
			Ext.CLIJ2_create3D(crossCorrelationCropImage, cropSizeCC, cropSizeCC, frames-1, 32);
		}
		for(f=1; f<frames; f++) {
			Ext.CLIJ2_copySlice(crossCorrelationImage, crossCorrelationFrame, f-1);
			Ext.CLIJ2_crop2D(crossCorrelationFrame, crossCorrelationFrameCrop, Xpos_max[f-1]-floor(cropSizeCC/2), Ypos_max[f-1]-floor(cropSizeCC/2), cropSizeCC, cropSizeCC);
			Ext.CLIJx_normalize(crossCorrelationFrameCrop, crossCorrelationFrameCropNormalized);
			Ext.CLIJ2_centerOfMass(crossCorrelationFrameCropNormalized);
			Xpos_max[f-1] = Xpos_max[f-1] + getResult("MassX") - floor(cropSizeCC/2);
			Ypos_max[f-1] = Ypos_max[f-1] + getResult("MassY") - floor(cropSizeCC/2);
			if(show_CC) Ext.CLIJ2_copySlice(crossCorrelationFrameCropNormalized, crossCorrelationCropImage, f-1);
		}
		if(show_CC) Ext.CLIJ2_pull(crossCorrelationCropImage);
		setBatchMode("show");
	}
	Ext.CLIJ2_clear();
	
	//Calculate translations
	translate_x = newArray(frames);
	translate_y = newArray(frames);
	for(f=1; f<frames; f++) {
		selectWindow(crossCorrelationImage);
		if(registration_against == "Previous frame") {
			if(f==1) {
				translate_x[f-1] = Xpos_max[f-1]-dim/2;
				translate_y[f-1] = Ypos_max[f-1]-dim/2;
			}
			if(f>=2) {
				translate_x[f-1] = translate_x[f-2] + Xpos_max[f-1]-dim/2;
				translate_y[f-1] = translate_y[f-2] + Ypos_max[f-1]-dim/2;
			}
		}
		else {
			translate_x[f-1] = Xpos_max[f-1]-dim/2;
			translate_y[f-1] = Ypos_max[f-1]-dim/2;
		}
		if(show_CC) {
			setSlice(f);
			makePoint(Xpos_max[f-1], Ypos_max[f-1], "medium red dot add");
		}
	}
	
	if(show_CC == false) close(crossCorrelationImage);
	
	if(edge_detect || channels>1) close(reg_image);
	close("Intensity_for_registration");
	
	selectWindow(image);
	run("Select None");
	run("Duplicate...", "title=["+image+"_driftcorr] duplicate");
	//print("Translation (pixels):");
	for(c=1; c<=channels; c++) {
		Stack.setChannel(c);
		for(f=1; f<frames; f++) {
			showProgress(f + (c-1)*frames, frames*channels);
			showStatus("Translating channel "+c+", frame "+f+"/"+frames);
			selectWindow(image+"_driftcorr");
			if(registration_against != "Last frame") Stack.setFrame(f+1);
			else Stack.setFrame(f);
			//print("channel "+c+", frame "+f+": "+translate_x[f-1]+", "+translate_y[f-1]);
			if(subpixel) run("Translate...", "x="+translate_x[f-1]+" y="+translate_y[f-1]+" interpolation=Bicubic slice");
			else run("Translate...", "x="+translate_x[f-1]+" y="+translate_y[f-1]+" interpolation=None slice");
		}
	}
	run("Clear Results");
	selectWindow(image);
	rename(image+"_NOT_drift_corrected");
	selectWindow(image+"_driftcorr");
	rename(image);
	setBatchMode("show");
	return image;
}


//expand the canvas to the nearest power of 2 and return the dimension
function pad_image_edges(image, location) {
	selectWindow(image);
	width = getWidth();
	height = getHeight();
	w=1;
	h=1;
	while(width>pow(2,w)) w++;
	while(height>pow(2,h)) h++;
	dim = pow(2,maxOf(w,h));
	//print("dimension of FFT: "+dim);
	run("Canvas Size...", "width="+dim+" height="+dim+" position="+location+" zero");
	return dim;
}


function correct_bidirectional_phase(image) {
//Adjust the phase of bidirectional confocal images

	image = getTitle();
	getDimensions(orgWidth, orgHeight, channels, slices, frames);
	getVoxelSize(pixelWidth, pixelHeight, depth, unit);
	run("Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1.0000");	//Remove pixel calibration (if any)

	//Change image size if required (has to be a power of two for the FFT)
	original = getTitle();
	run("Duplicate...", "title=[For_FFT] duplicate");
	dim = pad_image_edges("For_FFT", "Top-Left");
	getDimensions(width, height, channels, slices, frames);
	image = getTitle();

	run("Reslice [/]...", "output=1 start=Left avoid");
	run("Deinterleave", "how=2 keep");
	
	//Determine pixel shift with cross-correlation
	selectWindow("Reslice of "+image+" #1");
	run("Reslice [/]...", "output=1 start=Top rotate avoid");
	run("Canvas Size...", "width="+width+" height="+width+" position=Center zero");
	rename("Odd");
	getDimensions(_width, _height, _channels, _slices, _frames);
	if(_slices > 1) {
		run("Z Project...", "projection=[Sum Slices]");
		close("Odd");
	}
	rename("Odd");
	
	selectWindow("Reslice of "+image+" #2");
	run("Reslice [/]...", "output=1 start=Top rotate avoid");
	run("Canvas Size...", "width="+width+" height="+width+" position=Center zero");
	rename("Even");
	getDimensions(_width, _height, _channels, _slices, _frames);
	if(_slices > 1) {
		run("Z Project...", "projection=[Sum Slices]");
		close("Even");
	}
	rename("Even");
	
	image1 = "Odd";
	image2 = "Even";
	run("FD Math...", "image1=["+image1+"] operation=Correlate image2=["+image2+"] result=CC do");
	getStatistics(area, mean, min, max, std, histogram);
	run("Subtract...", "value="+min);
	run("Divide...", "value="+max-min);
	resetMinAndMax();
	
	//fit profile to a rectangle
	makeRectangle(width/2-16, width/2-16, 32, 32);
	run("Crop");
	run("Rotate 90 Degrees Right");
	run("Select All");
	profile = getProfile();
	x_points = Array.getSequence(32);
	Fit.doFit("Gaussian", x_points, profile);
	//Fit.plot;
	y0 = Fit.p(2);
	shift = (y0 - 15);
	print("Correcting bidirectional pixel shift: "+shift+", Rsq = "+Fit.rSquared);
	
	//Correct and recombine
	selectWindow("Reslice of "+image+" #1");
	//run("Translate...", "x="+shift+" y=0 interpolation=Bicubic stack");	//Doesn't allow subpixel translation on images 1 pixels high
	run("TransformJ Translate", "x-distance="+shift/2+" y-distance=0.0 z-distance=0.0 interpolation=[Cubic Convolution] background=0.0");
	selectWindow("Reslice of "+image+" #2");
	run("TransformJ Translate", "x-distance="+-shift/2+" y-distance=0.0 z-distance=0.0 interpolation=[Cubic Convolution] background=0.0");
	run("Interleave", "stack_1=[Reslice of "+image+" #1 translated] stack_2=[Reslice of "+image+" #2 translated]");
	run("Reslice [/]...", "output=1.000 start=Top rotate avoid");
	
	close("For_FFT");
	close("Reslice of "+image);
	close("Reslice of "+image+" #1");
	close("Reslice of "+image+" #1 translated");
	close("Reslice of "+image+" #2");
	close("Reslice of "+image+" #2 translated");
	close("Combined Stacks");
	close("Even");
	close("Odd");
	close("CC");

	corrected_image = image + "_phase_corrected";
	rename(corrected_image);
	setVoxelSize(pixelWidth, pixelHeight, depth, unit);
	makeRectangle(0, 0, orgWidth, orgHeight);
	run("Crop");
	run("Select None");
	
	getDimensions(width, height, channels, slices, frames);
	if(slices>1 && frames==1) run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	return corrected_image;
}


//Overlay intensity image with the smoothed lifetime image, and add a calibration bar 
function overlay_intensity(intensity_image, lifetime_stack, saveName, smoothRadiusOverlayXY, smoothRadiusOverlayTime) {
	selectWindow(lifetime_stack);
	getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(smoothRadiusOverlayXY > 0 || smoothRadiusOverlayTime > 0) {
		lifetime_stack_for_smoothing = lifetime_stack+"_noNaNs";
		run("Duplicate...", "title=["+lifetime_stack_for_smoothing+"] duplicate");
		//Create NaNs from zeros and then 'pad' them before smoothing on the GPU
		getDimensions(width, height, channels, slices, frames);
		if(frames>1) {
			for (f = 1; f <= frames; f++) {
		    	Stack.setFrame(f);
				changeValues(0, 0, NaN);
			}
		}
		else changeValues(0, 0, NaN);
		run("Remove NaNs...", "radius="+smoothRadiusOverlayXY+" stack");
		Ext.CLIJ2_push(lifetime_stack_for_smoothing);
		close(lifetime_stack_for_smoothing);
		if(frames>1) Ext.CLIJ2_mean3DSphere(lifetime_stack_for_smoothing, lifetime_stack_filtered, smoothRadiusOverlayXY, smoothRadiusOverlayXY, smoothRadiusOverlayTime);
		else Ext.CLIJ2_mean2DSphere(lifetime_stack_for_smoothing, lifetime_stack_filtered, smoothRadiusOverlayXY, smoothRadiusOverlayXY);
		Ext.CLIJ2_pull(lifetime_stack_filtered);
		run(lut);
		setMinAndMax(minLifetime, maxLifetime);
		Ext.CLIJ2_release(lifetime_stack_for_smoothing);
		Ext.CLIJ2_release(lifetime_stack_filtered);
	}
	else run("Duplicate...", "title=splits duplicate");

	run("Calibration Bar...", "location=["+calibrationBarPosition+"] fill=Black label=White number=5 decimal=1 font=12 zoom=1 overlay");
	Overlay.copy;
	rename("splits");
	run("RGB Color");
	run("Split Channels");
	selectWindow(intensity_image);

	run("Enhance Contrast", "saturated="+RGB_brightness);	//TO DO: another way to reliably set the B&C settings
	run("Conversions...", "scale");
	run("16-bit");
	showStatus("Generating overlay...");
	imageCalculator("Multiply 32-bit stack", "splits (red)", intensity_image);
	rename("Red");
	setMinAndMax(0, 65536*256);
	imageCalculator("Multiply 32-bit stack", "splits (green)", intensity_image);
	rename("Green");
	setMinAndMax(0, 65536*256);
	imageCalculator("Multiply 32-bit stack", "splits (blue)", intensity_image);
	rename("Blue");
	setMinAndMax(0, 65536*256);

	run("Merge Channels...", "c1=Red c2=Green c3=Blue");
	rename(lifetime_stack + " (RGB overlay)");
	Overlay.paste;
	rename("RGB_overlay");
	if(frames>1) Stack.setSlice(1);
	setBatchMode("show");
	run("Properties...", "channels=1 slices=1 frames="+frames+" pixel_width="+pixelWidth+" pixel_height="+pixelHeight+" voxel_depth=1.0000 frame="+frameInterval);

	close("splits (red)");
	close("splits (green)");
	close("splits (blue)");

	return "RGB_overlay";
}


function close_windows(name) {
	windowList = getList("window.titles");
	for(i=0 ; i<windowList.length ; i++) {
		if(matches(windowList[i],".*"+name+".*")) {
			selectWindow(windowList[i]);
			run("Close");
		}
	}
}

//Lower threshold an image at a certain percentile
function percentile_threshold(percentile) {
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	total = 0;
	bin=0;
	while (total < nPixels*percentile) {
		total += histogram[bin];
		bin++;
	} 
	setThreshold(bin-1, max);
}


//Rename table headers (remove "Mean([cell_xxx])")
function renameTableHeaders(table) {
	headings = Table.headings(table);
	headers = split(headings, "\t");
	for(i=0;i<headers.length;i++) {
		newHeader = substring(headers[i+1],5,lengthOf(headers[i+1])-1);
		Table.renameColumn(headers[i+1], newHeader);
	}
	Table.update;
}


//Adds two arrays of equal length element-wise
function addArrays(array1, array2) {
	added_array=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		added_array[a]=array1[a] + array2[a];
	}
	return added_array;
}


//Subtracts two arrays of equal length element-wise
function subtractArrays(array1, array2) {
	subtracted_array=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		subtracted_array[a]=array1[a] - array2[a];
	}
	return subtracted_array;
}


//Returns the differentiated array. The first element is set to zero
function differentiateArray(array) {
	diffArray = newArray(array.length-1);
	for (i = 1; i < array.length-1; i++) {
		diffArray[i] = array[i] - array[i-1];
	}
	diffArray[0] = 0;
	return diffArray;
}


//Adds a scalar to all elements of an array
function addScalarToArray(array, scalar) {
	added_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		added_array[a] = (array[a]) + (scalar);
	}
	return added_array;
}


//Divides the elements of two arrays and returns the new array
function divideArrays(array1, array2) {
	divArray=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		divArray[a]=array1[a]/array2[a];
	}
	return divArray;
}


//Multiplies all elements of an array with a scalar
function multiplyArraywithScalar(array, scalar) {
	multiplied_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		multiplied_array[a]= (array[a]) * (scalar);
	}
	return multiplied_array;
}


//Divides all elements of an array by a scalar
function divideArraybyScalar(array, scalar) {
	divided_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		divided_array[a]=array[a]/scalar;
	}
	return divided_array;
}


//Returns the indices at which a value occurs within an array
function indexOfArray(array, value) {
	count=0;
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			count++;
		}
	}
	if (count>0) {
		indices=newArray(count);
		count=0;
		for (a=0; a<lengthOf(array); a++) {
			if (array[a]==value) {
				indices[count]=a;
				count++;
			}
		}
		return indices;
	}
}
