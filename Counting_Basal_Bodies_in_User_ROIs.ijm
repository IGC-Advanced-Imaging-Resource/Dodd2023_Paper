/*
 * Macro counting basal body per user defined ROI
 * 
 * Expected Input:
 * The expected input for this script are non-timelapse stacks in lif format where the second channel represents the basal bodies that need quantified 
 * At the start of macro, user selects a parent directory where there images can be found (nested folders dealt with) and a directory where they wish the results to be saved
 * 
 * Description of script workflow:
 * Maximum intensity projection is created of image
 * User is directed to draw around cells and save the ROIs to the ROI manager. 
 * Each ROI is used to count basal bodies per cell (using the ImageJ Find Maxima function) 
 * Optional setting to only select images with "LNG" in title, to allow for analysis of files that have been processed using the Leica LIGTNING processing
 * 
 * Output:
 * An RGB image is saved with the user drawn cell outlines flattened
 * An ImageJ ROI set (.zip) of the user drawn cell outlines
 * A composite image with the basal body channel and the second image represents the points of the detected basal bodies
 * A .csv file containing file name, cell, area of cell and number of basal bodies found in cell
 * Sub-directory with x-y co-ordinates of basal bodies per cell
 * 
 * Installation Notes: 
 * Download FIJI here - https://imagej.net/software/fiji/downloads                                                
 * How to use an ImageJ macro from Github - https://github.com/IGC-Advanced-Imaging-Resource/ImageJ_Macros
 * This script is hosted here with an example image - https://github.com/IGC-Advanced-Imaging-Resource/Dodd2023_paper       
 * 
 * Written by Laura Murphy (laura.murphy@ed.ac.uk)                                                                                                                                                                                                                                                                                                                             
 * IGC Advanced Imaging Resource (https://www.ed.ac.uk/institute-genetics-cancer/facilities/advanced-imaging-resource)
 * First written: November February 2023  Last updated July 2023                                                                                                       				
*/

//--------------------------------//-------------------------------------------------------------------------------------------------------
//-- Part 0: Preparation steps: get directories from user, set up arrays to store results that are being added to and some other setup
//--------------------------------//-------------------------------------------------------------------------------------------------------

// -- Get user to select folder with images for input and the output folder where they want results saved
dir = getDirectory("Select a directory containing one or several .LIF files.");
saveloc = getDirectory("Choose the directory for saving");

// -- Create subdirectory for xy coordinates
xysaveloc = saveloc + "X_Y_Coords" + File.separator;
File.makeDirectory(xysaveloc);

// -- Get LIF files in input directory - function at end of macro
dirList = newArray();
dirList = getFileTree(dir, dirList);

// -- Magic numbers for basal body segmentation
maxima_tolerance = 1; // For analysis in ependymal cells, 800 was used
top_hat_radius = 8; // Sigma for top hat filter
median_radius = 1; // Sigma for median filter

// -- Only Lightning images
choose_lng = 1 // Change to 0 to turn off

run("Bio-Formats Macro Extensions");
setOption("ScaleConversions", true);

// Exmpty arrays to save results 
Filename = newArray();
Cell = newArray();
No_Spots = newArray();
ROI_Area = newArray();

// Start looping through images in input folders
for(f=0; f<dirList.length; f++) {
	path = dirList[f];
	Ext.setId(path);
	Ext.getCurrentFile(file);
	Ext.getSeriesCount(seriesCount);
	
	for (s = 0; s < seriesCount; s++) { 
		Ext.setSeries(s);
		Ext.getSeriesName(seriesName);
		seriesName = replace(seriesName, "/", "_");
		
		if (choose_lng == 1){
			if(matches(seriesName, ".*Lng.*")){
				seriesNameLng = seriesName;
			} else {
			}
			t = s+1;
        	run("Bio-Formats Importer", "open=&path autoscale color_mode=Composite view=Hyperstack stack_order=XYCZT series_" + t);
			getDimensions(width, height, channels, slices, frames);
			run("Z Project...", "projection=[Max Intensity]");
			ori = getTitle();
			print("The name of the file being processed is " + ori);
			print("The name results will be saved in is " + seriesNameLng);
			Stack.setChannel(2);
			run("Duplicate...", "title=green_channel");
			run("Duplicate...", "title=measure_green");
			
			// -- User draws ROIs
			selectWindow(ori);
			Stack.setChannel(2);
			run("Enhance Contrast", "saturated=0.35");
			setTool("freehand");
			waitForUser("Draw around your cells, clicking 't' once each one is finished to add to the ROI Manager\n\nIf you want to skip this image just click 'OK' without drawing any cells");
			
			if (roiManager("Count") == 0){
				run("Close All");
				roiManager("Reset");
				continue;
			}
			
			// Run filters on basal body channels
			selectWindow("green_channel");
			run("Top Hat...", "radius=" + top_hat_radius);
			run("Median...", "radius=" + median_radius);
			setAutoThreshold("Otsu dark");
			
			// Loop through ROIs to count maxima per cell - basal bodies
			for (r = 0; r < roiManager("Count"); r++){ 
			    selectWindow("green_channel");
			    roiManager("Select", r);    
			   	getStatistics(area, mean, min, max, std, histogram);
			   	run("Find Maxima...", "prominence=" + maxima_tolerance + " above output=[List]");
				spot_count = nResults();
	 		    roiManager("rename", "Cell_" + r+1);
			    roiManager("Select", r);     
			    saveAs("Results", xysaveloc + seriesNameLng + "_Cell" + r + 1 + ".csv");
				run("Clear Results");				
				
				// Save results for filename, cell, no of spots per cell and area per cell
				Filename = Array.concat(Filename, seriesNameLng);
				Cell = Array.concat(Cell, r+1);
				No_Spots = Array.concat(No_Spots, spot_count);
				ROI_Area = Array.concat(ROI_Area, area);	
							
			}					
			
			// Rerun filters and find maxima on channel to get single points for saving image
			selectWindow("measure_green");
			run("Top Hat...", "radius=" + top_hat_radius);
			run("Median...", "radius=" + median_radius);
			setAutoThreshold("Otsu dark");
		    roiManager("Deselect");
		    roiManager("Combine");
			run("Find Maxima...", "prominence=" + maxima_tolerance + " above output=[Single Points]");
			run("8-bit");
			selectWindow("measure_green");
			run("8-bit");
			
			run("Merge Channels...", "c1=measure_green c2=[measure_green Maxima] create");
			saveAs("Tiff", saveloc + seriesNameLng + "_spots.tiff");
			run("Close");
				
			// Save images with ROIs flattened
			selectWindow(ori);
			run("RGB Color");
			roiManager("Show All with labels");
			run("Flatten");
			
			saveAs("Tiff", saveloc + seriesNameLng + "_CellOverlay.tiff");
			roiManager("Save", saveloc + seriesNameLng + "_RoiSet.zip");
		
			Array.show(Filename, Cell, No_Spots, ROI_Area);
			saveAs("Results", saveloc + "Full_Results.csv"); 
			run("Close");
					
			run("Close All");
			run("Clear Results");
			roiManager("Reset");
		}
	}        
}	

setBatchMode(false);
Dialog.create("Progress");
Dialog.addMessage("Saving Complete!");
Dialog.show;


//--------------------------------//-----------------------------------------------------------------------
//-- Function for finding all files with .lif extension in directory ad subdirectories
//--------------------------------//-----------------------------------------------------------------------

function getFileTree(dir , fileTree){
	list = getFileList(dir);

	for(f = 0; f < list.length; f++){
		if (matches(list[f], "(?i).*\\.(lif)$"))
			fileTree = Array.concat(fileTree, dir + File.separator + list[f]);
		if(File.isDirectory(dir + File.separator + list[f]))
			fileTree = getFileTree(dir + list[f],fileTree);
	}
	return fileTree;
}
