// Ask the user to select the input folder
inputDir = getDirectory("Choose the folder containing the input files");

// Get a list of all TIFF files in the input folder
list = getFileList(inputDir);

// Loop through each file in the folder
for (i = 0; i < list.length; i++) {
    // Check if the file is a TIFF image (.tif or .tiff)
    if (!endsWith(list[i], ".tif") && !endsWith(list[i], ".tiff")) {
        continue; // Skip non-image files
    }
    
    // Open the image file
    open(inputDir + list[i]);
    
    // Extract the filename without extension
    fileName = substring(list[i], 0, lastIndexOf(list[i], "."));
    
    // Ask the user to apply threshold manually
    waitForUser("Image: " + fileName, 
        "Apply the threshold manually and click OK when done.");
    
    // Optionally, these steps can be enabled for further preprocessing:
    // setOption("BlackBackground", true);
    // run("Convert to Mask");
    // run("Fill Holes");
    // run("Watershed");
    
    // Analyze particles and add ROIs to ROI Manager
    run("Analyze Particles...", "size=50-500 display exclude clear add composite");
    
    // Close the processed image
    close();  
    
    // Reopen the original image
    open(inputDir + list[i]);
    roiManager("show all");
    waitForUser;
    run("Measure");

    // Save the ROI set
    roiManager("Save", inputDir + fileName + "_RoiSet.zip");
    
    // Save the measurement results
    saveAs("Results", inputDir + fileName + "_Results.csv");
    
    // Save a copy of the processed image
    saveAs("Tiff", inputDir + fileName + "_processed.tif");
    
    // Close the processed image
    close();
    
    // Clear the results table
    run("Clear Results");
    
    // Reset the ROI Manager
    roiManager("Reset");
}

// Close results table if it's still open
if (isOpen("Results")) {
    close("Results");
}
