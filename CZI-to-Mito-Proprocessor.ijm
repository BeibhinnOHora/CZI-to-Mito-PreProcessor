/*
 * CZI to Mitochondrial Analyzer Preprocessor
 * A Fiji/ImageJ macro for batch processing of confocal microscopy data
 * 
 * Author: Béibhinn O'Hora
 * Affiliation: Virus and Cellular Stress Unit, Institut Pasteur, Paris, France
 * Contact: beibhinn.ohora@pasteur.fr; jean-pierre.vartanian@pasteur.fr
 * DOI: 
 * Date: April 2026
 * Version: 1.0
 *
 * Purpose: Batch process Zeiss CZI microscopy files to extract mitochondrial 
 *          regions of interest (ROIs) for analysis with Mitochondrial Analyzer 
 *          or similar mitochondrial analysis tools (e.g., MiNa).
 * 
 * Workflow:
 * 1. Opens .czi file 
 * 2. Splits channels
 * 3. Creates Max Intensity Projection of the structural channel for ROI assignment
 * 4. Interactive ROI drawing (user draws the cell boundaries)
 * 5. Saves ROI coordinates and morphology metrics (area, perimeter)
 * 6. Extracts corresponding mitochondrial channel regions as 8-bit TIFFs
 * 
 * Outputs (per input file):
 *   - [filename]_ROIs.zip          : ROI coordinates 
 *   - [filename]_metrics.csv   : Area and perimeter measurements
 *   - [filename]_ROI_1.tif, _2.tif : 8-bit mitochondrial stacks by each ROI
 * 
 * Requirements:
 * - ImageJ v1.53c or later (or Fiji) https://fiji.sc
 * - Bio-Formats plugin (bundled with Fiji)
 * 
 * Usage:
 * 1. Set MITO_CHANNEL_PREFIX and ROI_CHANNEL_PREFIX below detailing in which channel is your staining of interest (e.g., "C2-" for mitochondrial marker, "C3-" for structural marker)
 * 2. Run macro and select input (CZI files) and output folders
 * 3. For each image, draw ROIs around cells using ROI Manager (press 'T' to add)
 * 4. Click OK when finished with each image
 * 
 * License: MIT
 */

// ================= USER CONFIGURATION =================

var MITO_CHANNEL_PREFIX = "C2-";   // channel with mitochondrial marker
var ROI_CHANNEL_PREFIX  = "C3-";   // channel with structural marker

// set directories for your input and output files
var INPUT_DIR  = ""; // add input directory
var OUTPUT_DIR = ""; // add output directory

// ======================================================

print("\\Clear");
run("Close All");

// setting up directories
function getDirectoryPath(prompt, defaultPath) {
    if (defaultPath!="" && File.exists(defaultPath))
        return defaultPath;
    return getDirectory(prompt);
}

inputDir  = getDirectoryPath("Select INPUT folder containing .czi files", INPUT_DIR);
outputDir = getDirectoryPath("Select OUTPUT folder for results", OUTPUT_DIR);

if (!File.exists(inputDir) || !File.exists(outputDir))
    exit("ERROR: Input or output directory not found.");

print("Processing folder: "+inputDir);
print("Saving results to: "+outputDir);
print("Mito channel: "+MITO_CHANNEL_PREFIX);
print("ROI channel : "+ROI_CHANNEL_PREFIX);

// sanity check for files
function countFiles(dir, ext) {
    list = getFileList(dir);
    n = 0;
    for (i=0; i<list.length; i++)
        if (indexOf(toLowerCase(list[i]), toLowerCase(ext))>=0)
            n++;
    return n;
}

function findWindow(prefix) {
    titles = getList("image.titles");
    for (i=0; i<titles.length; i++)
        if (indexOf(titles[i], prefix)==0)
            return titles[i];
    return "";
}

function zeroPad(num, width) {
    s = ""+num;
    while (lengthOf(s)<width) s = "0"+s;
    return s;
}

function closeAllImages() {
    titles = getList("image.titles");
    for (i=0;i<titles.length;i++) {
        selectWindow(titles[i]);
        close();
    }
}

// pre-processing function 
fileList = getFileList(inputDir);
processed = 0; errors = 0;

for (i=0; i<fileList.length; i++) {
    name = fileList[i];
    if (!endsWith(name, ".czi")) continue;
    print("\n--- Processing: "+name+" ---");

    base = replace(name, ".czi", "");
    if (File.exists(outputDir+base+"_ROIs.zip")) {
        print("  Output already exists, skipping.");
        continue;
    }
    run("Bio-Formats Importer",
        "open=["+inputDir+name+"] color_mode=Composite stack_order=XYCZT");
    orig = getTitle();
    getDimensions(w,h,c,z,t);
    run("Split Channels");
    close(orig);

    mitoWin = findWindow(MITO_CHANNEL_PREFIX);
    roiWin  = findWindow(ROI_CHANNEL_PREFIX);

    if (mitoWin=="" || roiWin=="") {
        print("  ERROR: channel not found.");
        closeAllImages(); errors++; continue;
    }

    // adding the max projection for ROI drawing
    selectWindow(roiWin);
    if (z>1) {
    run("Z Project...", "projection=[Max Intensity]");
    projTitle = "MAX_" + roiWin;
    selectWindow(projTitle);
    close(roiWin);
    roiWin = projTitle;
}
maxTitle = roiWin;
    run("ROI Manager...");
    roiManager("Reset");
    roiManager("Show All with labels");

    msg = "Draw cell ROIs on: " + maxTitle +
          "\n\n1. Select the Freehand or Polygon selection tool." +
          "\n2. Draw around every cell you wish to include." +
          "\n3. Click **Add** in the ROI Manager after each ROI." +
          "\n4. When finished adding all ROIs, click OK below." +
          "\n5. Click Cancel to skip this image.";

    waitForUser("Define ROIs", msg);

    roiCount = roiManager("count");
    if (roiCount==0) {
        print("  No ROIs were added; skipping image.");
        closeAllImages(); continue;
    }
    print("  "+roiCount+" ROI(s) defined.");

    // save and measure area and perimeter metrics
    roiPath = outputDir+base+"_ROIs.zip";
    roiManager("Save", roiPath);
    run("Set Measurements...", "area perimeter redirect=None decimal=3");
    run("Clear Results");
    for (r=0; r<roiCount; r++) {
        roiManager("select", r);
        roiManager("Rename", "Cell_"+(r+1));
        run("Measure");
    }
    saveAs("Results", outputDir+base+"_metrics.csv");
    run("Clear Results");

    // extract mitochondrial regions to covert to tiff files
    selectWindow(mitoWin);
    for (r=0; r<roiCount; r++) {
        roiManager("select", r);
        run("Duplicate...", "title=TempROI duplicate");
        selectWindow("TempROI");
        run("Clear Outside", "stack");
        run("8-bit");
        saveAs("Tiff", outputDir+base+"_ROI_"+zeroPad(r+1,2)+".tif");
        close();
    }

    closeAllImages();
    if (isOpen("ROI Manager")) {
        selectWindow("ROI Manager"); run("Close");
    }
    print("  Success.");
    processed++;
}

// Recap on macro actions
print("Done: "+processed+" file(s) processed, "+errors+" error(s).");
showMessage("Batch Complete",
            "Processed "+processed+" image(s).\nSee Log window for details.");