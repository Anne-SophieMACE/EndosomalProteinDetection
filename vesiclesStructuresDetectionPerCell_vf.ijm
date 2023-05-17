/*Tested with ImageJ 1.53c on Windows
 * macro_vesicle_distribution_v7 
 * DR
 * 
 */

setOption("ExpandableArrays",true);
// cleans everything
run("ROI Manager...");
roiManager("reset");
run("Clear Results");
run("Close All");

/////////////////////////////////////////////////////////////
////// begining of parameters customizable by the user //////
/////////////////////////////////////////////////////////////
// minimum area for vesicles, in pixels
minAreaVes_pix = 30;
// minimum area for structures around vesicles, in pixels
minSizePart_EnlargedROI_pix = 10;
/////////////////////////////////////////////////////////////
//////// end of parameters customizable by the user /////////
/////////////////////////////////////////////////////////////

// asks the user to choose an image
file_name = File.openDialog("Choose your file");
run("Bio-Formats", "open=["+file_name+"] color_mode=Composite open_files view=Hyperstack stack_order=XYCZT");

img_tit = getTitle();
dir_img = getDirectory("image");
ext_beg = indexOf(file_name, "."); // to get only the name without extension
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel"); // removes scale
getDimensions(width, height, channels, slices, frames);
// time to check the channel and slice for analysis
waitForUser("Check on your image the best plan(s) for vesicles detection");

// maximum number of channel for ring detection
nb_chan_enlarged = 2;
chan_enlarged_tab = newArray(2);
// Dialog box to choose slices and channels of analysis
choice_ch = newArray(channels);
for (i = 0; i < channels; i++)
	choice_ch[i] = parseInt(i+1); 
choice_ch = append(choice_ch,"None");
Dialog.create("Parameters for analysis")
Dialog.addChoice("Channel for vesicles detection",choice_ch, "4.0");
Dialog.addChoice("1st channel for structure detection (enlarged ROIs of the vesicles)",choice_ch, "3.0");
Dialog.addChoice("2nd channel for structure detection (enlarged ROIs of the vesicles)",choice_ch, "None");
Dialog.addMessage("Please choose the plans for the study, between 1 and "+slices);
Dialog.addMessage("The number of additional plan will be added in each direction");
Dialog.addMessage("and the maximum projection is applied on the stack with these plans");
Dialog.addMessage("(ex: if 0, only the focal plan, if 1, 1 plan above and 1 below added)");
Dialog.addNumber("Number of additional plans",1);
Dialog.addNumber("Focal plan", 1);
Dialog.addNumber("Number of pixels for enlarge of the vesicles ROIs", 7);
Dialog.show();
chan_vesicles = Dialog.getChoice();
chan_enlarged_tab[0] = Dialog.getChoice();
chan_enlarged_tab[1] = Dialog.getChoice();
number_add_plan = Dialog.getNumber();
slice_best_vesic = Dialog.getNumber();
width_enlarged = Dialog.getNumber();

if( chan_enlarged_tab[1] == "None")
	nb_chan_enlarged = 1;

if( chan_enlarged_tab[0] == "None")
	nb_chan_enlarged = 0;

subname_acq = substring(file_name,lastIndexOf(file_name, File.separator)+1,ext_beg);
// checks if there is already a ROI result
if( File.exists(dir_img+File.separator+subname_acq+"_cellROIs.zip" ) ){
	roiManager("open",dir_img+File.separator+subname_acq+"_cellROIs.zip");
	roiManager("Show All");
}

// asks the user to draw (modifiy if loaded) the ROIs
selectWindow(img_tit);
run("Duplicate...", "duplicate channels="+chan_vesicles);
run("Z Project...", "projection=[Max Intensity]");
setTool("freehand");
run("Enhance Contrast", "saturated=0.70");
waitForUser("Draw (or modify) you cells on the image (draw and add to the Manager)");
nbCell = roiManager("count");

// rename the ROIs with "Cell" prefix
for (i = 0; i < nbCell; i++) {
	roiManager("select", i);
	roiManager("rename", "Cell"+i);	
}

// saves ROIs
roiManager("deselect");
roiManager("save",dir_img+File.separator+subname_acq+"_cellROIs.zip" );

// measures ROIs to compute the ellipse fit and keeps centroids + angles 
run("Set Measurements...", "centroid fit stack redirect=None decimal=5");
roiManager("measure");

centroidX = newArray(nbCell);
centroidY = newArray(nbCell);
AngleAxe = newArray(nbCell);

for (index_analyse = 0; index_analyse < nbCell; index_analyse++) {
	centroidX[index_analyse] = getResult("X", index_analyse);
	centroidY[index_analyse] = getResult("Y", index_analyse);
	AngleAxe[index_analyse] = getResult("Angle", index_analyse);

	if( AngleAxe[index_analyse] > 90 )
		AngleAxe[index_analyse] = AngleAxe[index_analyse] -180;
	if( AngleAxe[index_analyse] < -90 )
		AngleAxe[index_analyse] = AngleAxe[index_analyse] + 180;
}

roiManager("show all");
roiManager("show none");
roiManager("reset");

// duplicate of the channel of vesicles -> careful to check if the slices do exist
selectWindow(img_tit);
run("Duplicate...", "title=vesicle_chan duplicate channels="+chan_vesicles+" slices="+maxOf(1,slice_best_vesic-number_add_plan)+"-"+minOf(slices, slice_best_vesic+number_add_plan));

run("Clear Results"); 
IJ.renameResults("Results","Results_vesicles");

updateResults(); // to make an (empty) Results table appear
IJ.renameResults("Results","Results_enlargedROIs");


run("Set Measurements...", "area mean centroid integrated redirect=None decimal=5");
for (index_analyse = 0; index_analyse < nbCell; index_analyse++) {
	// table will contain the number of vesicles in the vesicle channel + in each enlarged channels
	nbVesicles = newArray(nb_chan_enlarged+1);
	
	// rotation + translation of the image: the cell of interest is centered and its major axi at the horizontal
	createRotatedImage("vesicle_chan",index_analyse,width,height,centroidX[index_analyse],centroidY[index_analyse],AngleAxe[index_analyse]);
	
	// translation/rotation of the ROI of the concerned cell
	roiManager("open",dir_img+File.separator+subname_acq+"_cellROIs.zip");
	roiManager("select",index_analyse);
	roiManager("Add");
	roiManager("select",roiManager("count")-1);
	roiManager("translate", width/2-centroidX[index_analyse], height/2-centroidY[index_analyse]);
	roiManager("select",roiManager("count")-1);
	run("Rotate...", "rotate angle="+AngleAxe[index_analyse]);
	roiManager("Update");
	
	// removes the ROI of the cells (loaded)
	for (i = 0; i < nbCell; i++) { 
		roiManager("select", 0);
		roiManager("delete");
	}
	
	roiManager("show all");
	roiManager("show none");
	// z projection, pre-processing and automatic detection
	if( nSlices > 1)
		run("Z Project...", "projection=[Max Intensity]");
	
	zproj_name = getTitle();
	createMask(index_analyse,true);
	name_img_mask = getTitle();

	roiManager("select", roiManager("count")-1);
	run("Analyze Particles...", "size="+minAreaVes_pix+"-Infinity pixel add");

	selectWindow(name_img_mask);
	close();
	
	if( nSlices > 1){
		selectWindow("vesicle_chan_rotTr"+index_analyse);
		close();
	}
	
	// offers the possibility to the user to modify the vesicles ROIs
	selectWindow(zproj_name);
	resetMinAndMax();
	roiManager("Show All without labels");
	waitForUser("Modify the detected vesicles if necessary");
	roiManager("deselect");
	// first ROI is the cell, following are detected vesicles
	for (i_r = 1; i_r < roiManager("count"); i_r++) {
		roiManager("select", i_r);
		roiManager("rename", "Vesicle"+i_r);
	}
	// nb of vesicles in the vesicles channel
	nbVesicles[0] = roiManager("count")-1;

	// Image containing only the (rotated/translated) mask of the concerned cell 
	selectWindow(img_tit);
	run("Duplicate...", "title=mask_cell duplicate channels=1 slices=1");
	selectWindow("mask_cell");
	run("8-bit");
	run("Multiply...", "value=0.00000");
	roiManager("select",0);
	run("Add...", "value=255");
	// enlarge each vesicle: enlarged ROIs added in the Manager
	for (i_v = 1; i_v <= nbVesicles[0]; i_v++) { 
		roiManager("select", i_v);
		run("Enlarge...", "enlarge="+width_enlarged);
		roiManager("add");
	}

	// combination of all enlarged particles
	roiManager("deselect");
	tab_ind_bigROI = newArray(nbVesicles[0]);
	for (i_v = 0; i_v < nbVesicles[0]; i_v++) 
		tab_ind_bigROI[i_v] = nbVesicles[0]+i_v+1;
	
	roiManager("select", tab_ind_bigROI);
	roiManager("Combine");
	roiManager("add");
	roiManager("deselect");
	roiManager("select", tab_ind_bigROI);
	roiManager("delete");

	// rename the ROI so that the user does not remove it
	index_ROI_enlarged = roiManager("count")-1;
	roiManager("select", index_ROI_enlarged);
	roiManager("rename", "ROI_enlarged_DONOTREMOVE");
	roiManager("deselect");

	// for each channel specified as "enlarged" in the interface
	zproj_name_structChan = newArray(nb_chan_enlarged);
	for (i_rg = 0; i_rg < nb_chan_enlarged; i_rg++) {
		nbROIbefore = roiManager("count");
		selectWindow(img_tit);
		run("Duplicate...", "title=chan_structures duplicate channels="+chan_enlarged_tab[i_rg]+" slices="+maxOf(1,slice_best_vesic-number_add_plan)+"-"+minOf(slices, slice_best_vesic+number_add_plan));
		selectWindow("chan_structures");
		
		// rotation of the image
		createRotatedImage("chan_structures",i_rg,width,height,centroidX[index_analyse],centroidY[index_analyse],AngleAxe[index_analyse]);

		// creation of the mask
		if( nSlices > 1)
			run("Z Project...", "projection=[Max Intensity]");
		zproj_name_structChan[i_rg] = getTitle();
		createMask(index_analyse,false);
		name_img_mask = getTitle();
		
		// searches for structures in the enlarged ROIs
		roiManager("select", index_ROI_enlarged);
		run("Analyze Particles...", "size="+minSizePart_EnlargedROI_pix+"-Infinity pixel add");

		selectWindow(zproj_name_structChan[i_rg]);
		roiManager("show all");
		waitForUser("Check the detected particles, modify them if necessary");

		// renames the new ROIs
		for (i = nbROIbefore; i < roiManager("count"); i++) {
			roiManager("select", i);
			roiManager("rename", "Ring"+i_rg+1+"_part"+i-nbROIbefore+1);
		}
		nbVesicles[i_rg+1] = roiManager("count")-nbROIbefore;
		
		selectWindow("chan_structures");
		close();

		selectWindow("img_mask_rot"+index_analyse);
		close();
	}
	// removes the combined ROI (was only for the Analyze particles)
	roiManager("select", nbVesicles[0]+1);
	roiManager("delete");
	
	// measures characteristics of all vesicles
	roiManager("deselect");
	run("Clear Results");
	roiManager("measure");
	
	ves_X = newArray(nResults-1);
	ves_Y = newArray(nResults-1);
	R = newArray(nResults-1);
	theta = newArray(nResults-1);
	Rnormalized = newArray(nResults-1);
	areaVes = newArray(nResults-1);
	meanVes = newArray(nResults-1);
	rawIntVes = newArray(nResults-1);
	ves_assoc = newArray(nResults-1);
	
	for (i = 0; i < nResults-1; i++) { // 1st ROI is the cell
		ves_X[i] = getResult("X", i+1);
		ves_Y[i] = getResult("Y", i+1);
		areaVes[i] = getResult("Area", i+1);
		// distance to the middle = centroid of the cell
		R[i] = sqrt((width/2-ves_X[i])*(width/2-ves_X[i])+(height/2-ves_Y[i])*(height/2-ves_Y[i]));
	}

	findClosestVesicle(ves_X,ves_Y,nbVesicles[0],ves_assoc);

	run("Clear Results");
	// angle between vesicles and cell centroid
	for (i = 0; i < lengthOf(ves_X); i++) 
		theta[i] = - atan2(ves_Y[i]-height/2,ves_X[i]-width/2)*180/PI;

	// draw the line between the center of the cell and the vesicle
	// and looks for intersection with the cell contour
	for (i = 0; i < lengthOf(ves_X); i++) {
		if( width/2 !=ves_X[i] ){ // y=ax+b
			// a= (Y2-Y1)/(X2-X1)
			// (X2,Y2)= center of the cell -> middle of the image
			// (X1,Y1)= vesicle coordinates
			a = (height/2-ves_Y[i])/(width/2-ves_X[i]);
			b = height/2-a*width/2;

			// direction to explore (left or right): depends on the location
			// of the vesicle according to the middle of the cell
			if( ves_X[i]-width/2 > 0)
				sign = 1;
			else 
				sign = -1;
			
			selectWindow("mask_cell");
			// advances on line y=ax+b until meeting the "outside" = 0 
			i_tmp = 0;
			x = width/2+sign*i_tmp*2;
			while( getPixel(x,a*x+b)==255 ) {
				i_tmp++;
				x = width/2+sign*i_tmp*0.5;
			}
			// intersection point
			makePoint(x, a*x+b);
			roiManager("add");
			// line between middle and intersection point
			makeLine(width/2, height/2, x, a*x+b);
			roiManager("add");
			// R normalized is the ratio between the distance center of cell/vesicle and 
			// center of cell/border of the cell
			Rnormalized[i] = R[i]/sqrt((width/2-x)*(width/2-x)+(height/2-(a*x+b))*(height/2-(a*x+b)));
		}
		else { // x = cste
			xconst = width/2;
			
			// direction to explore (right or down): depends on the location
			// of the vesicle according to the middle of the cell
			if( ves_Y[i]-height/2 > 0)
				sign = 1;
			else 
				sign = -1;

			selectWindow("mask_cell");
			
			// advances on line x=cste until meeting the "outside" = 0 
			i_tmp = 0;
			y = height/2+sign*i_tmp*2;
			while( getPixel(xconst,y)==255 ) {
				i_tmp++;
				y = height/2+sign*i_tmp*0.5;
			}
			// intersection point
			makePoint(xconst, y);
			roiManager("add");
			// line between middle and intersection point
			makeLine(width/2, height/2, xconst, y);
			roiManager("add");
			// R normalized is the ratio between the distance center of cell/vesicle and 
			// center of cell/border of the cell
			Rnormalized[i] = R[i]/sqrt((width/2-xconst)*(width/2-xconst)+(height/2-y)*(height/2-y));
		}
	}
	
	run("Clear Results");
	// vesicles results
	selectWindow(zproj_name);
	for (i = 1; i <= nbVesicles[0]; i++) { // ROI 1= studied cell
		roiManager("select", i);
		roiManager("measure");
		meanVes[i-1] = getResult("Mean",i-1);
		rawIntVes[i-1] = getResult("RawIntDen",i-1); 
	}
	
	IJ.renameResults("Results_vesicles","Results");
	for (i = 0; i < nbVesicles[0]; i++) {
		val_rowRes = nResults;
		setResult("Name", val_rowRes, "cell"+index_analyse+1+"_vesicle"+i+1);
		setResult("Channel",val_rowRes,parseInt(chan_vesicles));
		setResult("Area",val_rowRes,areaVes[i]);
		setResult("Mean value",val_rowRes,meanVes[i]);
		setResult("Raw Int Den",val_rowRes,rawIntVes[i]);
		setResult("Distance center", val_rowRes, R[i]);
		setResult("Normalized distance center", val_rowRes, Rnormalized[i]);
		setResult("Angle", val_rowRes,theta[i]);
	}
	IJ.renameResults("Results","Results_vesicles");

	// index end of detected vesicles
	// = index of beginning of structures detected around vesicles
	beg_ind = nbVesicles[0]+1;

	for (i_rg = 0; i_rg < nb_chan_enlarged; i_rg++) {
		run("Clear Results");
		selectWindow(zproj_name_structChan[i_rg]);
		for (i = beg_ind; i < beg_ind+nbVesicles[i_rg+1]; i++) {
			roiManager("select", i);
			roiManager("measure");
			meanVes[i-1] = getResult("Mean",i-beg_ind);
			rawIntVes[i-1] = getResult("RawIntDen",i-beg_ind); 
		}
		
		IJ.renameResults("Results_enlargedROIs","Results");
		for (i = beg_ind-1; i < beg_ind+nbVesicles[i_rg+1]-1; i++) {
			val_rowRes = nResults;
			setResult("Name", val_rowRes, "cell"+index_analyse+1+"_channel"+i_rg+1+"_part"+i-beg_ind+2);
			setResult("Channel",val_rowRes,chan_enlarged_tab[i_rg]);
			setResult("Area",val_rowRes,areaVes[i]);
			setResult("Mean value",val_rowRes,meanVes[i]);
			setResult("Raw Int Den",val_rowRes,rawIntVes[i]);
			setResult("Distance center", val_rowRes, R[i]);
			setResult("Normalized distance center", val_rowRes, Rnormalized[i]);
			setResult("Angle", val_rowRes,theta[i]);
			setResult("Vesicle",val_rowRes,"cell"+index_analyse+1+"_vesicle"+ves_assoc[i]);
		}
		IJ.renameResults("Results","Results_enlargedROIs");
		beg_ind += nbVesicles[i_rg+1];
	}
	
	roiManager("reset");

	// closes all intermediate images
	selectWindow("mask_cell");
	close();
	
	for (i_rg = 0; i_rg < nb_chan_enlarged; i_rg++) {
		selectWindow(zproj_name_structChan[i_rg]);
		close();
	}
	
	selectWindow(zproj_name);
	close();
}

// saves the results by vesicle
IJ.renameResults("Results_vesicles","Results");
saveAs("Results", dir_img+File.separator+"Results_"+subname_acq+"_vesicles_distributions_chan"+parseInt(chan_vesicles)+"_focPlan_"+slice_best_vesic+"add_plan"+number_add_plan+".xls");

// saves the results by structures around the vesicles, if chosen option
IJ.renameResults("Results_enlargedROIs","Results");
if( nb_chan_enlarged > 0 ){
	if( nb_chan_enlarged == 2 )
		name_file = "Results_"+subname_acq+"_enlargedVes_distributions_chans"+parseInt(chan_enlarged_tab[0])+"_"+parseInt(chan_enlarged_tab[1])+"_focPlan_"+slice_best_vesic+"add_plan"+number_add_plan+".xls";
	else
		name_file = "Results_"+subname_acq+"_enlargedVes_distributions_chan"+parseInt(chan_enlarged_tab[0])+"_focPlan_"+slice_best_vesic+"add_plan"+number_add_plan+".xls";
	saveAs("Results", dir_img+File.separator+name_file);
}

run("Close All");

// adds a cell in the table arr, containing value
function append(arr, value){
 	arr2 = newArray(arr.length+1);
 	for (i=0; i<arr.length; i++)
 		arr2[i] = arr[i];
 	arr2[arr.length] = value;
 	return arr2;
}

// this function creates a binary image of the selected image
// depending of the image, 2 different auto-threshold methods
// can be applied
function createMask(ind,ves){
	run("Duplicate...", "title=img_mask_rot"+ind+" duplicate");
	run("Subtract Background...", "rolling=10");
	run("Median...", "radius=2");
	run("8-bit");
	
	if( ves )
		run("Auto Threshold", "method=Otsu white");
	else
		run("Auto Threshold", "method=MaxEntropy white");
}

// rotation and translation of the image so that the cell is in the middle, aligned on the "Major axis":
// translation is computed so that (centrX,centrY) is located in the middle and rotation is of angle_rot
function createRotatedImage(win_name,ind,width,height,centrX,centrY,angle_rot){
	selectWindow(win_name);
	run("Duplicate...", "title="+win_name+"_rotTr"+ind+" duplicate");
	newMiddleX = width/2-centrX;
	newMiddleY = height/2-centrY;
	run("Translate...", "x="+newMiddleX+" y="+newMiddleY+" interpolation=None stack");
	run("Rotate... ", "angle="+angle_rot+" grid=1 interpolation=Bilinear stack");
	
}

// vesiculeiX,vesiculeiY: tables with coordinates of all vesicules; the nbVesicles_1stChan first
// correspond to the 1st detection; vesicle ring contains for each vesicle in the enlarged ROI to which
// original vesicle it corresponds to
function findClosestVesicle(ves_X,ves_Y,nbVesicles_1stChan,vesicle_ring){
	for (i = nbVesicles_1stChan; i < lengthOf(ves_X); i++) {
		closest_ves = -1 ;
		min_dist = 10e5;

		for(j= 0 ; j<nbVesicles_1stChan; j++) {
			if( sqrt( (ves_X[i]-ves_X[j])*(ves_X[i]-ves_X[j])+(ves_Y[i]-ves_Y[j])*(ves_Y[i]-ves_Y[j]) ) < min_dist ){
				min_dist = sqrt( (ves_X[i]-ves_X[j])*(ves_X[i]-ves_X[j])+(ves_Y[i]-ves_Y[j])*(ves_Y[i]-ves_Y[j]) );
				closest_ves = j;
			}
		}
		vesicle_ring[i] = closest_ves+1;
	}
}
