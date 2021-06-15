import os,itertools
import DSM_functions
import numpy as np
import pyvista as pv

# -------- VARIABLE DEFINITONS -------- #
ROI_Name = 'RECTUM'		#what contour to make a DSM for
DCM_folderpath = 'DCMfolder/'		#folder where all dicoms of a given patient are stored
#The below three lines are useful if you have CT,RTStruct,RTDose, and RTPlan Dicoms in the same folder and don't know their exact names
  #If you know the exact file names you can simply provide them instead
Allfiles = os.listdir(DCM_folderpath)
RS_file = next((s for s in Allfiles if 'RS.' in s), None)	#searches and returns first RT Structure DICOM in the directory
RD_file = next((s for s in Allfiles if 'RD.' in s), None)	#searches and returns first RT Dose DICOM in the directory

###### -----------------
## EXECUTION CODE
###### -----------------

# ------ STAGE 1: Getting and constructing the contour mesh
# get the vertices of the contour mesh from the dicom file
	#at the same time, the function will also return a list of the centers of mass for each slice of the contour
xs,ys,zs, COM = DSM_functions.get_contourpoint_cloud('RECTUM', DCM_folderpath+RS_file)
#create a pyvista point cloud object
cloud = pv.PolyData( np.array([xs,ys,zs]).T )
#Apply delaunay triangulation to create a volume from the mesh, then request just the surface shell of the volume
volume = cloud.delaunay_3d(alpha=2.)
shell = volume.extract_geometry()

# ------ STAGE 2: Definte what planes you will sample your contour with to create your DSM
###		FOR A DSM THAT JUST USES THE AXIAL SLICES AS THEY ARE DEFINED IN THE DICOM FILE
	#we can just use the list of CoMs
centerline = np.array(COM)
samppoints = np.arange(len(COM))	#this tells the sampling function to use all slices (alternatively, could request every second slice) 
	#NOTE: If you want to sample slices of the DSM differently than the basic axial Buettner approach please check the non-axial DSM example

#Next we generate a dictionary of points along the centerline to center sampling planes and the normals of these planes
	#this step is technically optional for planar DSMs, as all planes have the normal [0,0,1]
	#if you skip it, replace the second input to the sample_uniangular function with [0,0,1] and third with centerline[x]
Slices = DSM_functions.sample_contour_surface(volume,centerline,samppoints,NonPlanar=False)


# ------ STAGE 3: Retrive the dose grid information from the RT Dose DICOM
	#dosegrid is the 3D dose grid, in Gy (indexes of the array are z=0,y=1,x=2)
	#the _pos variables give the positions (in mm) along the three axes where the 
	#voxel centers lie and are used to find the dose at each point on the DSM
dosegrid, xpos,ypos,zpos = DSM_functions.get_RT_Dose_info(DCM_folderpath,RD_file)

# ------ STAGE 4: Calculate the dose surface map
	#in this example we will sample dose to the rectal wall at 30 points spaced 12 degrees from each other

dp_list = []	#list where the coordinates of all points on the DSM will be stored
for x in samppoints:
    print('--------',x) #gives the slice index being worked on
    dps,fails = DSM_functions.sample_uniangular(shell,Slices[x]['planenormal'],Slices[x]['point'],30) 
    #the fails variable is used to track which points ray tracing failed to find a point on the contour surface for 
    #and the alternate (slower) interpolation method was used instead
    dp_list.append(dps)
    
DSM = DSM_functions.sample_dose(dp_list,dosegrid,xpos,ypos,zpos)	#performs 3D interpolation on the dose grid
DSM = np.squeeze(DSM)	#removes extra dimensions (as sometimes sample_dose returns a 1 slice 3D array rather than the desired 2D)

DSM_functions.plot_single_DSM(DSM,'Planar DSM','DSMpreview')	#DSM to plot, title on plot, save name for the file
